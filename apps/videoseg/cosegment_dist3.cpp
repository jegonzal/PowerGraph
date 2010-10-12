#define cimg_display 0
#define cimg_use_jpeg 1
#include "CImg.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <fstream>


#include <graphlab.hpp>
#include <graphlab/schedulers/multiqueue_priority_scheduler.hpp>
#include <graphlab/distributed/graph/distributed_graph.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <boost/program_options.hpp>
#include <graphlab/macros_def.hpp>
using namespace cimg_library;
using namespace graphlab;
const size_t SHARED_FACTOR_ID = 0;
const size_t GAUSSIAN_CLUSTERS = 1;
const size_t GAUSSIAN_CLUSTER_VERSION = 2;
const size_t BOUND_ID = 3;
const size_t DAMPING_ID = 4;

const double FEATURESCALE = 1;
const size_t INTERFRAME_POTENTIAL = 6;
const size_t INTRAFRAME_POTENTIAL = 4;
/*
std::string imagepath = "/mnt/bigbrofs/usr3/ylow/data/soccer/";
std::string inputfile = "partsoccer2.bin";
size_t width = 640/2;
size_t height = 480/2;
size_t numimgs = 107;*/

std::string imagepath = "";
std::string inputfile = "partnsh2.bin";
bool outputimages = false;
size_t width = 480;
size_t height = 640;
size_t numimgs = 480;
std::string labeledbin = "";
std::string labeledimg = "";
double multiplicity = 1.0;
double smooth=5E-3;


size_t tovertexid(size_t img, size_t x, size_t y) {
  return img * width * height + x * height + y;
}

void fromvertexid(size_t vid, size_t &img, size_t &x, size_t &y) {
  img = vid / (width * height);
  vid -= img * width * height;
  x = vid / height;
  vid -= x * height;
  y = vid;
}


struct edge_data {
  graphlab::unary_factor message;
  double edgepot;
  void save(graphlab::oarchive &oarc) const{
    oarc << message;
    oarc << edgepot;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> message;
    iarc >> edgepot;
  }
  
};


struct vertex_data {
  graphlab::unary_factor potential;
  graphlab::unary_factor belief;
  std::vector<float> features;
  size_t numpixels;
  size_t potential_lastversion;
  
  void save(graphlab::oarchive &oarc) const{
    oarc << potential;
    oarc << belief;
    oarc << features;
    oarc << numpixels;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> potential;
    iarc >> belief;
    iarc >> features;
    iarc >> numpixels;
    potential_lastversion = 0;
  }
};




struct gaussian{
  std::vector<double> mean;
  std::vector<double> variance;
  double totalweight;
  
  double loglikelihood(std::vector<float> &feature) {
    if (totalweight == 0) return -100000000000000000;
    double ll = log(totalweight);
    for (size_t i = 0;i < feature.size(); ++i) {
      ll += log(1/sqrt(variance[i])) 
                        -(feature[i] - mean[i]) * (feature[i] - mean[i]) / (2*variance[i]);
    }
    if (ll < -10000) ll = -10000;
    return ll;
  }
  
  
  void count(double weight, const std::vector<float> &features) {
    for (size_t i = 0; i < features.size(); ++i) {
      mean[i] = mean[i] + weight * features[i];
      variance[i] = variance[i] + weight * features[i] * features[i];
    }
    totalweight += weight;
  }
  void normalize_to(double w = 1.0) {
    for (size_t i = 0; i < mean.size(); ++i) {
      mean[i] = w * mean[i] / totalweight;
      variance[i] = w * variance[i] / totalweight;
    }
    totalweight = w;
  }

  void save(graphlab::oarchive &oarc) const{
    oarc << mean;
    oarc << variance;
    oarc << totalweight;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> mean;
    iarc >> variance;
    iarc >> totalweight;
  }
};

typedef std::vector<gaussian> gaussian_cluster_type;

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;


CImg<unsigned char> RGBtoGrayScale2(const CImg<unsigned char> &im) {
  CImg<unsigned char> grayImage(im.width(),im.height(),im.depth(),3,0);
  if (im.spectrum() == 3) {
    cimg_forXYZ(im,X,Y,Z) {
      grayImage(X,Y,Z,0) = (unsigned char)(0.299*im(X,Y,Z,0) + 0.587*im(X,Y,Z,1) + 0.114*im(X,Y,Z,2));
      grayImage(X,Y,Z,1) = grayImage(X,Y,Z,0);
      grayImage(X,Y,Z,2) = grayImage(X,Y,Z,0);
    }
  }
  return grayImage;
}

double compute_agreementstrength(const std::vector<float> &features1, 
                                 const std::vector<float> &features2,
                                 size_t arity,
                                 double agg_scale_factor) {
  float diff = 0;
  float sum1 = 0;
  float sum2 = 0;
  for (size_t i  =0;i < features1.size(); ++i) {
    sum1 += features1[i] * features1[i];
    sum2 += features2[i] * features2[i];
    diff += features1[i] * features2[i];
  }
  return  agg_scale_factor * diff;    
}
graphlab::binary_factor create_edge_factor(size_t arity, double strength) {
  graphlab::binary_factor bf(0, arity,0,arity);
  bf.set_as_agreement(strength);
  return bf;
}

void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {
  //  std::cout << scope.vertex();;
  //  std::getchar();
  assert(shared_data != NULL);

  // Get the shared data
  double bound = shared_data->get_constant(BOUND_ID).as<double>();
  double damping = shared_data->get_constant(DAMPING_ID).as<double>();

  // Grab the state from the scope
  // ---------------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  
  // see if we need to update the potential
  size_t versionnumber = shared_data->get(GAUSSIAN_CLUSTER_VERSION).as<size_t>();
  bool versionchange = false;
  if (versionnumber > v_data.potential_lastversion) {
    gaussian_cluster_type gaussians = shared_data->get(GAUSSIAN_CLUSTERS).as<gaussian_cluster_type>();
    for (size_t i = 0;i < v_data.potential.arity(); ++i) {
      v_data.potential.logP(i) = gaussians[i].loglikelihood(v_data.features);
      //*                                                  (1.0 + 3.0/(1+exp(-double(v_data.numpixels))));
      versionchange = true;
      v_data.potential_lastversion=versionnumber;
    }
    v_data.potential.normalize();
    /*for (size_t i = 0;i < v_data.potential.arity(); ++i) {
      if (v_data.potential.logP(i) < -10) v_data.potential.logP(i) = -10;
    }*/
    v_data.potential.normalize();
    if (scope.vertex() == 1) {
      std::cout << v_data.potential << "\n";
    }
  }
  // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();


  // Compute the belief
  // ---------------------------------------------------------------->
  // Initialize the belief as the value of the factor
  v_data.belief = v_data.potential;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the message
    const edge_data& e_data = scope.const_edge_data(ineid);
    // Notice we now use the old message since neighboring vertices
    // could be changing the new messages
    v_data.belief.times( e_data.message );
  }
  v_data.belief.normalize(); // finally normalize the belief
  
  // Compute outbound messages
  // ---------------------------------------------------------------->
  /*const graphlab::binary_factor edge_factor =
      shared_data->get_constant(SHARED_FACTOR_ID).as<graphlab::binary_factor>();
*/
  // Send outbound messages
  graphlab::unary_factor cavity, tmp_msg;
  for(size_t i = 0; i < in_edges.size(); ++i) {
    // Get the edge ids
    graphlab::edge_id_t outeid = out_edges[i];
    graphlab::edge_id_t ineid = in_edges[i];
    // CLEVER HACK: Here we are expoiting the sorting of the edge ids
    // to do fast O(1) time edge reversal
    assert(scope.target(outeid) == scope.source(ineid));
    // Get the in and out edge data
    const edge_data& in_edge = scope.const_edge_data(ineid);
    edge_data& out_edge = scope.edge_data(outeid);
    
    // Compute cavity
    cavity = v_data.belief;
    cavity.divide(in_edge.message); // Make the cavity a cavity
    cavity.normalize();

    const graphlab::binary_factor edge_factor = create_edge_factor(cavity.arity(),out_edge.edgepot);
    // convolve cavity with the edge factor storing the result in the
    // temporary message

    tmp_msg.resize(out_edge.message.arity());
    tmp_msg.var() = out_edge.message.var();
    tmp_msg.convolve(edge_factor, cavity);
    tmp_msg.normalize();

    // Damp the message
    tmp_msg.damp(out_edge.message, damping);
    
    // Compute message residual
    double residual = tmp_msg.residual(out_edge.message);
    
    // Assign the out message
    out_edge.message = tmp_msg;
    if (scope.vertex() % 100000 == 0) {
      std::cout << scope.vertex() << " " << residual << std::endl;
    }
    if(residual > bound) {
      gl_types::update_task task(scope.target(outeid), bp_update);      
      scheduler.add_task(task, residual);
    }    
    
  }
/*  if (versionnumber < 10) {
    gl_types::update_task task(scope.vertex(), bp_update);      
    scheduler.add_task(task, 10);
  }*/

} // end of BP_update



/*std::vector<float> featuresquash(const std::vector<float> &infeatures, 
                                      size_t numfeatures) {
  std::vector<float> output;
  size_t numfeaturesperdim = size_t(round(pow(double(numfeatures),(1.0/3.0))));
  
  size_t featperbin = 8 / numfeaturesperdim;
  output.resize(numfeaturesperdim*numfeaturesperdim*numfeaturesperdim);
  for (size_t i = 0;i < infeatures.size(); ++i) {
    size_t c = i;
    size_t r = c / (8*8);
    c = c % (8*8);
    size_t g = c / 8;
    size_t b = c % 8;
    
    r /= featperbin;
    g /= featperbin;
    b /= featperbin;
    size_t o = r * numfeaturesperdim * numfeaturesperdim  + g * numfeaturesperdim + b;
    assert(o < output.size());
    output[o] += FEATURESCALE*infeatures[i];
  }
  return output;
}*/
                   
/**
  In this function, we construct the grid graph
*/
template <typename Gtype, typename SDMType>
void create_graph(std::string archivefile,
             Gtype& g,
             SDMType & sdm,
             std::vector<uint32_t> &pix2part,
             size_t arity, size_t* featurearity,
              double agg_scale_factor) {
  std::ifstream fin(archivefile.c_str());
  graphlab::iarchive iarc(fin);
  
  std::vector<std::vector<size_t> > adjlist;
  std::vector<size_t> part2frame;
  std::vector<std::vector<float> > features;
  std::cout << "deserializing pix2part ... " << std::endl;
  iarc >> pix2part;
  std::cout << "deserializing adjacency list ... " << std::endl;
  iarc >> adjlist;

  std::cout << "deserializing features... " << std::endl;
  iarc >> features;
  part2frame = std::vector<size_t>(features.size(), -1);
  std::vector<size_t> numpixperpart;
  numpixperpart.resize(features.size());
  for (size_t i = 0;i < pix2part.size(); ++i) {
    size_t img,x,y;
    numpixperpart[pix2part[i]]++;
    if (part2frame[pix2part[i]] == size_t(-1)) {
      fromvertexid(i, img, x,y); 
      part2frame[pix2part[i]] = img;
    }
    else {
      fromvertexid(i, img, x,y); 
      //assert(img == part2frame[pix2part[i]]);
    }
  }
  
    // create vertices
  vertex_data vtx;
  vtx.potential.resize(arity);
  vtx.belief.resize(arity);
  vtx.potential_lastversion = 0;
  for (size_t i = 0; i < features.size(); ++i){
    for (size_t j = 0;j < arity; ++j) vtx.potential.logP(j) = gl_types::random::rand01() + 1E-10;
    vtx.belief = vtx.potential;
    vtx.belief.normalize();
    vtx.features = features[i];
    vtx.numpixels = numpixperpart[i];
    g.add_vertex(vtx);
  }
  //create edges
  
  edge_data edata;
  //size_t next_edge_id = 10;
  edata.message.resize(arity);
  edata.message.uniform(0);
  edata.message.normalize();
  for (size_t i = 0; i < adjlist.size(); ++i){
    for (size_t j = 0; j < adjlist[i].size(); ++j){
      double effscalefactor = agg_scale_factor;
      if (part2frame[i] != part2frame[adjlist[i][j]]) {
        effscalefactor = INTERFRAME_POTENTIAL;
      }
      edata.edgepot = compute_agreementstrength(g.vertex_data(i).features, 
                                                        g.vertex_data(adjlist[i][j]).features,
                                                        arity,
                                                        effscalefactor);
      g.add_edge(i, adjlist[i][j], edata);
    }
  }

  g.finalize();
  (*featurearity) = features[0].size();
}







void reduce_gaussian_clusters(size_t index,
                              const gl_types::ishared_data& shared_data,
                              gl_types::iscope& scope,
                              graphlab::any& accumulator) {
  gaussian_cluster_type &cursum = accumulator.as<gaussian_cluster_type>();
  
  const vertex_data& vtx = scope.const_vertex_data();
  for (size_t i = 0; i < vtx.belief.arity(); ++i) {
    cursum[i].count(double(vtx.numpixels) * exp(vtx.belief.logP(i)), vtx.features);
  }
}


void increment(size_t index,
              const gl_types::ishared_data& shared_data,
              graphlab::any& current_data,
              const graphlab::any& new_data) {
  current_data.as<size_t>()++;
}

void merge_gaussian_clusters(size_t index,
              const gl_types::ishared_data& shared_data,
              graphlab::any& current_data,
              const graphlab::any& new_data) {
  gaussian_cluster_type &curval = current_data.as<gaussian_cluster_type>();
  const gaussian_cluster_type & newval = new_data.as<gaussian_cluster_type>();
  for (size_t cluster = 0; cluster < curval.size(); ++cluster) {
    curval[cluster].totalweight += newval[cluster].totalweight;
    for (size_t i = 0; i < curval[cluster].mean.size(); ++i) {
      curval[cluster].mean[i] += newval[cluster].mean[i];
      curval[cluster].variance[i] += newval[cluster].variance[i];
    }
  }
  
}

                                       
void apply_gaussian_clusters(size_t index,
                            const gl_types::ishared_data& shared_data,
                            graphlab::any& current_data,
                            const graphlab::any& new_data) {
  
  gaussian_cluster_type &curval = current_data.as<gaussian_cluster_type>();
  const gaussian_cluster_type & newval = new_data.as<gaussian_cluster_type>();

  double allweight = 0;
  for (size_t cluster = 0; cluster < newval.size(); ++cluster) {
    allweight += newval[cluster].totalweight;
  }
    
  for (size_t cluster = 0; cluster < curval.size(); ++cluster) {
    curval[cluster].totalweight = newval[cluster].totalweight / allweight;
    std::cout << "cluster " << cluster << ": " << curval[cluster].totalweight << std::endl;
    if (curval[cluster].totalweight > 1E-5) {
      for (size_t i = 0; i < curval[cluster].mean.size(); ++i) {
        curval[cluster].mean[i] = newval[cluster].mean[i] / newval[cluster].totalweight;
        curval[cluster].variance[i] = smooth + newval[cluster].variance[i] / newval[cluster].totalweight 
                                            - curval[cluster].mean[i]*curval[cluster].mean[i];
        //std::cout << "\t" << curval[cluster].mean[i] << "\t" << curval[cluster].variance[i] << "\n";
      }
    }
    else {
      for (size_t i = 0; i < curval[cluster].mean.size(); ++i) {
        curval[cluster].mean[i] = 100000000000;
        curval[cluster].variance[i] = 1E-100;
      }
      curval[cluster].totalweight = 0;
    }
  }
  const_cast<gl_types::ishared_data&>(shared_data).atomic_apply(GAUSSIAN_CLUSTER_VERSION, increment, graphlab::any());
}


template <typename T1, typename T2>
float l1dist(const std::vector<T1> &a, const std::vector<T2> &b) {
  float ret = 0;
  for (size_t i = 0;i < a.size(); ++i) {
    ret += fabs(a[i] - b[i]);
  }
  return ret;
}




gaussian_cluster_type 
  read_labelled_image(std::string featarc, std::string labelimgfile, double multiplicity) {
    std::cout << "Reading labelled image from " << labelimgfile << std::endl;
  CImg<unsigned char> labels(labelimgfile.c_str());
  std::ifstream fin(featarc.c_str());
  graphlab::iarchive iarc(fin);
  
  std::vector<uint32_t> pix2part;
  std::vector<std::vector<size_t> > adjlist;
  std::vector<std::vector<float> > features;
  iarc >> pix2part;
  iarc >> adjlist;
  iarc >> features;
  // it should be a one frame image
  //labels should be grayscale
  // if not, just grab any one channel.
  if (labels.depth() == 3) labels.channel(0);
  ASSERT_EQ(labels.width() * labels.height(), pix2part.size());
  
  size_t color2gauss[256];
  size_t nextlbl = 0;
  for (size_t i = 0;i < 256; ++i) color2gauss[i] = size_t(-1);
  
  // each class gets a weighted average of the features of the superpixels in it

  // count the labels
  for (int i = 0; i < labels.width(); ++i) {
    for (int j = 0; j < labels.height(); ++j) {
      if (color2gauss[labels(i,j)] == size_t(-1)) {
        color2gauss[labels(i,j)] = nextlbl;
        nextlbl++;
      }
    }
  }
  gaussian_cluster_type ret;
  std::cout << nextlbl << " colors detected" << std::endl;
  ret.resize(nextlbl);
  for (size_t i  =0;i < ret.size(); ++i) {
    ret[i].mean.resize(features[0].size(), 0);
    ret[i].variance.resize(features[0].size(), 0);
  }
  // one more time, now accumulating 
  std::vector<std::vector<size_t> > pixcounts;
  pixcounts.resize(ret.size());
  for (size_t i = 0;i < pixcounts.size(); ++i) {
    pixcounts[i].resize(features.size(), 0);
  }
  for (int i = 0;i < labels.width(); ++i) {
    for (int j = 0;j < labels.height(); ++j) {
      size_t lbl = color2gauss[labels(i,j)];
      //ASSERT_LT(labels.height() * i + j, pix2part.size());
      size_t part = pix2part[labels.height() * i + j];
      pixcounts[lbl][part]++;
      // add the means
      //ASSERT_LT(lbl, ret.size());
      //ASSERT_LT(part, features.size());      
    }
  }
  for (size_t lbl = 0;lbl < pixcounts.size();++lbl) {
    for (size_t part = 0;part < pixcounts[lbl].size();++part) {
      if (pixcounts[lbl][part] > 0) ret[lbl].count(pixcounts[lbl][part], features[part]);
    }
  }
  // renormalize and scale the weights
  for (size_t i  =0;i < ret.size(); ++i) {
    std::cout << "cluster " << i << ": " << ret[i].totalweight << std::endl;
    for (size_t j = 0; j < ret[i].mean.size(); ++j) {
      std::cout << "\t" << ret[i].mean[j] << "\t" << ret[i].variance[j] << "\n";
    }
  }
  return ret;
}

/**
  Here we create the shared data table
*/
template <typename Gtype, typename SDMType>
void create_shared_data(Gtype &g,
                        SDMType & sdm,
                        size_t arity, size_t featurearity, double agreementstrength) {
  
  
  graphlab::binary_factor edge_potential(0, arity,
                                         0, arity);
  edge_potential.set_as_agreement(agreementstrength);
  std::cout << edge_potential;
  sdm.set_constant(SHARED_FACTOR_ID, edge_potential);
  gaussian_cluster_type vg,vg0;
  
  if (labeledimg != "") {
    vg0 = read_labelled_image(labeledbin, labeledimg, multiplicity);
  }
  else {
    vg0.resize(arity);
    for (size_t i = 0; i < arity; ++i) {
      vg0[i].mean.resize(featurearity);
      vg0[i].variance.resize(featurearity);
      // initialize them to be as far out as possible
      for (size_t j = 0; j < featurearity; ++j) {
        vg0[i].mean[j] = 0;
        vg0[i].variance[j] = 0;
      }
      vg0[i].totalweight = 0;
    }
  }
  ASSERT_EQ(vg0.size(), arity);
  for (size_t i = 0;i < arity; ++i) {
    ASSERT_EQ(vg0[i].mean.size(), featurearity);
    ASSERT_EQ(vg0[i].variance.size(), featurearity);
  }
  
  sdm.set_tree_sync(GAUSSIAN_CLUSTERS,     
               reduce_gaussian_clusters,
               apply_gaussian_clusters, 
               merge_gaussian_clusters,
               vg0,
               100);                 


  vg.resize(arity);  
  
  if (labeledimg != "") {
    vg = vg0;
    for (size_t i = 0; i < arity; ++i) {
      vg[i].normalize_to(1.0);
      vg[i].totalweight = 1.0/arity;
      for (size_t j = 0;j < vg[i].variance.size(); ++j) {
        vg[i].variance[j] += smooth;
      }
    }
  }
  else{
    // randomly pick the first vector
    vg[0].mean.resize(featurearity);
    vg[0].variance.resize(featurearity, 0.2);
    
    for (size_t i = 0; i < arity; ++i) {
      vg[i].mean.resize(featurearity);
      vg[i].variance.resize(featurearity, 0.2);
      // initialize the rest to be as far out as possible
      size_t rvec = rand() % g.num_vertices();
      for (size_t j  =0;j < featurearity; ++j) {
        vg[i].mean[j] = g.vertex_data(rvec).features[j];
      }
      vg[i].totalweight = 1.0/arity;
    } 
  }
  
  for (size_t i = 0; i < arity; ++i) {
    for (size_t j  =0;j < featurearity; ++j) {
      std::cout << vg[i].mean[j] << " ";
    }
    std::cout << "\n";
  } 

  sdm.atomic_set(GAUSSIAN_CLUSTERS,vg);
  sdm.create_atomic(GAUSSIAN_CLUSTER_VERSION, size_t(1));
  sdm.set_constant(BOUND_ID, graphlab::any(1E-3));
  sdm.set_constant(DAMPING_ID, graphlab::any(0.3));
}



template <typename GType>
void makepictures(GType &g, std::vector<uint32_t> pix2part, size_t arity) {
  std::vector<std::vector<unsigned char> > colors;
  colors.resize(arity);
  srand(time(NULL));
  CImg<unsigned char> superpixels(width,height,numimgs, 3,0);
  for (size_t p = 0; p < arity; ++p) {
    colors[p].resize(3);
    colors[p][0] = rand() % 255;
    colors[p][1] = rand() % 255;
    colors[p][2] = rand() % 255;
  }
  for (size_t v = 0; v < pix2part.size(); ++v) {
    size_t i,j,k;
    fromvertexid(v, i,j,k);
    assert(i < numimgs);
    assert(j < width);
    assert(k < height);
    size_t color = g.vertex_data(pix2part[v]).belief.max_asg();
    superpixels(j,k,i,0) = colors[color][0];
    superpixels(j,k,i,1) = colors[color][1];
    superpixels(j,k,i,2) = colors[color][2];
  }
  std::cout << "writing files" << std::endl;
  for (size_t i = 0;i < numimgs; ++i) {
    if (imagepath.size() > 0) {
      char cfilename[1024];
      sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),int(i+1));
      std::cout << "Loading " << cfilename << std::endl;
      CImg<unsigned char> temp(cfilename);
      CImg<unsigned char> temp2 = RGBtoGrayScale2(temp);
      CImg<unsigned char> suppixslice = superpixels.get_slice(i);
      //suppixslice.resize_doubleXY();
      char fname[1024];
      sprintf(fname,"part%d.jpg", int(i+1));
      temp2 = temp2 / 2.0 + suppixslice / 2.0;
      temp2.save(fname);
    }
    else{
      CImg<unsigned char> suppixslice = superpixels.get_slice(i);
      char fname[1024];
      sprintf(fname,"part%d.jpg", int(i+1));
      suppixslice.save(fname);
    }
  }
}

struct options {
  size_t ncpus;
  bool prepart;
  size_t arity;
  std::string scope;
  std::string scheduler;
  std::string basefile;
};

bool parse_command_line(int argc, char** argv, options& opts) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
  // Set the program options
  desc.add_options()
    ("ncpus",  boost_po::value<size_t>(&(opts.ncpus))->default_value(2),
     "Number of cpus to use.")
    ("arity",  boost_po::value<size_t>(&(opts.arity))->default_value(4),
     "Number of initial clusters.")
    ("prepartition",  boost_po::value<bool>(&(opts.prepart))->default_value(false),
     "prepartiton")
    ("scope",
     boost_po::value<std::string>(&(opts.scope))->default_value("edge"),
     "Options are {vertex, edge, full}")
    ("scheduler",
     boost_po::value<std::string>(&(opts.scheduler))->default_value("multiqueue_fifo"),
     "Options are multiqueue_fifo/multiqueue_priority")
    ("infile",  boost_po::value<std::string>(&(inputfile))->default_value(""),
     "infile")
    ("width",  boost_po::value<size_t>(&(width))->default_value(0),
     "width")
    ("height",  boost_po::value<size_t>(&(height))->default_value(0),
     "height")
    ("numimgs",  boost_po::value<size_t>(&(numimgs))->default_value(0),
     "numimgs")
    ("labeledbin",  boost_po::value<std::string>(&(labeledbin))->default_value(""),
     "labeledbin")
    ("labeledimg",  boost_po::value<std::string>(&(labeledimg))->default_value(""),
     "labeledimg")
    ("labeledweight",  boost_po::value<double>(&(multiplicity))->default_value(1.0),
     "labeledweight")
    ("imagepath",  boost_po::value<std::string>(&(imagepath))->default_value(""),
     "imagepath")
     ("outputimages",  boost_po::value<bool>(&(outputimages))->default_value(false),
     "outputimages")
     ("basefile",
     boost_po::value<std::string>(&(opts.basefile))->default_value(""),
     "Base input/output graphname");
  // Parse the arguments
  boost_po::variables_map vm;
  boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
  boost_po::notify(vm);
  if(vm.count("help")) {
    std::cout << desc << std::endl;
    return false;
  }
  return true;
} // end of parse command line arguments


int main(int argc,  char *argv[]) {
  srand(time(NULL));
  //sets the logging level of graphlab
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::distributed_control dc(&argc, &argv);
  dc.init_message_processing(4);
  dc.barrier();
  // create the distributed data structures
  
  options opts; 
  bool success = parse_command_line(argc, argv, opts);
  ASSERT_TRUE(success);
  ASSERT_GT(opts.basefile.length(), 0);  
  

  if (opts.prepart) {
    ASSERT_EQ(dc.numprocs(), 1);
    
    if (fopen(graph_type::filename_for_part(opts.basefile, opts.ncpus-1, opts.ncpus), "r") != NULL) {
        std::cout << "Partition already exists, abort! N of partitions=" << opts.ncpus << std::endl;
        return 0;
    }
    
    // prepartition
    // create graph   
    graph<vertex_data, edge_data> g;
    thread_shared_data<graph<vertex_data, edge_data> > shared_data;
    std::vector<uint32_t> pix2part;
    size_t featurearity;
    create_graph(inputfile, g, shared_data, pix2part, opts.arity, &featurearity, INTRAFRAME_POTENTIAL);
    graph_type::partition_graph_tofile(g, opts.ncpus, partition_method::PARTITION_METIS, opts.basefile);
    return 0;
  }

  // actually do stuff

  gl_types::graph graph; dc.barrier();
  graphlab::distributed_shared_data<gl_types::graph> shared_data(dc); dc.barrier();

  graph.load(opts.basefile, dc);
  dc.barrier();
  gl_types::iengine* engine = NULL;
  
  if (opts.scheduler == "multiqueue_fifo") {
    graphlab::distributed_engine<gl_types::graph,
    graphlab::distributed_scheduler_wrapper<gl_types::graph, 
    graphlab::multiqueue_fifo_scheduler<gl_types::graph> > > *tengine = new graphlab::distributed_engine<gl_types::graph,
    graphlab::distributed_scheduler_wrapper<gl_types::graph, 
    graphlab::multiqueue_fifo_scheduler<gl_types::graph> > >(dc, graph, opts.ncpus);
    tengine->set_caching(true);
    engine=tengine;
  }
  else {
    graphlab::distributed_engine<gl_types::graph,
    graphlab::distributed_scheduler_wrapper<gl_types::graph, 
    graphlab::multiqueue_priority_scheduler<gl_types::graph> > >* tengine = new graphlab::distributed_engine<gl_types::graph,
    graphlab::distributed_scheduler_wrapper<gl_types::graph, 
    graphlab::multiqueue_priority_scheduler<gl_types::graph> > >(dc, graph, opts.ncpus);
    tengine->set_caching(true);
    engine=tengine;
  }
  dc.barrier();
  engine->set_shared_data_manager(&shared_data);  
  dc.barrier();
  if (dc.procid() == 0) {
    std::vector<uint32_t> pix2part;
    //size_t featurearity = 0;
    create_shared_data(graph,shared_data, opts.arity, 
                       graph.vertex_data(graph.my_vertices()[0]).features.size(), 
                       INTRAFRAME_POTENTIAL);
  }
  dc.barrier();

  
 
  if (opts.scope == "edge") {
    engine->set_default_scope(graphlab::scope_range::EDGE_CONSISTENCY);
  }
  else if (opts.scope == "vertex") {
    engine->set_default_scope(graphlab::scope_range::VERTEX_CONSISTENCY);
  }
  else if (opts.scope == "full") {
    engine->set_default_scope(graphlab::scope_range::FULL_CONSISTENCY);
  }
  else {
    ASSERT_TRUE(false);
  }


  graphlab::timer ti;
  
  if (dc.procid() == 0) {
    std::vector<graphlab::vertex_id_t> vec;
    for (graphlab::vertex_id_t i = 0; i < graph.num_vertices(); ++i) {
      vec.push_back(i);
    }
    std::random_shuffle(vec.begin(), vec.end());
    engine->get_scheduler().add_tasks(vec, bp_update, 100.0);
    shared_data.trigger_sync(GAUSSIAN_CLUSTERS);
  }
  dc.barrier();
  ti.start();  
  // Starte the engine
  engine->start();
  double runtime = ti.current_time();
  dc.barrier();
  graph.send_vertices_to_proczero();
  if (dc.procid() == 0) {
    // output the time
    std::cout << "Completed in " << runtime << " seconds" << std::endl;

    
    std::vector<uint32_t> pix2part;
    std::ifstream fin(inputfile.c_str());
    graphlab::iarchive iarc(fin);
    std::cout << "deserializing pix2part ... " << std::endl;
    iarc >> pix2part;

    // make pictures
    if (outputimages) makepictures(graph,pix2part, opts.arity);
  }
  dc.barrier();
}

