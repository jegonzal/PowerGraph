#define cimg_display 0
#include "CImg.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <fstream>

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>
using namespace cimg_library;

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

std::string imagepath = "/mnt/bigbrofs/usr3/ylow/backups/Desktop/make3d/videos/t7/";
std::string inputfile = "partnsh2.bin";
size_t width = 480;
size_t height = 640;
size_t numimgs = 480;


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
  graphlab::unary_factor old_message;
  size_t edge_factor_id;
  void save(graphlab::oarchive &oarc) const{
    oarc << message;
    oarc << old_message;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> message;
    iarc >> old_message;
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
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> potential;
    iarc >> belief;
    iarc >> features;
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
    //if (ll < -10) ll = -10;
    return ll;
  }
  
  
  void count(double weight, const std::vector<float> &features) {
    for (size_t i = 0; i < features.size(); ++i) {
      mean[i] = mean[i] + weight * features[i];
      variance[i] = variance[i] + weight * features[i] * features[i];
    }
    totalweight += weight;
  }
  void save(graphlab::oarchive &oarc) const{
    oarc << mean;
    oarc << variance;
  }
  void load(graphlab::iarchive &iarc) {
    iarc >> mean;
    iarc >> variance;
  }
};

typedef std::vector<gaussian> gaussian_cluster_type;
typedef graphlab::graph<vertex_data, edge_data> graph_type;
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

  // Flip the old and new messages to improve safety when using the
  // vertex scope
  foreach(graphlab::edge_id_t ineid, in_edges) {   
    // Get the in and out edge data
    edge_data& in_edge = scope.edge_data(ineid);
    // Since we are about to receive the current message make it the
    // old message
    in_edge.old_message = in_edge.message;
  }

  // Compute the belief
  // ---------------------------------------------------------------->
  // Initialize the belief as the value of the factor
  v_data.belief = v_data.potential;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the message
    const edge_data& e_data = scope.const_edge_data(ineid);
    // Notice we now use the old message since neighboring vertices
    // could be changing the new messages
    v_data.belief.times( e_data.old_message );
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
    cavity.divide(in_edge.old_message); // Make the cavity a cavity
    cavity.normalize();

    const graphlab::binary_factor edge_factor =
      shared_data->get_constant(out_edge.edge_factor_id).as<graphlab::binary_factor>();
    // convolve cavity with the edge factor storing the result in the
    // temporary message

    tmp_msg.resize(out_edge.message.arity());
    tmp_msg.var() = out_edge.message.var();
    tmp_msg.convolve(edge_factor, cavity);
    tmp_msg.normalize();

    // Damp the message
    tmp_msg.damp(out_edge.message, damping);
    
    // Compute message residual
    double residual = tmp_msg.residual(out_edge.old_message);
    
    // Assign the out message
    out_edge.message = tmp_msg;
    if(residual > bound || versionchange == true) {
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
                   


graphlab::binary_factor create_edge_factor(std::vector<float> features1, 
                                 std::vector<float> features2,
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
  diff  = diff/(sqrt(sum1) * sqrt(sum2));
  graphlab::binary_factor bf(0, arity,0,arity);
  bf.set_as_agreement(agg_scale_factor * diff);    
/*  float diff = 0;
  float sum1 = 0, sum2 = 0;
  for (size_t i  =0;i < features1.size(); ++i) {
    diff += fabs(features1[i] - features2[i]);
  }
  graphlab::binary_factor bf(0, arity,0,arity);
  bf.set_as_agreement(agg_scale_factor / (diff + 1E-1));   */
/*  std::cout << diff << "\n";
  std::cout << bf << "\n";
  getchar();*/
  return bf;
}
/**
  In this function, we construct the grid graph
*/
void create_graph(std::string archivefile,
             graph_type& g,
             gl_types::thread_shared_data& sdm,
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
  
  std::map<std::pair<size_t, size_t>, size_t> edge_factors;
  edge_data edata;
  size_t next_edge_id = 10;
  edata.message.resize(arity);
  edata.message.uniform(0);
  edata.message.normalize();
  edata.old_message = edata.message;
  for (size_t i = 0; i < adjlist.size(); ++i){
    for (size_t j = 0; j < adjlist[i].size(); ++j){
      std::pair<size_t, size_t> ed = std::make_pair(std::min(i, adjlist[i][j]),
                                                    std::max(i, adjlist[i][j]));
      if (edge_factors.find(ed) != edge_factors.end()) {
        edata.edge_factor_id = edge_factors[ed];
        g.add_edge(i, adjlist[i][j], edata);
      }
      else {
        double effscalefactor = agg_scale_factor;
        if (part2frame[i] != part2frame[adjlist[i][j]]) {
          effscalefactor = INTERFRAME_POTENTIAL;
        }
        sdm.set_constant(next_edge_id, create_edge_factor(g.vertex_data(i).features, 
                                                          g.vertex_data(adjlist[i][j]).features,
                                                          arity,
                                                          effscalefactor));
        edge_factors[ed] = next_edge_id;
        ++next_edge_id;
        edata.edge_factor_id = edge_factors[ed];
        //edata.edge_factor_id = SHARED_FACTOR_ID;
        g.add_edge(i, adjlist[i][j], edata);
      }
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
        curval[cluster].variance[i] = 1E-3 + newval[cluster].variance[i] / newval[cluster].totalweight 
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

/**
  Here we create the shared data table
*/
void create_shared_data(graph_type &g,
                        gl_types::thread_shared_data& sdm,
                        size_t arity, size_t featurearity, double agreementstrength) {
  
  
  graphlab::binary_factor edge_potential(0, arity,
                                         0, arity);
  edge_potential.set_as_agreement(agreementstrength);
  std::cout << edge_potential;
  sdm.set_constant(SHARED_FACTOR_ID, edge_potential);
  gaussian_cluster_type vg,vg0;
  
  
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

  sdm.set_sync(GAUSSIAN_CLUSTERS,     
               reduce_gaussian_clusters,
               apply_gaussian_clusters, 
               vg0,
               1000);                 


  vg.resize(arity);  
  
  
  // randomly pick the first vector
  vg[0].mean.resize(featurearity);
  vg[0].variance.resize(featurearity, 0.2);
/*
  std::vector<float> randfeat = g.vertex_data(rand() % (g.num_vertices())).features;
  for (size_t j  =0;j < featurearity; ++j) {
    vg[0].mean[j] = randfeat[j];
    std::cout << vg[0].mean[j] << " ";
  }
  vg[0].totalweight = 1.0/arity;
    
  for (size_t i = 1; i <= arity; ++i) {
    vg[i % arity].mean.resize(featurearity);
    vg[i % arity].variance.resize(featurearity, 0.2);
    // initialize the rest to be as far out as possible
    size_t rvec = 0;
    float bestdist = 0;
    for (size_t r = 0;r < g.num_vertices(); ++r) {
      float mindist = 100000000;
      for (size_t j = 0;j < i; ++j) {
        mindist = std::min(l1dist(g.vertex_data(r).features, vg[j].mean), mindist);
      }
      if (mindist > bestdist) {
        bestdist = mindist;
        rvec = r;
      }
    }
    for (size_t j  =0;j < featurearity; ++j) {
      vg[i % arity].mean[j] = g.vertex_data(rvec).features[j];
      std::cout << vg[i % arity].mean[j] << " ";
    }
    std::cout << "\n";
    vg[i % arity].totalweight = 1.0/arity;
  }*/
  
  for (size_t i = 0; i < arity; ++i) {
    vg[i].mean.resize(featurearity);
    vg[i].variance.resize(featurearity, 0.2);
    // initialize the rest to be as far out as possible
    size_t rvec = rand() % g.num_vertices();
    for (size_t j  =0;j < featurearity; ++j) {
      vg[i].mean[j] = g.vertex_data(rvec).features[j];
      std::cout << vg[i].mean[j] << " ";
    }
    std::cout << "\n";
    vg[i].totalweight = 1.0/arity;
  } 

  sdm.atomic_set(GAUSSIAN_CLUSTERS,vg);
  sdm.create_atomic(GAUSSIAN_CLUSTER_VERSION, size_t(1));
  sdm.set_constant(BOUND_ID, graphlab::any(1E-3));
  sdm.set_constant(DAMPING_ID, graphlab::any(0.4));
}




void makepictures(graph_type &g, std::vector<uint32_t> pix2part, size_t arity) {
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
    char cfilename[1024];
    sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),int(i+1));
    std::cout << "Loading " << cfilename << std::endl;
    CImg<unsigned char> temp(cfilename);
    CImg<unsigned char> temp2 = RGBtoGrayScale2(temp);
    CImg<unsigned char> suppixslice = superpixels.get_slice(i);
    //suppixslice.resize_doubleXY();
    char fname[1024];
    sprintf(fname,"nshcut%d.jpg", int(i+1));
    temp2 = temp2 / 2.0 + suppixslice / 2.0;
    temp2.save(fname);
  }
}




int main(int argc,  char *argv[]) {
  srand(time(NULL));
  //sets the logging level of graphlab
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // create the graph and shared data objects
  graph_type graph;
  gl_types::thread_shared_data shared_data;
  std::vector<uint32_t> pix2part;
  size_t arity = 4;
  double agreementstrength = INTRAFRAME_POTENTIAL;
  size_t featurearity;
  create_graph(inputfile, graph, shared_data, pix2part, arity ,&featurearity, agreementstrength);

  create_shared_data(graph,shared_data, arity, featurearity, agreementstrength);


  graphlab::command_line_options opts;
  opts.parse(argc, argv);
  opts.print();
  
  // create the engine
  gl_types::iengine *engine = opts.create_engine(graph);
  assert(engine != NULL);
  // set the shared data object
  engine->set_shared_data_manager(&shared_data);
  //shared_data.sync(graph, GAUSSIAN_CLUSTERS); 
  // since we are using a task scheduler, we need to
  // to create tasks. otherwise the engine will just terminate immediately
  // there are DIM * DIM vertices
  // the GraphLab timer object is convenient for high resolution timing.
  graphlab::timer ti;
  ti.start();
  // begin executing the engine
  // this function will only return when all tasks are done
  engine->set_timeout(120);
  for (size_t iter = 0;iter < 10; ++iter) {
    if (iter > 0) shared_data.trigger_sync(GAUSSIAN_CLUSTERS);
    for (size_t i = 0;i < graph.num_vertices(); ++i) {
      engine->add_task(gl_types::update_task(i, bp_update), 10.0);
    }

    engine->start();
    
  }
  // output the time
  std::cout << "Completed in " << ti.current_time() << " seconds" << std::endl;

  // make pictures
  makepictures(graph,pix2part, arity);
}

