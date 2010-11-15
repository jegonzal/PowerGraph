#define __STDC_LIMIT_MACROS
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <graphlab.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <boost/program_options.hpp>
#define cimg_display 0
#define cimg_use_jpeg 1
#include "CImg.h"

#include <stdint.h>
#include <graphlab/macros_def.hpp>
// Include the macro for the for each operation

using namespace cimg_library;

const size_t TEMPERATURE = 1;

struct vertexdata{ 
  vertexdata(double vweight, uint64_t id)
      :id(id), counter(0) {}
  uint64_t id;
  uint16_t counter;
};
typedef graphlab::graph<vertexdata, float> graph_type;
typedef graphlab::types<graph_type> gl_types;


std::string imagepath = "/mnt/bigbrofs/usr3/ylow/backups/Desktop/make3d/videos/t7/";
size_t width = 480/2;
size_t height = 640/2;
size_t numimgs = 480;
double newidprior = 1.0;
std::string outfile = "partnsh.bin";
double tempnum = 2.0;
double tempbase = 1.5;
size_t iterations = 5;
size_t countermax = 3;
double threshold = 10;
bool diagonals = false;
/*std::string imagepath = "/mnt/bigbrofs/usr3/ylow/data/soccer/";
size_t width = 640/2;
size_t height = 480/2;

size_t numimgs = 107; */


size_t step = 1;
//size_t numparts = numimgs*20;
size_t numparts = numimgs*20;



bool parse_command_line(graphlab::command_line_options& opt, int argc, char** argv) {
  size_t partfactor;
  // Set the program options
  opt.attach_option("imagepath", &imagepath,std::string(""),"");
  opt.attach_option("width", &width,size_t(0),"");
  opt.attach_option("height", &height,size_t(0),"");
  opt.attach_option("numimgs", &numimgs,size_t(0),"");
  opt.attach_option("outfile", &outfile,std::string(""),"");
  opt.attach_option("partfactor", &partfactor,size_t(20),"");
  opt.attach_option("prior", &newidprior,double(1.0),"");
  opt.attach_option("tempbase", &tempbase,double(1.5),"");
  opt.attach_option("tempnum", &tempnum,double(1.0),"");
  opt.attach_option("iterations", &iterations,size_t(5),"");
  opt.attach_option("countermax", &countermax,size_t(3),"");
  opt.attach_option("threshold", &threshold,double(10),"");
  opt.attach_option("diagonals", &diagonals,bool(false),"");
  opt.parse(argc, argv);
  opt.print();
  width/=2;
  height/=2;
  numparts = numimgs * partfactor;
  std::cout << "w: " << width << std::endl;
  std::cout << "h: " << height << std::endl;
  std::cout << "n: " << numimgs << std::endl;
  std::cout << "outfile: " << outfile << std::endl;
  std::cout << "nparts: " << numparts << std::endl;
  return true;
} // end of parse command line arguments



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


CImg<unsigned char> RGBtoGrayScale(const CImg<unsigned char> &im) {
  CImg<unsigned char> grayImage(im.width(),im.height(),im.depth(),1,0);
  if (im.spectrum() == 3) {
    cimg_forXYZ(im,X,Y,Z) {
      grayImage(X,Y,Z,0) = (unsigned char)(0.299*im(X,Y,Z,0) + 0.587*im(X,Y,Z,1) + 0.114*im(X,Y,Z,2));
    }
  }
  return grayImage;
}

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

void proj(float f){ 
  if (f >= threshold) f = 100.0;
  if (f < threshold) f = 100 * (f / threshold);
  f = 100.0 - f;
}


void tf(float &f) {
  f = -log(1 + 200 * exp(-f/10.0));
  //return std::max(100 - f, 1.0);
}

void construct_graph(gl_types::graph& g) {
  CImg<float> everything(width,height,numimgs, 1);

  for (size_t i = 0;i < numimgs; ++i) {
    char cfilename[1024];
    sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),(i+1)*step);
    std::cout << "Loading " << cfilename << std::endl;
    CImg<unsigned char> temp(cfilename);
    temp.resize_halfXY();
    CImg<unsigned char> temp2 = RGBtoGrayScale(temp);
//    temp2.blur(5);
    everything.get_shared_plane(i) = temp2;

    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {
        uint64_t uint64id = gl_types::random::rand_int(UINT64_MAX - 1);;
        //uint64_t uint64id = 0;
        g.add_vertex(vertexdata(0.0, uint64id));
      }
    }
  }
  
  std::cout << g.num_vertices() << " vertices added " << std::endl;
  std::cout << "Convolving" << std::endl;
  everything.blur_bilateral(5,5,5,5,8,8,8,8);
    CImgList<float> xyz = everything.get_gradient("xyz", 2);
  
  CImg<float>& vert = xyz(0);
  CImg<float>& horz = xyz(1);
  CImg<float>& frame = xyz(2);
  vert.sharpen(10);
  horz.sharpen(10);
  frame.sharpen(10);
  vert.abs().normalize(0,100);
  horz.abs().normalize(0,100);
  frame.abs().normalize(0,100);
  for (size_t i = 0;i < numimgs; ++i) {
    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {
        tf(vert(j,k,i));
        tf(horz(j,k,i));
        tf(frame(j,k,i));
        //frame(j,k,i) = frame(j,k,i) / 2;
      }
    }
  }
  vert.save("vert.jpg");
  horz.save("horz.jpg");
  frame.save("frame.jpg");
  

  std::cout << vert.width() << " " << vert.height() << " " << vert.depth()<< std::endl;
  std::cout << horz.width() << " " << horz.height() << " " << horz.depth()<< std::endl;
  std::cout << frame.width() << " " << frame.height() << " " << frame.depth()<< std::endl;
  if (diagonals) {
    for (size_t i = 0;i < numimgs; ++i) {
      if (i % 10 == 0) std::cout << i << std::endl;
      for (size_t j = 0;j < width; ++j) {
        for (size_t k = 0;k < height; ++k) {
          
          for (size_t ii = (i>0?i-1:0) ; ii <= (i<numimgs-1?i+1:i); ++ii) {
            for (size_t jj = (j>0?j-1:0) ; jj <= (j<width-1?j+1:j); ++jj) {
              for (size_t kk = (k>0?k-1:0) ; kk <= (k<height-1?k+1:k); ++kk) {
                if (i==ii && j==jj && k==kk) continue;
                else if (tovertexid(i,j,k) > tovertexid(ii,jj,kk)){
                  float d = 0;
                  //if (i != ii) d += frame(j,k,i) * frame(j,k,i) ;
                  if (i != ii) d += 3;
                  if (j != jj) d += vert(j,k,i) * vert(j,k,i) ;
                  if (k != kk) d += horz(j,k,i) * horz(j,k,i) ;
                  g.add_edge(tovertexid(i,j,k),tovertexid(ii,jj,kk), 2 * sqrt(d));
                }
              }
            }
          }
        }
      }
    } 
  }
  else {
    for (size_t i = 0;i < numimgs; ++i) {
      for (size_t j = 0;j < width; ++j) {
        for (size_t k = 0;k < height; ++k) {
          
          for (size_t ii = (i>0?i-1:0) ; ii <= (i<numimgs-1?i+1:i); ++ii) {
            for (size_t jj = (j>0?j-1:0) ; jj <= (j<width-1?j+1:j); ++jj) {
              for (size_t kk = (k>0?k-1:0) ; kk <= (k<height-1?k+1:k); ++kk) {
                if ((i==ii) + (j==jj) + (k==kk) == 2 
                  && (tovertexid(i,j,k) > tovertexid(ii,jj,kk))){
                  float d = 0;
                  //if (i != ii) d += frame(j,k,i) * frame(j,k,i) ;
                  if (i != ii) d += 3;
                  if (j != jj) d += vert(j,k,i) * vert(j,k,i) ;
                  if (k != kk) d += horz(j,k,i) * horz(j,k,i) ;
                  g.add_edge(tovertexid(i,j,k),tovertexid(ii,jj,kk), sqrt(d));
                }
              }
            }
          }
        }
      }
    } 
  }
  g.finalize();
} // End of construct graph


// energy method
void cluster_update(gl_types::iscope& scope, 
                    gl_types::icallback& scheduler,
                    gl_types::ishared_data* shared_data) {
  
  vertexdata &curvertexdata = scope.vertex_data();
  if (curvertexdata.counter >= countermax) return;
  
  uint64_t oldid = curvertexdata.id;
    // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
  double temperature = shared_data->get_constant(TEMPERATURE).as<double>();
  // grab 
  std::map<uint64_t, double> assignmentpreferences;
  std::set<uint64_t> possibleasgs;
  foreach(graphlab::edge_id_t eid, in_edges) {   
    graphlab::vertex_id_t nbrv = scope.source(eid);
    const vertexdata &nbrdat = scope.const_neighbor_vertex_data(nbrv);
    assignmentpreferences[nbrdat.id] = 0;
    possibleasgs.insert(nbrdat.id);
  }
  
  foreach(graphlab::edge_id_t eid, out_edges) {
    graphlab::vertex_id_t nbrv = scope.target(eid);
    const vertexdata &nbrdat = scope.const_neighbor_vertex_data(nbrv);
    assignmentpreferences[nbrdat.id] = 0;
    possibleasgs.insert(nbrdat.id);
  }
  uint64_t newid;
  if (assignmentpreferences.find(oldid) == assignmentpreferences.end()) {
    newid = oldid;
  }else {
    newid = gl_types::random::rand_int(UINT64_MAX - 1);
  }
  possibleasgs.insert(newid);
  assignmentpreferences[newid] = 0;  
  
  
  // the cost of assigning to a particular color
  // is the sum over all the edge weights which are of a different color.
  // To compute that, we first sum over all the edge weights which are of the
  // same color, and subtract
  
  double totalmass = 0;
  foreach(uint64_t asgid, possibleasgs) {
    foreach(graphlab::edge_id_t eid, in_edges) {   
      graphlab::vertex_id_t nbrv = scope.source(eid);
      const vertexdata &nbrdat = scope.const_neighbor_vertex_data(nbrv);
      if (nbrdat.id == asgid) {
        float edgeweight = scope.const_edge_data(eid);
        assignmentpreferences[nbrdat.id] += -edgeweight / temperature;
        totalmass += -edgeweight / temperature;
      }
      else {
        assignmentpreferences[nbrdat.id] += -2/ temperature;
        totalmass += -2 / temperature;
      }
    }
    
    foreach(graphlab::edge_id_t eid, out_edges) {
      graphlab::vertex_id_t nbrv = scope.target(eid);
      const vertexdata &nbrdat = scope.const_neighbor_vertex_data(nbrv);
      if (nbrdat.id == asgid) {
        float edgeweight = scope.const_edge_data(eid);
        assignmentpreferences[nbrdat.id] += -edgeweight / temperature;
        totalmass += -edgeweight / temperature;
      }
      else {
        assignmentpreferences[nbrdat.id] += -2/ temperature;
        totalmass += -2 / temperature;
      }
    }
  }
  // if oldid now has mass
  // also consider picking a totally new color
//  uint64_t newid;
  
  //totalmass = totalmass;
  // invert so that we have the sum over edge weights of a different color
  std::map<uint64_t, double>::iterator i = assignmentpreferences.begin();
  while(i != assignmentpreferences.end()) {
    i->second = totalmass - i->second;
    ++i;
  }
  //assignmentpreferences[newid] += std::log(numparts - assignmentpreferences.size() + 1);
  
  // normalize by subtract the maximum element
  {
    double m = -1E100;
    i = assignmentpreferences.begin();
    while(i != assignmentpreferences.end()) {
      if (i->second > m) m = i->second;
      ++i;
    }
    // reset the normalizer and sum again
    totalmass = 0;
    i = assignmentpreferences.begin();
    while(i != assignmentpreferences.end()) {
      i->second -= m;
      i->second = std::exp(i->second);
      //std::cout << i->second << " ";
      totalmass += i->second;
      ++i;
    }
    //std::cout << std::endl;
  }
  
  // select the new value
  double r = gl_types::random::rand01() * totalmass;
  
  i = assignmentpreferences.begin();
  double c = 0;
  while(i != assignmentpreferences.end()) {
    c += i->second;
    if (c >= r) {
      curvertexdata.id = i->first;
      break;
    }
    ++i;
  }

  curvertexdata.counter++;
  
  if (curvertexdata.id != oldid) {
    foreach(graphlab::edge_id_t eid, out_edges) {
      scheduler.add_task(gl_types::update_task(scope.target(eid), cluster_update), 10.0); 
    }
    foreach(graphlab::edge_id_t eid, in_edges) {
      scheduler.add_task(gl_types::update_task(scope.source(eid), cluster_update), 10.0); 
    }
  }

}





size_t uf_find(std::vector<size_t> &unionfind, size_t a) {
  if (unionfind[a] == a) return a;
  else {
    unionfind[a] = uf_find(unionfind, unionfind[a]);
    return unionfind[a];
  }
}
void uf_union(std::vector<size_t> &unionfind, size_t a, size_t b) {
  size_t aroot = uf_find(unionfind, a);
  uf_find(unionfind, b);
  unionfind[b] = aroot;
}
void renumber_graph(gl_types::graph &graph) {
  std::vector<size_t> unionfind(graph.num_vertices());
  
  for (size_t v = 0; v < graph.num_vertices(); ++v) {
    unionfind[v] = v;
  }
  for (size_t v = 0; v < graph.num_vertices(); ++v) {
    uint64_t vpartid = graph.vertex_data(v).id;
    graphlab::edge_list in_edges = graph.in_edge_ids(v);
    graphlab::edge_list out_edges = graph.out_edge_ids(v);
    foreach(graphlab::edge_id_t ineid, in_edges) {   
      uint64_t nbrid = graph.source(ineid);
      if (vpartid == graph.vertex_data(nbrid).id) {
        uf_union(unionfind, v, nbrid);
      }
    }
  }
  for (size_t v = 0; v < graph.num_vertices(); ++v) {
    graph.vertex_data(v).id = uf_find(unionfind, v);
  }

}

int main(int argc, char** argv) {
  // create graph
  graphlab::command_line_options opts;
  
  if (parse_command_line(opts, argc,argv) == false) return 0;
  gl_types::graph graph;
  construct_graph(graph);
  std::cout << "starting partitioning..." << std::endl;

  
  gl_types::thread_shared_data shared_data;
  shared_data.set_constant(TEMPERATURE, double(1.0));
  
  // create a permutation of vertices so that we don't update in a 
  // fixed swwp
  std::vector<size_t> perm;
  for (size_t i = 0;i < graph.num_vertices(); ++i) perm.push_back(i);
  std::random_shuffle(perm.begin(), perm.end());
    
  
  // create the engine
  gl_types::iengine *engine = opts.create_engine(graph);
  
  assert(engine != NULL);
  // set the shared data object
  engine->set_shared_data_manager(&shared_data);
  graphlab::timer ti;
  ti.start();
  for (size_t i  =0;i < iterations; ++i) {
    if (i > 0) renumber_graph(graph);
    double temp = double(tempnum / std::pow(tempbase,double(i)));
    //newidprior = newidprior * 2;
    shared_data.set_constant(TEMPERATURE, temp);
    std::set<uint64_t> numparts;
    
    for (size_t i = 0;i < perm.size(); ++i) {
      graph.vertex_data(perm[i]).counter = 0;
      numparts.insert(graph.vertex_data(perm[i]).id);
      engine->get_scheduler().add_task(gl_types::update_task(perm[i], cluster_update), 10.0);
    }
    std::cout << "Starting to Anneal at temperature: " << temp << "\n";
    std::cout << "Currently has " << numparts.size() << " partitions" << std::endl;
    engine->start();
  }
  
    
  renumber_graph(graph);
  
  
  std::cout << "partitioning complete!" << std::endl;
  CImg<unsigned char> superpixels(width,height,numimgs, 3,0);
  CImg<uint32_t> superpixels_idx(width,height,numimgs,1,0);
  std::map<uint64_t, uint32_t> renumber;
  

  // renumber from id to a sequential number
  for (size_t v = 0; v < graph.num_vertices(); ++v) {
    if (renumber.find(graph.vertex_data(v).id) == renumber.end()) {
      size_t d = renumber.size();
      renumber[graph.vertex_data(v).id] = d;
    }
  }
  numparts = renumber.size();
  std::cout << "partitioned into " << numparts << std::endl;
  
  std::vector<std::vector<unsigned char> > colors;
  colors.resize(numparts);
  for (size_t p = 0; p < numparts; ++p) {
    colors[p].resize(3);
    colors[p][0] = rand() % 255;
    colors[p][1] = rand() % 255;
    colors[p][2] = rand() % 255;
  }


  for (size_t v = 0; v < graph.num_vertices(); ++v) {
    uint32_t partid = renumber[graph.vertex_data(v).id];
    size_t i,j,k;
    fromvertexid(v, i,j,k);
    assert(i < numimgs);
    assert(j < width);
    assert(k < height);
    superpixels_idx(j,k,i,0) = partid;
    superpixels(j,k,i,0) = colors[partid][0];
    superpixels(j,k,i,1) = colors[partid][1];
    superpixels(j,k,i,2) = colors[partid][2];
  }
  std::cout << "making adjacency matrix." << std::endl;
  width *= 2;
  height *= 2;
  superpixels_idx.resize(width, height);
  CImg<uint32_t> superpixels_newidx(width,height,numimgs,1,0);
  superpixels.normalize(0,255);
  
  
  
  std::cout << "writing files" << std::endl;
  for (size_t i = 0;i < numimgs; ++i) {
    
    char cfilename[1024];
    sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),(i+1)*step);
    std::cout << "Loading " << cfilename << std::endl;
    CImg<unsigned char> temp(cfilename);
    CImg<unsigned char> temp2 = RGBtoGrayScale2(temp);
    CImg<unsigned char> suppixslice = superpixels.get_slice(i);
    suppixslice.resize_doubleXY();
    char fname[1024];
    sprintf(fname,"scut%d.jpg", int(i+1));
    temp2 = temp2 / 2.0 + suppixslice / 2.0;
    temp2.save(fname);
  }


  // renumber partitions
  std::map<size_t, size_t> newpartid2oldpartid;
  std::map<size_t, size_t> newpartid2frame;
  // frame, oldpartid pair to newpartid
  std::map<std::pair<size_t, size_t>, size_t > partitions;
  size_t nextpartid = 0;
  
  std::vector<uint32_t> newpart(width*height*numimgs, 0);
  for (size_t v = 0; v < width*height*numimgs; ++v) {
    size_t i,j,k;
    fromvertexid(v, i,j,k);
    std::pair<size_t, size_t> idx = std::make_pair(i, superpixels_idx(j,k,i));
    if (partitions.find(idx) != partitions.end()) {
      newpart[v] = partitions[idx];
    }
    else {
      partitions[idx] = nextpartid;
      newpartid2frame[nextpartid] = i;
      newpartid2oldpartid[nextpartid] = superpixels_idx(j,k,i);
      newpart[v] = nextpartid;
      ++nextpartid;
    }
    superpixels_newidx(j,k,i,0) = newpart[v];
  }
  // save the partitioning
  std::ofstream fout(outfile.c_str(), std::ofstream::binary);
  graphlab::oarchive oarc(fout);
  oarc << newpart;

  numparts  = nextpartid;



  // generate the superpixels adjacency structure
  std::vector<std::vector<size_t> > adjlist;
  adjlist.resize(numparts);
  
  std::vector<std::set<size_t> > allnbrs;
  allnbrs.resize(numparts);
  for(size_t f = 0;f < numimgs; ++f) {
    CImg<uint32_t> prevfrm;
    if (f > 0) prevfrm = superpixels_newidx.get_shared_plane(f-1);
    CImg<uint32_t> frm = superpixels_newidx.get_shared_plane(f);
    CImg<uint32_t> nextfrm;
    if (f <numimgs-1) nextfrm= superpixels_newidx.get_shared_plane(f+1);
    std::set<size_t> nbrs;
    cimg_forXY(frm,X,Y) {
      size_t i = frm(X,Y);
      if (X>0) allnbrs[i].insert(frm(X-1,Y));
      if (Y>0) allnbrs[i].insert(frm(X,Y-1));
      if (X<superpixels_newidx.width()-1) allnbrs[i].insert(frm(X+1,Y));
      if (Y<superpixels_newidx.height()-1) allnbrs[i].insert(frm(X,Y+1));
      if (f>0) allnbrs[i].insert(prevfrm(X,Y));
      if (f<numimgs-1) allnbrs[i].insert(nextfrm(X,Y));
    }
  }
  
  for (size_t i = 0; i < numparts;++i) {
    size_t frame = newpartid2frame[i];
    CImg<uint32_t> frm = superpixels_newidx.get_shared_plane(frame);
    std::set<size_t> nbrs = allnbrs[i];
    nbrs.erase(i);
    // add the cross frames
  /*  std::pair<size_t, size_t> pr = std::make_pair(frame + 1, newpartid2oldpartid[i]);
    if (partitions.find(pr) != partitions.end()) {
      nbrs.insert(partitions[pr]);
    }
    pr = std::make_pair(frame - 1, newpartid2oldpartid[i]);
    if (partitions.find(pr) != partitions.end()) {
      nbrs.insert(partitions[pr]);
    }*/
    std::copy(nbrs.begin(), nbrs.end(), std::back_inserter(adjlist[i]));
  }
/*
  for (size_t i = 0;i < numparts; ++i) {
    std::cout << i << ": ";
    for (size_t j = 0; j < adjlist[i].size(); ++j) {
      std::cout << adjlist[i][j] << " ";
    }
    std::cout << std::endl;
  }*/
  // generate the superpixel features
    // histogram is compacted to 8*8*8 bins
  std::vector<std::vector<float> > features;
  features.resize(numparts);
  for (size_t i = 0;i < numparts; ++i) {
    features[i].resize(1, 0); 
  }

  oarc << adjlist;
  oarc << features;
  // write
  fout.close();
}

