#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <graphlab.hpp>
#include <serialization/serialization_includes.hpp>
#define cimg_display 0
#include "CImg.h"
// Include the macro for the for each operation

using namespace cimg_library;


struct rgb{ 
  rgb(unsigned char r, unsigned char g, unsigned char b, double vweight):r(r),g(g),b(b), vweight(vweight) {}
  unsigned char r;
  unsigned char g;
  unsigned char b;
  float vweight;
};
typedef graphlab::graph<rgb, float> graph_type;
typedef graphlab::types<graph_type> gl_types;


size_t BLOCKSIZE = 16;

std::string imagepath = "/mnt/bigbrofs/usr3/ylow/backups/Desktop/make3d/videos/t7/";
size_t width = 480;
size_t height = 640;
size_t numimgs = 480;
/*std::string imagepath = "/mnt/bigbrofs/usr3/ylow/data/soccer/";
size_t width = 640/2;
size_t height = 480/2;

size_t numimgs = 107; */

std::string outfile = "partnsh.bin";

size_t step = 1;
size_t numparts = numimgs*20;

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

size_t rgb2idx_hist(unsigned char r, unsigned char g, unsigned char b) {
  return size_t(double(r) / 32.0) * 8 * 8 +
        size_t(double(g) / 32.0) * 8 +
        size_t(double(b) / 32.0);
}

size_t rgb2idx_hist(const rgb &pix) {
  return rgb2idx_hist(pix.r, pix.g, pix.b);
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

float tf(float f) {
  return 1 + 200 * exp(-f/10.0);
  //return std::max(100 - f, 1.0);
}

float tf2(float f) {
  return 1 + 200.0/(1+exp(-f));
  //return std::max(100 - f, 1.0);
}

/*void construct_graph(gl_types::graph& g) {
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
        g.add_vertex(rgb(temp(j,k,0,0),temp(j,k,0,1),temp(j,k,0,2), 0.0));
      }
    }
  }
  
  std::cout << g.num_vertices() << " vertices added " << std::endl;
  std::cout << "Convolving" << std::endl;
  CImgList<float> xyz = everything.get_gradient("xyz", 2);
  CImg<float>& vert = xyz(0);
  CImg<float>& horz = xyz(1);
  CImg<float>& frame = xyz(2);
  vert.abs().normalize(0,100);
  horz.abs().normalize(0,100);
  frame.abs().normalize(0,100);
  
  
  std::cout << vert.width() << " " << vert.height() << " " << vert.depth()<< std::endl;
  std::cout << horz.width() << " " << horz.height() << " " << horz.depth()<< std::endl;
  std::cout << frame.width() << " " << frame.height() << " " << frame.depth()<< std::endl;
  for (size_t i = 0;i < numimgs; ++i) {
    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {
        // add the edges in increasing directions only
        if (i < numimgs - 1) {
          //g.add_edge(tovertexid(i,j,k),tovertexid(i+1,j,k), frame(j,k,i));
          //g.add_edge(tovertexid(i,j,k),tovertexid(i+1,j,k), std::max(tf(horz(j,k,i) + vert(j,k,i))/5, float(1.0)));
          g.add_edge(tovertexid(i,j,k),tovertexid(i+1,j,k), 5);
        }
        if (j < width - 1) {
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j+1,k), 1+tf2(horz(j,k,i)));
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j+1,k), 20+tf(frame(j,k,i) + vert(j,k,i)));
          g.add_edge(tovertexid(i,j,k),tovertexid(i,j+1,k), tf(vert(j,k,i)));
        }
        if (k < height - 1) {
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j,k+1), 1+tf2(vert(j,k,i)));
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j,k+1), 20+tf(frame(j,k,i) + horz(j,k,i)));
          g.add_edge(tovertexid(i,j,k),tovertexid(i,j,k+1), tf(horz(j,k,i)));
        }
      }
    } 
  }
 
} // End of construct graph

*/


void construct_graph2(gl_types::graph& g) {
  CImg<float> vert(width,height,numimgs, 1,0);
  CImg<float> horz(width,height,numimgs, 1,0);

  for (size_t i = 0;i < numimgs; ++i) {
    char cfilename[1024];
    sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),int((i+1)*step));
    std::cout << "Loading " << cfilename << std::endl;
    CImg<unsigned char> temp(cfilename);
    temp.resize_halfXY();
    for (size_t c = 0; c < 3; ++c) {
       CImgList<float> xy = temp.get_channel(c).get_gradient("xy", 2);
       xy(0).abs(); xy(1).abs();
       vert.get_shared_plane(i).max(xy(0));
       horz.get_shared_plane(i).max(xy(1));
    }

    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {
        double vw = sqrt(vert(j,k,i)*vert(j,k,i)+horz(j,k,i)*horz(j,k,i));
        g.add_vertex(rgb(temp(j,k,0,0),temp(j,k,0,1),temp(j,k,0,2),1+100*vw));
      }
    }
  }
  
  std::cout << g.num_vertices() << " vertices added " << std::endl;

  std::cout << vert.width() << " " << vert.height() << " " << vert.depth()<< std::endl;
  std::cout << horz.width() << " " << horz.height() << " " << horz.depth()<< std::endl;
  for (size_t i = 0;i < numimgs; ++i) {
    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {
        // add the edges in increasing directions only
        if (i < numimgs - 1) {
          //g.add_edge(tovertexid(i,j,k),tovertexid(i+1,j,k), frame(j,k,i));
          //g.add_edge(tovertexid(i,j,k),tovertexid(i+1,j,k), std::max(tf(horz(j,k,i) + vert(j,k,i))/5, float(1.0)));
          g.add_edge(tovertexid(i,j,k),tovertexid(i+1,j,k), 3);
        }
        if (j < width - 1) {
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j+1,k), 1+tf2(horz(j,k,i)));
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j+1,k), 20+tf(frame(j,k,i) + vert(j,k,i)));
          g.add_edge(tovertexid(i,j,k),tovertexid(i,j+1,k), tf(vert(j,k,i)));
        }
        if (k < height - 1) {
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j,k+1), 1+tf2(vert(j,k,i)));
          //g.add_edge(tovertexid(i,j,k),tovertexid(i,j,k+1), 20+tf(frame(j,k,i) + horz(j,k,i)));
          g.add_edge(tovertexid(i,j,k),tovertexid(i,j,k+1), tf(horz(j,k,i)));
        }
      }
    } 
  }
 
} // End of construct graph


float weight_function(const float &f) {
  return f;
}

float vertex_weight_function(const rgb &vw) {
  return vw.vweight;
}

std::vector<std::vector<size_t> > solve_adjacencies(CImg<uint32_t> &suppix, size_t parts) {
  std::vector<std::vector<size_t> > ret;
  ret.resize(parts);
  for (size_t i = 0;i < parts; ++i) {
    std::set<size_t> nbrs;
    cimg_forXYZ(suppix,X,Y,Z) {
      if (suppix(X,Y,Z,0) == i) {
        if (X>0) nbrs.insert(suppix(X-1,Y,Z,0));
        if (Y>0) nbrs.insert(suppix(X,Y-1,Z,0));
        if (Z>0) nbrs.insert(suppix(X,Y,Z-1,0));
        if (X<suppix.width()-1) nbrs.insert(suppix(X+1,Y,Z,0));
        if (Y<suppix.height()-1) nbrs.insert(suppix(X,Y+1,Z,0));
        if (Z<suppix.depth()-1) nbrs.insert(suppix(X,Y,Z+1,0));
      }
    }
    nbrs.erase(i);
    std::copy(nbrs.begin(), nbrs.end(), std::back_inserter(ret[i]));
  }
  return ret;
}


int main(int argc, char** argv) {
  // create graph
  //gl_types::graph g;
  //construct_graph2(g);
  std::vector<uint32_t> ret_part;
  std::cout << "starting partitioning..." << std::endl;

  //g.metis_edge_weighted_partition(numparts, ret_part, vertex_weight_function, weight_function);
  size_t nw = width / BLOCKSIZE;
  size_t nh = height / BLOCKSIZE;
  ret_part.resize(width*height*numimgs);
  numparts = nw * nh * numimgs;
  for (size_t i = 0;i < numimgs; ++i) {
    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {
        size_t vidx = tovertexid(i,j,k);
        size_t partidx  = i * nw * nh + (j / BLOCKSIZE) * nh + k / BLOCKSIZE;
        assert(partidx < numparts);
        ret_part[vidx] = partidx;
      }
    }
  }
  std::cout << "partitioning complete!" << std::endl;

  // save the partitioning
  std::ofstream fout(outfile.c_str(), std::ofstream::binary);
  graphlab::oarchive oarc(fout);
  oarc << ret_part;



  // generate the superpixels adjacency structure
  std::vector<std::vector<size_t> > adjlist;
  adjlist.resize(numparts);
  for (size_t v = 0;v < numparts; ++v) {
    // it is a grid!
    size_t vidx = v;
    size_t i = vidx / (nw * nh);
    vidx = vidx % (nw * nh);
    size_t j = vidx / nh;
    size_t k = vidx % nh;
    if (i > 0)  adjlist[v].push_back((i-1) * nw * nh + (j) * nh + k);
    if (i < numimgs-1)  adjlist[v].push_back((i+1) * nw * nh + (j) * nh + k);
    if (j > 0)  adjlist[v].push_back((i) * nw * nh + (j-1) * nh + k);
    if (j < nw-1)  adjlist[v].push_back((i) * nw * nh + (j+1) * nh + k);
    if (k > 0)  adjlist[v].push_back((i) * nw * nh + j * nh + k-1);
    if (k < nh-1)  adjlist[v].push_back((i) * nw * nh + j * nh + k+1);
  }
  for (size_t i = 0;i < numparts; ++i) {
    std::cout << i << ": ";
    for (size_t j = 0; j < adjlist[i].size(); ++j) {
      std::cout << adjlist[i][j] << " ";
    }
    std::cout << std::endl;
  }
  // generate the superpixel features
    // histogram is compacted to 8*8*8 bins
  std::vector<std::vector<float> > features;
  features.resize(numparts);
  for (size_t i = 0;i < numparts; ++i) {
    features[i].resize(2*2*2, 0); 
  }

  oarc << adjlist;
  oarc << features;
  // write
  fout.close();
}

