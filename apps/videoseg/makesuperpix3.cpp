#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <graphlab.hpp>
#include <serialization/serialization_includes.hpp>
#include <boost/program_options.hpp>
#define cimg_display 0
#define cimg_use_jpeg 1
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


std::string imagepath = "/mnt/bigbrofs/usr3/ylow/backups/Desktop/make3d/videos/t7/";
size_t width = 480/2;
size_t height = 640/2;
size_t numimgs = 480;
std::string outfile = "partnsh.bin";
/*std::string imagepath = "/mnt/bigbrofs/usr3/ylow/data/soccer/";
size_t width = 640/2;
size_t height = 480/2;

size_t numimgs = 107; */


size_t step = 1;
//size_t numparts = numimgs*20;
size_t numparts = numimgs*20;



bool parse_command_line(int argc, char** argv) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
    
  size_t partfactor;
  // Set the program options
  desc.add_options()
    ("help",   "produce this help message")
    ("imagepath",  boost_po::value<std::string>(&(imagepath))->default_value(""),
     "imagepath")
    ("width",  boost_po::value<size_t>(&(width))->default_value(0),
     "width")
    ("height",  boost_po::value<size_t>(&(height))->default_value(0),
     "height")
    ("numimgs",  boost_po::value<size_t>(&(numimgs))->default_value(0),
     "numimgs")
    ("outfile",  boost_po::value<std::string>(&(outfile))->default_value(""),
     "outfile")
    ("partfactor",  boost_po::value<size_t>(&(partfactor))->default_value(20),
     "partfactor") ;

// Parse the arguments
  boost_po::variables_map vm;
  boost_po::store(boost_po::parse_command_line(argc, argv, desc), vm);
  boost_po::notify(vm);
  if(vm.count("help")) {
    std::cout << desc << std::endl;
    return false;
  }
  
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
    sprintf(cfilename,"%s%d.jpg", imagepath.c_str(), int((i+1)*step));
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
  if (parse_command_line(argc,argv) == false) return 0;
  gl_types::graph g;
  construct_graph2(g);
  std::vector<uint32_t> ret_part;
  std::cout << "starting partitioning..." << std::endl;

  g.metis_weighted_partition(numparts, ret_part, 
                             vertex_weight_function, 
                             weight_function);
  std::cout << "partitioning complete!" << std::endl;
  CImg<unsigned char> superpixels(width,height,numimgs, 3,0);
  CImg<uint32_t> superpixels_idx(width,height,numimgs,1,0);

  std::vector<std::vector<unsigned char> > colors;
  colors.resize(numparts);
  for (size_t p = 0; p < numparts; ++p) {
    colors[p].resize(3);
    colors[p][0] = rand() % 255;
    colors[p][1] = rand() % 255;
    colors[p][2] = rand() % 255;
  }
  for (size_t v = 0; v < ret_part.size(); ++v) {
    size_t i,j,k;
    fromvertexid(v, i,j,k);
    assert(i < numimgs);
    assert(j < width);
    assert(k < height);
    superpixels_idx(j,k,i,0) = ret_part[v];
    superpixels(j,k,i,0) = colors[ret_part[v]][0];
    superpixels(j,k,i,1) = colors[ret_part[v]][1];
    superpixels(j,k,i,2) = colors[ret_part[v]][2];
  }
  width *= 2;
  height *= 2;
  superpixels_idx.resize(width, height);
  CImg<uint32_t> superpixels_newidx(width,height,numimgs,1,0);
  superpixels.normalize(0,255);
//   std::cout << "writing files" << std::endl;
//   for (size_t i = 0;i < numimgs; ++i) {
//     
//     char cfilename[1024];
//     sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),(i+1)*step);
//     std::cout << "Loading " << cfilename << std::endl;
//     CImg<unsigned char> temp(cfilename);
//     CImg<unsigned char> temp2 = RGBtoGrayScale2(temp);
//     CImg<unsigned char> suppixslice = superpixels.get_slice(i);
//     suppixslice.resize_doubleXY();
//     char fname[1024];
//     sprintf(fname,"scut%d.jpg", int(i+1));
//     temp2 = temp2 / 2.0 + suppixslice / 2.0;
//     temp2.save(fname);
//   }


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
  ret_part = newpart;
  numparts  = nextpartid;

  // save the partitioning
  std::ofstream fout(outfile.c_str(), std::ofstream::binary);
  graphlab::oarchive oarc(fout);
  oarc << ret_part;


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
    features[i].resize(1, 0); 
  }

  oarc << adjlist;
  oarc << features;
  // write
  fout.close();
}

