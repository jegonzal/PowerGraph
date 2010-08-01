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
// Include the macro for the for each operation

using namespace cimg_library;


struct rgb{ 
  std::vector<float> features;
};
typedef graphlab::graph<rgb, float> graph_type;
typedef graphlab::types<graph_type> gl_types;


std::string imagepath = "/mnt/bigbrofs/usr3/ylow/backups/Desktop/make3d/videos/t7/";
size_t width = 480;
size_t height = 640;
std::string inputfile = "partnsh.bin";
std::string outputfile = "partnsh2.bin";
size_t numimgs = 480;
bool HALF = false;
/*
std::string imagepath = "/mnt/bigbrofs/usr3/ylow/data/soccer/";
std::string inputfile = "partsoccer.bin";
std::string outputfile = "partsoccer2.bin";
size_t width = 640/2;
size_t height = 480/2;
size_t numimgs = 107;
*/
size_t step = 1;
//size_t numparts = numimgs*20;

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

float tf(float f) {
  return 1 + 200 * exp(-f/10.0);
  //return std::max(100 - f, 1.0);
}

float tf2(float f) {
  return 1 + 200.0/(1+exp(-f));
  //return std::max(100 - f, 1.0);
}

std::vector<float> pix2feature(float h, float s, float v, float ve, float he) {
  std::vector<float> ret;
  const size_t CDIM = 4;
  const size_t EDGEDIM = 8;
  
  ret.resize(CDIM*CDIM*CDIM + EDGEDIM, 0);
  size_t hidx = size_t((h * CDIM) / 360);
  if (hidx >=CDIM) hidx = CDIM-1;
  size_t sidx = size_t((s * CDIM) / 1);
  if (sidx >=CDIM) sidx = CDIM-1;
  size_t vidx = size_t((v * CDIM) / 1);
  if (vidx >=CDIM) vidx = CDIM-1;

  ret[hidx*CDIM * CDIM + sidx * CDIM + vidx] = 1;
  if (hidx > 0) ret[(hidx - 1) *CDIM * CDIM + sidx * CDIM + vidx ] = 0.05;
  if (hidx < CDIM-1) ret[(hidx +1) *CDIM * CDIM + sidx * CDIM + vidx ] = 0.05;
  
  if (sidx > 0) ret[(hidx) *CDIM * CDIM + (sidx-1) * CDIM + vidx ] = 0.05;
  if (sidx < CDIM-1) ret[hidx  *CDIM * CDIM + (sidx+1) * CDIM + vidx ] = 0.05;
  
  if (vidx > 0) ret[(hidx) *CDIM * CDIM + sidx * CDIM + (vidx -1)] = 0.05;
  if (vidx < CDIM-1) ret[hidx  *CDIM * CDIM + sidx * CDIM + (vidx +1)] = 0.05;
  
  
  size_t eidx = size_t(EDGEDIM * sqrt(ve *ve + he * he));
  if (eidx  >=EDGEDIM) eidx  = EDGEDIM-1;
  ret[CDIM*CDIM*CDIM + eidx] = 1;
  if (eidx > 0) ret[CDIM*CDIM*CDIM + eidx - 1] = 0.1;
  if (eidx < EDGEDIM - 1) ret[CDIM*CDIM*CDIM + eidx + 1] = 0.1;

  //std::cout << h << " " << s << " " << val << " " << ve << " " << he << " " << eidx << std::endl;
  //getchar();
  return ret;
/*  
  
  ret.resize(CDIM*CDIM*CDIM+2 * EDGEDIM, 0);
  size_t hidx = size_t(h * CDIM / 360);
  if (hidx >=CDIM) hidx = CDIM-1;
  size_t sidx = size_t(s * CDIM / 100);
  if (sidx >=CDIM) sidx = CDIM-1;
  size_t validx = size_t(val * CDIM / 100);
  if (validx >=CDIM) validx = CDIM-1;
  
  ret[hidx*CDIM*CDIM + sidx * CDIM + validx] = 1;
  
  size_t veidx = size_t(ve);
  if (veidx >=EDGEDIM) veidx = EDGEDIM-1;
  size_t heidx = size_t(he);
  if (hidx >=EDGEDIM) heidx = EDGEDIM-1;
  ret[CDIM*CDIM*CDIM + veidx] = 1;
  ret[CDIM * CDIM * CDIM + EDGEDIM +heidx] = 1;
  return ret;*/
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

std::vector<float> pix2feature2(float h, float s, float v, float ve, float vh) {
  std::vector<float> ret;
  const size_t CDIM = 8;
  const size_t EDGEDIM = 9;
  
  ret.resize(CDIM + EDGEDIM + EDGEDIM, 0);
  size_t hidx = (CDIM * h) / 360;
  if (hidx >= CDIM) hidx = CDIM - 1;
  ret[hidx] = 1 ;
  
  size_t eidx = size_t(EDGEDIM * ve);
  if (eidx  >=EDGEDIM) eidx  = EDGEDIM-1;
  ret[CDIM + eidx] = 1;
  
  eidx = size_t(EDGEDIM * vh);
  if (eidx  >=EDGEDIM) eidx  = EDGEDIM-1;
  ret[CDIM +EDGEDIM+ eidx] = 1;

  //std::cout << h << " " << s << " " << val << " " << ve << " " << he << " " << eidx << std::endl;
  //getchar();
  return ret;
}
std::vector<float> pix2feature3(float rgbr, float rgbg, float rgbb,
                                float labl, float laba, float labb,
                                float hsvh, float hsvs, float hsvv, float ve, float he) {
  std::vector<float> ret;ret.resize(10);
  ret[0]=(rgbr/256.0);
  ret[1]=(rgbg/256.0);
  ret[2]=(rgbb/256.0);
  ret[3]=(labl/100.0);
  ret[4]=(laba/100.0);
  ret[5]=(labb/100.0);
  ret[6]=(hsvh/360.0);
  ret[7]=(hsvs/100.0);
  ret[8]=(hsvv/100.0);
  ret[9]=5*(sqrt(ve * ve + he * he));
  
  return ret;
}
void get_feature_counts(std::vector<uint32_t> pix2part, std::vector<std::vector<float> > &features) {
  CImg<float> vert;
  CImg<float> horz;
  std::vector<size_t> pixcount(features.size(), 0);
  size_t pixid = 0;
  for (size_t i = 0;i < numimgs; ++i) {
    char cfilename[1024];
    sprintf(cfilename,"%s%d.jpg", imagepath.c_str(),int((i+1)*step));
    std::cout << "Loading " << cfilename << std::endl;
    CImg<unsigned char> temp(cfilename);
    CImg<float> hsvtmp = temp.get_RGBtoHSV();
    CImg<float> labtmp = temp.get_RGBtoLab();
    CImgList<float> xy = ((labtmp.get_shared_plane(0)) / 100.0).get_gradient("xy", 3);
    vert = xy(0) * 4;
    horz = xy(1)* 4;
    vert.abs();
    horz.abs();
//    CImg<float> magnitude = (vert.get_mul(vert) + horz.get_mul(horz)).sqrt();
//    (CImg<unsigned char>(magnitude  * 256)).save("testv.jpg");

//    getchar();
    for (size_t j = 0;j < width; ++j) {
      for (size_t k = 0;k < height; ++k) {

        rgb vdata;
        vdata.features= pix2feature(hsvtmp(j,k,0,0),hsvtmp(j,k,0,1),hsvtmp(j,k,0,2),
                                    vert(j,k), horz(j,k));
        /*vdata.features= pix2feature3(temp(j,k,0,0),temp(j,k,0,1),temp(j,k,0,2),
                                     labtmp(j,k,0,0),labtmp(j,k,0,1),labtmp(j,k,0,2),
                                     hsvtmp(j,k,0,0),hsvtmp(j,k,0,1),hsvtmp(j,k,0,2),
                                     vert(j,k,0,1),horz(j,k,0,1));*/
        size_t superpixid = pix2part[pixid];
        assert(superpixid < pixcount.size());
        assert(superpixid < features.size());
        pixcount[superpixid]++;
        if (features[superpixid].size() == 0) features[superpixid] = vdata.features;
        else {
          assert(features[superpixid].size() == vdata.features.size());
          for (size_t j = 0;j < vdata.features.size(); ++j) {
            features[superpixid][j] += vdata.features[j];
          }
        }
        ++pixid;
      }
    }
  }
  
  for (size_t i = 0;i < features.size(); ++i) {
    for (size_t j = 0; j < features[i].size(); ++j) {
      features[i][j] /= pixcount[i];
//      std::cout << features[i][j] << " " ;
    }
//    std::cout << std::endl;
//    std::getchar();
  }
} // End of construct graph


float weight_function(const float &f) {
  return f;
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



bool parse_command_line(int argc, char** argv) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
    
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
    ("infile",  boost_po::value<std::string>(&(inputfile))->default_value(""),
     "infile")
    ("outfile",  boost_po::value<std::string>(&(outputfile))->default_value(""),
     "outfile");
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



int main(int argc, char** argv) {
  if (parse_command_line(argc,argv) == false) return 0;
  // create graph
  std::vector<uint32_t> ret_part;
  std::cout << "reading partitioning..." << std::endl;
  std::ifstream fin(inputfile.c_str());
  graphlab::iarchive iarc(fin);
  
  std::vector<std::vector<size_t> > adjlist;
  std::vector<std::vector<float> > features;
  iarc >> ret_part;
  iarc >> adjlist;
  iarc >> features;
  fin.close();
  // save the partitioning
  std::ofstream fout(outputfile.c_str());
  graphlab::oarchive oarc(fout);
  oarc << ret_part;
  oarc << adjlist;
  //oarc << features;
  //fout.close();
  //return 0 ;

  // generate the superpixels adjacency structure
 /* for (size_t i = 0;i < adjlist.size(); ++i) {
    std::cout << i << ": ";
    for (size_t j = 0; j < adjlist[i].size(); ++j) {
      std::cout << adjlist[i][j] << " ";
    }
    std::cout << std::endl;
  }*/
  // generate the superpixel features
    // histogram is compacted to 8*8*8 bins
  
  features.clear();
  features.resize(adjlist.size());

  get_feature_counts(ret_part, features);

  oarc << features;
  // write
  fout.close();
}

