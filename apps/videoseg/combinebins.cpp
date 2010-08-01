#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <graphlab.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <boost/program_options.hpp>
#define cimg_display 0
#include "CImg.h"
// Include the macro for the for each operation

using namespace cimg_library;


struct rgb{ 
  std::vector<float> features;
};
typedef graphlab::graph<rgb, float> graph_type;
typedef graphlab::types<graph_type> gl_types;


std::vector<std::string> inputfiles;
std::string outputfile = "partnsh2.bin";

bool parse_command_line(int argc, char** argv) {
  // Because typing is painful
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
    desc("Denoise a randomly generated image using Gibbs Sampling.");
    
  // Set the program options
  desc.add_options()
    ("help",   "produce this help message")
    ("infiles",  boost_po::value<std::vector<std::string> >(&(inputfiles))->default_value(std::vector<std::string>(),""),
     "infiles")
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
  std::vector<std::vector<size_t> > adjlist;
  std::vector<std::vector<float> > features;
  
  for (size_t i = 0 ;i  < inputfiles.size(); ++i) {
    std::vector<uint32_t> localret_part;
    std::vector<std::vector<size_t> > localadjlist;
    std::vector<std::vector<float> > localfeatures;
    size_t renumberfrom  = adjlist.size();
    std::cout << "reading " << inputfiles[i] << std::endl;
    std::ifstream fin(inputfiles[i].c_str());
    graphlab::iarchive iarc(fin);
  
    iarc >> localret_part;
    iarc >> localadjlist;
    iarc >> localfeatures;
    fin.close();
  
    for (size_t j  = 0;j < localret_part.size(); ++j) {
      localret_part[j] += renumberfrom;
    }
    for (size_t j  = 0;j < localadjlist.size(); ++j) {
      for (size_t k = 0;k < localadjlist[j].size(); ++k) {
        localadjlist[j][k] += renumberfrom;
      }
    }
    // insert into ret_part etc
    std::copy(localret_part.begin(), localret_part.end(), std::back_inserter(ret_part));
    std::copy(localadjlist.begin(), localadjlist.end(), std::back_inserter(adjlist));
    std::copy(localfeatures.begin(), localfeatures.end(), std::back_inserter(features));
  }
  
  // save the partitioning
  std::ofstream fout(outputfile.c_str());
  graphlab::oarchive oarc(fout);
  oarc << ret_part;
  oarc << adjlist;
  oarc << features;
  fout.close();
  std::cout << adjlist.size() << " frames processed" << std::endl;
}

