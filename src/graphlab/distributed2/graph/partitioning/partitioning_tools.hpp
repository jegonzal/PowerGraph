#ifndef GRAPHLAB_PARTITIONING_TOOLS
#define GRAPHLAB_PARTITIONING_TOOLS


#include <string>

namespace graphlab {


  namespace partitioning_tools {
    
    void construct_partitioning(int argc, char** argv, 
                                int numparts,
                                const std::string& path,
                                const double acceptance_probability);
  };
  



}; // end graphlab namspace


#endif

