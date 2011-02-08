#ifndef GRAPHLAB_ZOLTAN_ATOM_BUILDER
#define GRAPHLAB_ZOLTAN_ATOM_BUILDER


#include <string>

namespace graphlab {


  void construct_partitioning(int argc, char** argv, 
                              int numparts,
                              const std::string& path);
  



}; // end graphlab namspace


#endif

