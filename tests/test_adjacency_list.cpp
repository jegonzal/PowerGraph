
#include <iostream>
#include <string>
#include <graphlab.hpp>
#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>


int main(int argc, char** argv) {
  std::string alist_fname(argv[1]);
  
  graphlab::adjacency_list alist;
  alist.load(alist_fname);
  alist.save(argv[2], 0);

  return EXIT_SUCCESS;

}
