
#include <iostream>
#include <vector>
#include <string>

#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab.hpp>



using namespace graphlab;


int main(int argc, char** argv) {
  mpi_tools::init(argc, argv);
  
  size_t procid = mpi_tools::rank();
  size_t numprocs = mpi_tools::size();

  std::vector<std::string> messages(numprocs);
  
  for(size_t i = 0; i < numprocs; ++i) {
    messages[i] = 
      "Hello from " + tostr(procid)  + " to " + tostr(i) + ".";
  }

  std::vector<std::string> results;
  mpi_tools::all2all(messages, results);
 
  for(size_t i = 0; i < results.size(); ++i) {
    std::cout << procid << ": \"" << results[i] << "\"" << std::endl;
  }
  


  mpi_tools::finalize();
  return EXIT_SUCCESS;
};
