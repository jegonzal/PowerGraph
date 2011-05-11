#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
using namespace graphlab;


void print(std::string val) {
  std::cout << val << std::endl;
}

std::vector<int> add_one(std::vector<int> val) {
  val.push_back(1);
  return val;
}


int main(int argc, char ** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  
  if (mpi_tools::size() != 2) {
    std::cout<< "RPC Example 2: Asynchronous RPC with Built-in Serialization\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }

  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(param);
  
  
  dc.barrier();
  if (dc.procid() == 0) {
    dc.remote_call(1, print, "hello world!");
    /** Create a vector with a few elements */  
    std::vector<int> vec;
    vec.push_back(1); vec.push_back(2);
    /** Call the remote machine */  
    vec = dc.remote_request(1, add_one, vec);
    
    std::stringstream strm;
    /** Print the vector */  
    for (size_t i = 0; i < vec.size(); ++i) {
      strm << vec[i] << ", ";
    }
    strm << std::endl;
    strm.flush();
    dc.remote_call(1, print, strm.str());
  }
  dc.barrier();

  mpi_tools::finalize();
}
