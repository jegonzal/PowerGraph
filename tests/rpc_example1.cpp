#include <iostream>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
using namespace graphlab;


void print(int val) {
  std::cout << val << std::endl;
}

int add_one(int val) {
  return val + 1;
}


int main(int argc, char ** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  
  if (mpi_tools::size() != 2) {
    std::cout<< "RPC Example 1: Basic Synchronous RPC\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
  // set up parameters
  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  // not necessary. Set log level. Change to LOG_WARNING to get less output
  global_logger().set_log_level(LOG_INFO);
  // create distributed control
  distributed_control dc(param);
  
  if (dc.procid() == 0) {
    dc.control_call(1, print, 10);
    std::cout << "5 plus 1 is : " << dc.remote_request(1, add_one, 5) << std::endl;
    std::cout << "11 plus 1 is : " << dc.remote_request(1, add_one, 11) << std::endl;
  }
  dc.barrier();
  // terminate MPI
  mpi_tools::finalize();
}
