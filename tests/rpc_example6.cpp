#include <iostream>
#include <string>
#include <vector>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/util/generics/any.hpp>
using namespace graphlab;


void print(any val) {
  val.print(std::cout);
  std::cout << std::endl;
}


int main(int argc, char ** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  
  if (mpi_tools::size() != 2) {
    std::cout<< "RPC Example 6: Asynchronous RPC with any \n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }

  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(param);
  
  if (dc.procid() == 0) {
    dc.remote_call(1, print, 15);
    dc.remote_call(1, print, 10.5);
    dc.remote_call(1, print, std::string("hello world"));    
  }
    
  int i = dc.procid() == 0 ? 10 : 100;
  dc.broadcast(i, dc.procid() == 0);
  std::cout << i << std::endl;
  assert(i == 10);
  
  mpi_tools::finalize();
}
