#include <iostream>
#include <string>
#include <map>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
using namespace graphlab;

struct teststruct{
  int a;
  double b;
};
SERIALIZABLE_POD(teststruct);


void print(std::map<int, teststruct> &data1,  
           std::string data2) {
  std::cout << "1.a = " << data1[1].a << std::endl;
  std::cout << "10.b = " << data1[10].b << std::endl;
  std::cout << "string = " << data2 << std::endl;
}



int main(int argc, char ** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  
  
  if (mpi_tools::size() != 2) {
    std::cout<< "RPC Example 3: Asynchronous RPC with Struct POD Serialization\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
  
  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  global_logger().set_log_level(LOG_INFO);
  distributed_control dc(param);
  
    
  if (dc.procid() == 0) {
    std::map<int, teststruct> data;
    data[1].a = 10;
    data[2].b = 15.0;
    dc.remote_call(1, print, data, "hello world!");
  }
  dc.barrier();

  mpi_tools::finalize();
}
