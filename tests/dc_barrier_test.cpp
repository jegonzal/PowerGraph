#include <iostream>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>
#include <graphlab/rpc/portable.hpp>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;



void basic_call(size_t s) {
  std::cout << " ---  basic call with value " << s << std::endl;
}


int main(int argc, char ** argv) {
  
  global_logger().set_log_level(LOG_INFO);
  if (argc != 2) {
    std::cout << "insufficient args\n"; 
    return -1;
  }
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");
  //machines.push_back("127.0.0.1:10002");
  assert(machineid < machines.size());

  distributed_control dc(machines,"", machineid);
  dc_services services(dc);
  
  if (dc.procid() == 0) {
    dc.remote_call(1, basic_call, 123);
  }
  services.barrier();
  std::cout <<"!";
  if (dc.procid() == 0) {
    getchar();
  }
  services.barrier();
}
