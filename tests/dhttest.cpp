#include <iostream>
#include <util/generics/any.hpp>
#include <distributed/dc_tcp_comm.hpp>
#include <distributed/dc.hpp>
#include <distributed/dc_services.hpp>

#include <distributed/dht.hpp>
#include <distributed/portable.hpp>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;



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
  dht<std::string, std::string> testdht(dc);
  services.barrier();
  if (dc.procid() == 0) {
    testdht.set("hello", "world");
  }
  else {
    testdht.set("world", "hello");
  }
  services.barrier();
  std::cout << "hello --> " << testdht.get("hello").second << std::endl;
  std::cout << "world --> " << testdht.get("world").second << std::endl;
  services.barrier();
  
}
