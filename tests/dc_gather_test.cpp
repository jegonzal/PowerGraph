#include <iostream>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>
#include <graphlab/rpc/portable.hpp>
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
 
  assert(machineid < machines.size());

  distributed_control dc(machines,"", machineid);

  std::vector<std::string> data(dc.numprocs());
  std::string tempstring;
  tempstring.resize(10*1024*1024,'a');
  data[dc.procid()] = tempstring;

  // test gather
  std::vector<std::string> test1 = data;
  dc.services().gather(test1, 0);
  if (dc.procid() == 0) {
    for (size_t i = 0;i < test1.size(); ++i) {
      assert(data[dc.procid()] == tempstring);
    }
  }

  // test all gather
  test1 = data;
  dc.services().all_gather(test1);
  for (size_t i = 0;i < test1.size(); ++i) {
    assert(data[dc.procid()] == tempstring);
  }

  dc.services().barrier();
}
