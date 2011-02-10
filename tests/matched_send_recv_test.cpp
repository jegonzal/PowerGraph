#include <iostream>
#include <cstdio>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>
using namespace graphlab;


int main(int argc, char ** argv) {
  /** Initialization */
  global_logger().set_log_level(LOG_INFO);
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");

  distributed_control dc(machines,"", machineid, 8, TCP_COMM);
  std::string h = "hello";
  size_t i = 100;
  if (dc.procid() == 0) {
    dc.send_to(1, h);
    size_t r;
    dc.recv_from(1, r);
    assert(r == i);
    dc.send_to(1, h);
    dc.recv_from(1, r);
    assert(r == i);

  }
  else {
    std::string r;
    dc.recv_from(0, r);
    assert(r == h);
    dc.send_to(0, i);
    dc.recv_from(0, r);
    assert(r == h);
    dc.send_to(0, i);
  }
  dc.services().barrier();
}
