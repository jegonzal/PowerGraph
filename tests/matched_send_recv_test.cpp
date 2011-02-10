#include <iostream>
#include <cstdio>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>
using namespace graphlab;

class some_object {
 public:
  dc_dist_object<some_object> rmi;
  some_object(distributed_control &dc):rmi(dc, this) { }
  void test_stuff() {
    std::string h = "hello";
    size_t i = 100;
    if (rmi.procid() == 0) {
      rmi.send_to(1, h, true);
      size_t r;
      rmi.recv_from(1, r);
      assert(r == i);
      rmi.send_to(1, h);
      rmi.recv_from(1, r);
      assert(r == i);
  
    }
    else {
      std::string r;
      rmi.recv_from(0, r, true);
      assert(r == h);
      rmi.send_to(0, i);
      rmi.recv_from(0, r);
      assert(r == h);
      rmi.send_to(0, i);
    }
  }
};

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
    dc.services().send_to(1, h, true);
    size_t r;
    dc.services().recv_from(1, r);
    assert(r == i);
    dc.services().send_to(1, h);
    dc.services().recv_from(1, r);
    assert(r == i);

  }
  else {
    std::string r;
    dc.services().recv_from(0, r, true);
    assert(r == h);
    dc.services().send_to(0, i);
    dc.services().recv_from(0, r);
    assert(r == h);
    dc.services().send_to(0, i);
  }
  some_object so(dc);
  so.test_stuff();
  dc.services().full_barrier();
  std::cout << so.rmi.calls_received() << " calls received\n";
  std::cout << so.rmi.calls_sent() << " calls sent\n";
  dc.services().barrier();
}
