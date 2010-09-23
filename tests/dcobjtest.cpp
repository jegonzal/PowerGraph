#include <iostream>
#include <util/generics/any.hpp>
#include <distributed/dc_tcp_comm.hpp>
#include <distributed/dc.hpp>
#include <distributed/dc_dist_object.hpp>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;

class testclass{
 private:
   dc_dist_object<testclass> rpc;
   size_t id;
 public:
  testclass(distributed_control &dc, size_t id):rpc(dc, this), id(id) { }
  void hello() {
    std::cout << id<<"hello" << std::endl;
  }
  std::string gethello(size_t i) {
    std::cout << id<< ":gethello called" << std::endl;
    return "Hello";
  }
  void testcall() {
    std::cout << "calling id " << id << std::endl;
    rpc.remote_call(1, &testclass::hello);
    std::string ret = rpc.remote_request(1, &testclass::gethello, 10);
    std::cout << "replied with " << ret << std::endl;
  }
};


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
 
  std::vector<testclass*> tc;
  for (size_t i = 0;i < 10; ++i) {
    tc.push_back(new testclass(dc, i));
  }
  sleep(5);
  if (dc.procid() == 0) {
    for (size_t i = 0;i < 10; ++i) {
      tc[i]->testcall();
    }
  }
  sleep(10);
}
