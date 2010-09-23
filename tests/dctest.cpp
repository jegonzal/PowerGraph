#include <iostream>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/portable.hpp>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;



void basic_call(size_t s) {
  std::cout << " ---  basic call with value " << s << std::endl;
}

void remote_call(distributed_control &dc, procid_t src, size_t s) {
  std::cout << " ---  remote call with value " << s << std::endl;
}


std::string basic_request(std::string s) {
  std::cout << " ---  basic request with value " << s << std::endl;
  return "boo";
}


std::string any_request(any s) {
  std::cout << " ---  any_request call with value ";
  s.print(std::cout);
  std::cout << std::endl;
  return "boo";
}

std::string portable_call_test(distributed_control &dc, procid_t src, 
                          std::string s) {
  std::cout << " ---  portable call with value " << s << std::endl;
  return "pikachu";
}


struct teststruct {
  std::string s;
  void save(oarchive& oarc) const {
    oarc << s;
  }
 void load(iarchive& iarc) {
    iarc >> s;
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
  REGISTER_RPC(dc, portable_call_test);
  REGISTER_RPC(dc, basic_request);
  REGISTER_RPC(dc, basic_call);
  if (dc.procid() == 0) {
      dc.remote_call(1, printf, "hello %d = %f\n", 10, 20.0);
  /*  dc.remote_call(1, basic_call, 123);
    dc.remote_call(1, remote_call, 123);
        
    std::string ret = dc.remote_request(1, basic_request, "hello world");
    
  //  std::cout << ret << std::endl;
  //      dc.remote_request(1, remote_call, 123);

    std::cout << "portable request 1 returned: " << dc.remote_request(1, PORTABLE(portable_call_test), "hello") << std::endl;
    std::cout << "portable request 2 returned: " << dc.remote_request(1, PORTABLE(basic_request), "hello") << std::endl;
    dc.remote_call(1, PORTABLE(basic_request), "again! again!");
    dc.remote_request(1, PORTABLE(basic_call), 'a');
    
    //ret = dc.remote_request(1, any_request, 123);
  //  ret = dc.remote_request(1, any_request, 123.1);

  //  std::cout << ret << std::endl;*/
  }
   sleep(10);
}
