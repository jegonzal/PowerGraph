#include <iostream>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>

#include <graphlab/rpc/caching_dht.hpp>
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
  //machines.push_back("127.0.0.1:10002");
  assert(machineid < machines.size());

  distributed_control dc(machines,"", machineid);
  caching_dht<std::string, std::string> testdht(dc);
  
  dc.services().barrier();
  if (dc.procid() == 0) {
    testdht.set("hello", "world");
  }
  else {
    testdht.set("world", "hello");
  }
  dc.services().barrier();
  ASSERT_EQ(testdht.get("hello").second, std::string("world"));
  ASSERT_EQ(testdht.get("world").second, std::string("hello"));
  ASSERT_EQ(testdht.get("hello").second, std::string("world"));
  ASSERT_EQ(testdht.get("world").second, std::string("hello"));
  testdht.invalidate("hello");
  testdht.invalidate("world");
  ASSERT_EQ(testdht.get_cached("hello").second, std::string("world"));
  ASSERT_EQ(testdht.get_cached("world").second, std::string("hello"));
  ASSERT_EQ(testdht.get_cached("hello").second, std::string("world"));
  ASSERT_EQ(testdht.get_cached("world").second, std::string("hello"));
  dc.services().barrier();
  if (dc.procid() == 0) {
    testdht.set("hello", "pika");
  }
  dc.services().barrier();
  if (dc.procid() == 1) {
    ASSERT_EQ(testdht.get_cached("hello").second, std::string("world"));
    testdht.invalidate("hello");
  }
  dc.services().barrier();
  ASSERT_EQ(testdht.get_cached("hello").second, std::string("pika"));
  ASSERT_EQ(testdht.get_cached("hello").second, std::string("pika"));
  dc.services().barrier();
}
