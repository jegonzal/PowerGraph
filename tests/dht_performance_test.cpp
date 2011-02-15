#include <iostream>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
    
#include <graphlab/rpc/dht.hpp>
#include <graphlab/rpc/portable.hpp>
#include <graphlab/serialization/podify.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;

std::string randstring(size_t len) {
  std::string str;
  str.resize(len);
  const char *charset="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890,.";
  size_t charsetlen = 64;
  for (size_t i = 0;i < len; ++i) {
    str[i] = charset[rand()  % charsetlen];
  }
  return str;
}

int main(int argc, char ** argv) {
  
  global_logger().set_log_level(LOG_INFO);

  dc_init_param param;
  if (init_param_from_env(param) == false) {
    return 0;
  }
  
  global_logger().set_log_level(LOG_DEBUG);

  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
  distributed_control dc(param);
  std::cout << "I am machine id " << dc.procid() 
            << " in " << dc.numprocs() << " machines"<<std::endl;
  dht<std::string, std::string> testdht(dc);
  
  std::vector<std::pair<std::string, std::string> > data;
  const size_t NUMSTRINGS = 10000;
  // fill rate
  if (dc.procid() == 0) {
    for (size_t i = 0;i < NUMSTRINGS; ++i) {
      data.push_back(std::make_pair(randstring(8), randstring(10*1024)));
    }
    std::cout << "100k random strings generated" << std::endl;
    std::cout << "Starting set" << std::endl;
    timer ti;
    ti.start();
    for (size_t i = 0;i < NUMSTRINGS; ++i) {
      testdht.set(data[i].first, data[i].second);
      if (i % 100 == 0) {
        std::cout << ".";
        std::cout.flush();
      }
    }
    std::cout << "100k insertions in " << ti.current_time() << std::endl;
  }
    dc.full_barrier();
//  dc.barrier();
  // get rate
  if (dc.procid() == 0) {
    std::cout << "Starting get" << std::endl;

    timer ti;
    ti.start();
    for (size_t i = 0;i < NUMSTRINGS; ++i) {
      std::pair<bool, std::string> ret = testdht.get(data[i].first);
      assert(ret.first);
      if (i % 100 == 0) {
        std::cout << ".";
        std::cout.flush();
      }
    }
    std::cout << "100k reads in " << ti.current_time() << std::endl;
  }
  dc.barrier();
  testdht.print_stats();
}
