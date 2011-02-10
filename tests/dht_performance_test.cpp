#include <iostream>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc_tcp_comm.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_services.hpp>

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
  if (argc < 2) {
    std::cout << "insufficient args\n"; 
    return -1;
  }
  size_t machineid = atoi(argv[1]);
  std::vector<std::string> machines;
  for (size_t i = 2;i < (size_t)argc; ++i ) {
    machines.push_back(argv[i]);
  }
  //machines.push_back("127.0.0.1:10002");
  assert(machineid < machines.size());
  distributed_control dc(machines,"", machineid);
  std::cout << "I am machine id " << dc.procid() 
            << " in " << dc.numprocs() << " machines"<<std::endl;
  dc_services services(dc);
  dht<std::string, std::string> testdht(dc);
  services.barrier();
  
  std::vector<std::pair<std::string, std::string> > data;
  const size_t NUMSTRINGS = 100;
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
      if (i % 10 == 0) {
        std::cout << ".";
        std::cout.flush();
      }
    }
    std::cout << "100k insertions in " << ti.current_time() << std::endl;
  }
    services.comm_barrier();
  services.barrier();
  // get rate
  if (dc.procid() == 0) {
    std::cout << "Starting get" << std::endl;

    timer ti;
    ti.start();
    for (size_t i = 0;i < NUMSTRINGS; ++i) {
      std::pair<bool, std::string> ret = testdht.get(data[i].first);
      assert(ret.first);
      if (i % 1 == 0) {
        std::cout << ".";
        std::cout.flush();
      }
    }
    std::cout << "100k reads in " << ti.current_time() << std::endl;
  }
  services.barrier();
  testdht.print_stats();
}
