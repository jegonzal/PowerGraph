#include <iostream>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>    
#include <graphlab/rpc/dht.hpp>
#include <graphlab/logger/logger.hpp>
using namespace graphlab;

std::string randstring(size_t len) {
  std::string str;
  str.resize(len);
  const char *charset="ab";
  size_t charsetlen = 64;
  for (size_t i = 0;i < len; ++i) {
    str[i] = charset[rand()  % charsetlen];
  }
  return str;
}

int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;
  
  std::cout << "I am machine id " << dc.procid() 
            << " in " << dc.numprocs() << " machines"<<std::endl;
  dht<std::string, std::string> testdht(dc);
  
  std::vector<std::pair<std::string, std::string> > data;
  const size_t NUMSTRINGS = 10000;
  const size_t strlen[4] = {16, 128, 1024, 10240};
  // fill rate
  for (size_t l = 0; l < 4; ++l) {
    timer ti;
    if (dc.procid() == 0) {
      std::cout << "String Length = " << strlen[l] << std::endl;
      data.clear();
      for (size_t i = 0;i < NUMSTRINGS; ++i) {
        data.push_back(std::make_pair(randstring(8), randstring(strlen[l])));
      }
      std::cout << "10k random strings generated" << std::endl;
      std::cout << "Starting set" << std::endl;
     ti.start();
      for (size_t i = 0;i < NUMSTRINGS; ++i) {
        testdht.set(data[i].first, data[i].second);
        if (i % 100000 == 0) {
          std::cout << ".";
          std::cout.flush();
        }
      }
      std::cout << "10k insertions in " << ti.current_time();
    }
      dc.full_barrier();
      if (dc.procid() == 0) std::cout << "\t" << ti.current_time() << " " << double(strlen[l]*NUMSTRINGS)/ti.current_time()/1024/1024 <<  std::endl;
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
      std::cout << "10k reads in " << ti.current_time() << std::endl;
    }
    testdht.clear();
  }
  dc.barrier();
  testdht.print_stats();
  mpi_tools::finalize();
}
