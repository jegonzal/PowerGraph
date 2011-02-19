#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <graphlab/serialization/integer.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/logger/assertions.hpp>

using namespace graphlab;




int main(int argc, char** argv) {
  timer ti;
  ti.start();
  char c[10];
  int64_t u2;
  
  for (int64_t u = -10000000;u < 100000000; ++u) {
    unsigned char len = compress_int2(u, c);
    decompress_int2(c + 10 - len, u2);
  }
  std::cout << ti.current_time() << std::endl;
  
  ti.start();
  for (int64_t u = -10000000;u < 100000000; ++u) {
    unsigned char len = compress_int(u, c);
    decompress_int(c + 10 - len, u2);
  }
  std::cout << ti.current_time() << std::endl;

}
