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
  
  std::cout << "compressing and decompressing range [-10000000, 100000000]" << std::endl;
  ti.start();
  size_t totallen = 0;
  for (int64_t u = -10000000;u <= 100000000; ++u) {
    unsigned char len = compress_int(u, c);
    totallen += len;
    decompress_int(c + 10 - len, u2);
  }
  std::cout << "runtime: " << ti.current_time() << std::endl;
  std::cout << "accumulated length: " << totallen << std::endl;

}
