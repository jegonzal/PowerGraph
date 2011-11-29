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
  

  ti.start();
  size_t totallen = 0;
  int64_t rangelow = -100000000;
  int64_t rangehigh = 100000000;
  std::cout << "compressing and decompressing range [" << rangelow << ", "<< rangehigh << "]" << std::endl;
  for (int64_t u = rangelow;u <= rangehigh; ++u) {
    unsigned char len = compress_int(u, c);
    totallen += len;
    decompress_int(c + 10 - len, u2);
  }
  std::cout << "runtime: " << ti.current_time() << std::endl;
  std::cout << "Range compression rate: " << (8 * (rangehigh - rangelow) / ti.current_time())/(1024 * 1024) << " MBps" << std::endl;
  std::cout << "Compression Ratio: " << (double)totallen / (8 * (rangehigh - rangelow)) << std::endl;

}
