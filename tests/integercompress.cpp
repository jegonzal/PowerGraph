#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <graphlab/serialization/integer.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/serialization_includes.hpp>


using namespace graphlab;


void check(uint64_t src, uint64_t target, char c[10], size_t len) {
  if (src != target) {
    std::cout << "expected: " << src << " = " << std::hex << src << std::endl;
    std::cout << "got : " << target << " = " << std::hex << target << std::endl;
    for (size_t i = 0;i < 10; ++i) {
      if (i == 10 - len) std::cout << " | ";
      std::cout << std::hex << (unsigned int)(unsigned char)c[i] << " ";
    }
    std::cout << std::endl;
    ASSERT_EQ(src, target);
  }
}

int main(int argc, char** argv) {
  char c[10];
 
  for (size_t i = 0;i < 100000; ++i) {
    int64_t u = ((int64_t)(rand()) << 32) + rand();
    if (rand() % 2) u = -u;
    unsigned char len = compress_int(u, c);
    int64_t u2;
    decompress_int<int64_t>(c + 10 - len, u2);
    check(u, u2, c, len);
  }
  for (size_t i = 0;i < 100000; ++i) {
    int32_t u =  rand();
    if (rand() % 2) u = -u;
    unsigned char len = compress_int(u, c);
    int64_t u2;
    decompress_int<int64_t>(c + 10 - len, u2);
    check(u, u2, c, len);
  }
  
  for (int64_t u = -100000;u < 100000; ++u) {
    unsigned char len = compress_int(u, c);
    int64_t u2;
    ASSERT_LE((int)len, 10);
    decompress_int<int64_t>(c + 10 - len, u2);
    check(u, u2, c, len);
  }
  for (int32_t u = -100000;u < 100000; ++u) {
    unsigned char len = compress_int(u, c);
    int64_t u2;
    ASSERT_LE((int)len, 10);
    decompress_int<int64_t>(c + 10 - len, u2);
    check(u, u2, c, len);
  } 
}
