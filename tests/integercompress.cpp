#include <stdint.h>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <graphlab/serialization/integer.hpp>
#include <graphlab/logger/assertions.hpp>

using namespace graphlab;
int main(int argc, char** argv) {
  char c[10];
  for (size_t i = 0;i < 100000; ++i) {
    int64_t u = ((int64_t)(rand()) << 32) + rand();
    if (rand() % 2) u = -u;
    unsigned char len = compress_int(u, c);
    int64_t u2;
    decompress_int<int64_t>(c + 10 - len, u2);
    ASSERT_EQ(u2, u);
  }
  for (size_t i = 0;i < 100000; ++i) {
    int32_t u =  rand();
    if (rand() % 2) u = -u;
    unsigned char len = compress_int(u, c);
    int64_t u2;
    decompress_int(c + 10 - len, u2);
    ASSERT_EQ(u2, u);
  }
  
  for (int64_t u = -100000;u < 100000; ++u) {
    unsigned char len = compress_int(u, c);
    int64_t u2;
    ASSERT_LE((int)len, 10);
    decompress_int(c + 10 - len, u2);
    ASSERT_EQ(u2, u);
  }
  for (int32_t u = -100000;u < 100000; ++u) {
    unsigned char len = compress_int(u, c);
    int64_t u2;
    ASSERT_LE((int)len, 10);
    decompress_int(c + 10 - len, u2);
    ASSERT_EQ(u2, u);
  } 
}