#ifndef GRAPHLAB_UINT128_HPP
#define GRAPHLAB_UINT128_HPP
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
union gl_uint128_t {
  struct {
    uint64_t high;
    uint64_t low;  
  } ints;
  char bytes[16];
  
  gl_uint128_t() { }
  
  explicit gl_uint128_t(uint64_t val) {
    ints.high = 0;
    ints.low = val;
  }
};

inline gl_uint128_t fill_128b() {
  gl_uint128_t i;
  i.ints.high = (uint64_t)(-1);
  i.ints.low = (uint64_t)(-1);
  return i;
}

inline std::ostream& operator<<(std::ostream& out, const gl_uint128_t &val) {
  static char hexchar[17] = "0123456789abcdef";
  
  for (size_t i = 0;i < 16; ++i) {
    out << hexchar[(val.bytes[i] >> 4) & 15];
    out << hexchar[val.bytes[i] & 15];
  }
  return out;
}

}

SERIALIZABLE_POD(graphlab::gl_uint128_t);

#endif