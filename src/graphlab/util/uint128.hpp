/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_UINT128_HPP
#define GRAPHLAB_UINT128_HPP
#include <stdint.h>
#include <iostream>
#include <iomanip>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
/**
 * A 128 bit numeric type. This type is a union of a 16-byte character array (bytes),
 * and struct of two 64-bit integers (ints.high and ints.low).
 */
union gl_uint128_t {
  struct {
    uint64_t high;
    uint64_t low;  
  } ints;
  char bytes[16];
  
  gl_uint128_t() { }
  
  /**
   * Constructs a 128-bit type from a 64-bit value.
   * It simply clears the "high" 64 bits of the 128-bit integer, and sets
   * the low 64-bits to the input
   */
  explicit gl_uint128_t(uint64_t val) {
    ints.high = 0;
    ints.low = val;
  }
};

/**
 * Sets all 128bits of the the gl_uint128_t to 'true'.
 * Or the 128-bit integer representation of "-1"
 */
inline gl_uint128_t fill_128b() {
  gl_uint128_t i;
  i.ints.high = (uint64_t)(-1);
  i.ints.low = (uint64_t)(-1);
  return i;
}

/**
 * Prints the 128-bit integer as hexadecimal
 */
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