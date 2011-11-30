/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef INTEGER_HPP
#define INTEGER_HPP
#include <stdint.h>
#include <graphlab/logger/assertions.hpp>
namespace graphlab {
  
/**
  * Performs a variable bit length encode of the integer u into the output array
  * function returns the length of the output.
  * NOTE: The output is written into the array aligned TO THE RIGHT.
  * in other words, if compress_int returns the value 4.
  * The actual bytes used to store the output is in output[6..9]
  */
inline unsigned char compress_int(uint64_t u, char output[10]) {
  // if 1st bit of u is set. could be negative,
  // flip all the bits if it is
  uint64_t isneg = (uint64_t)((int64_t)(u) >> 63);
  // if first bit of u is set, isneg = -1.......
  // otherwise isneg = 0....
  u = (u ^ isneg) - isneg;
  
  // get the largest bit
  unsigned char nbits = 1;  // minimum of 1 bit. even if u == 0
  if (u != 0) nbits = (unsigned char)(64 - __builtin_clzll(u));
  
  // figure out how many bytes we really need
  unsigned char b = (nbits >> 3) + ((nbits & 7) > 0);
  
  // build the base byte.
  // storage will take 'b' bytes
  unsigned char signbit = (unsigned char)(isneg > 0);
  // bits [0-2] bits of first byte of output contains 1 less of the length of the 
  // integer. 
  // bit 3 is a sign bit
  
  output[10 - b - 1] = (b - 1) | (signbit << 3);
  switch(b) {
    case 8:
      output[2] = (char)(u >> 56);
    case 7:
      output[3] = (char)(u >> 48);
    case 6:
      output[4] = (char)(u >> 40);
    case 5:
      output[5] = (char)(u >> 32);
    case 4:
      output[6] = (char)((uint32_t)(u) >> 24);
    case 3:
      output[7] = (char)((uint32_t)(u) >> 16);
    case 2:
      output[8] = (char)((uint16_t)(u) >> 8);
    case 1:
      output[9] = (char)(unsigned char)(u);
  }
  return b + 1;
}



template <typename IntType>
inline void decompress_int(const char* arr, IntType &ret) {
  bool isneg = (arr[0] & 8);
  unsigned char len = (arr[0] & 7) + 1;
  ++arr;
  ret = 0;
  while (len) {
    ret = (ret << 8) | (unsigned char)(*arr);
    --len;
    ++arr;
  };
  if (isneg)  ret = -ret;
}

template <typename IntType>
inline void decompress_int_from_ref(const char* &arr, IntType &ret) {
  bool isneg = (arr[0] & 8);
  unsigned char len = (arr[0] & 7) + 1;
  ++arr;
  ret = 0;
  while (len) {
    ret = (ret << 8) | (unsigned char)(*arr);
    --len;
    ++arr;
  };
  if (isneg)  ret = -ret;
}


template <typename IntType>
inline void decompress_int(std::istream &strm, IntType &ret) {
  char c;
  strm.read(&c, 1);
  bool isneg = (c & 8);
  unsigned char len = (c & 7) + 1;
  char buf[10];
  strm.read(buf, len);
  
  char* ptr = buf;
  ret = 0;
  while (len) {
    ret = (ret << 8) | (unsigned char)(*ptr);
    --len;
    ++ptr;
  };
  if (isneg)  ret = -ret;
}


}



#endif

