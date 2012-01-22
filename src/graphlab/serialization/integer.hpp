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
  
  // ok. we align to only 2, 4 or 8 bytes.
  
  unsigned char nbytes = 8;
  // 16 bit max
  if (u <= 65535) nbytes = 2;
  else if (u <= 4294967295LL) nbytes = 4;
  
  output[10 - nbytes - 1] = nbytes | (unsigned char)(isneg > 0);
  // cast to desired length;
  
  switch(nbytes) {
    case 2: 
      { uint16_t* r = reinterpret_cast<uint16_t*>(output + 8);
        (*r) = u; 
        break;  }
    case 4:
      { uint32_t* r = reinterpret_cast<uint32_t*>(output + 6);
        (*r) = u; 
        break;  }
    case 8:
      { uint64_t* r = reinterpret_cast<uint64_t*>(output + 2);
        (*r) = u; 
        break;  }
    default:
      assert(false);
  }
  return nbytes + 1;
}




template <typename IntType>
inline void decompress_int(const char* arr, IntType &ret) {
  unsigned char len = arr[0] & 14; //(2 | 4 | 8) 
  bool isneg = arr[0] & 1;
  ret = 0;
  switch(len) {
    case 2:
      ret = *reinterpret_cast<const uint16_t*>(arr + 1); break;
    case 4:
      ret = *reinterpret_cast<const uint32_t*>(arr + 1); break;
    case 8:
      ret = *reinterpret_cast<const uint64_t*>(arr + 1); break;
    default:
      assert(false);
  }
  if (isneg)  ret = -ret;
}

template <typename IntType>
inline void decompress_int_from_ref(const char* &arr, IntType &ret) {
  unsigned char len = arr[0] & 14; //(2 | 4 | 8) 
  bool isneg = arr[0] & 1;
  ret = 0;
  switch(len) {
    case 2:
      ret = *reinterpret_cast<const uint16_t*>(arr + 1); break;
    case 4:
      ret = *reinterpret_cast<const uint32_t*>(arr + 1); break;
    case 8:
      ret = *reinterpret_cast<const uint64_t*>(arr + 1); break;
    default:
      assert(false);
  }
  if (isneg)  ret = -ret;
}


template <typename ArcType, typename IntType>
inline void decompress_int(ArcType &strm, IntType &ret) {
  if (strm.has_directbuffer()) {
    char c;
    c = strm.read_char();
    
    unsigned char len = c & 14; //(2 | 4 | 8) 
    bool isneg = c & 1;
  
    const char* arr = strm.get_direct_buffer(len);
      
    switch(len) {
      case 2:
        ret = *reinterpret_cast<const uint16_t*>(arr); break;
      case 4:
        ret = *reinterpret_cast<const uint32_t*>(arr); break;
      case 8:
        ret = *reinterpret_cast<const uint64_t*>(arr); break;
      default:
        assert(false);
    }
    if (isneg)  ret = -ret;
  }
  else {
    char c;
    strm.read(&c, 1);

        
    unsigned char len = c & 14; //(2 | 4 | 8) 
    bool isneg = c & 1;
    char arr[8];
    strm.read(arr, len);
    // hack to avoid "dereferencing type punned... warning"
    char* tmp = arr;
    switch(len) {
      case 2:
        ret = *reinterpret_cast<const uint16_t*>(tmp); break;
      case 4:
        ret = *reinterpret_cast<const uint32_t*>(tmp); break;
      case 8:
        ret = *reinterpret_cast<const uint64_t*>(tmp); break;
      default:
        assert(false);
    }
    if (isneg)  ret = -ret;
  }
}


}



#endif

