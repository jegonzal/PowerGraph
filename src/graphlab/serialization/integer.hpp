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
  
inline unsigned char compress_int(uint64_t u, char output[10]) {
  *reinterpret_cast<uint64_t*>(output + 2) = u;
  return 8;
}

template <typename IntType>
inline void decompress_int(const char* arr, IntType &ret) {
  ret = *reinterpret_cast<const IntType*>(arr);
}

template <typename IntType>
inline void decompress_int_from_ref(const char* &arr, IntType &ret) {
  ret = *reinterpret_cast<const IntType*>(arr);
}


template <typename StreamType, typename IntType>
inline void decompress_int(StreamType &strm, IntType &ret) {
  strm.read(reinterpret_cast<char*>(&ret), sizeof(IntType));
}


}



#endif

