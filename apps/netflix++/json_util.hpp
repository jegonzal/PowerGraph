/*  
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


#include <boost/format.hpp>
#include <boost/algorithm/string/find_format.hpp>
#include <boost/algorithm/string/join.hpp>
template<typename Map>
std::string map2json(const Map& map) {
  typedef typename Map::value_type kv_type;
  typedef typename Map::const_iterator iter_type;
  std::stringstream ss;
  ss << "{\n";
  iter_type iter = map.begin();
  if (iter != map.end()) {
    ss << "\"" << (iter->first) << "\": "
       << "\"" << (iter->second) << "\"";
    ++iter;
    while (iter != map.end()) {
      ss << ",\n"
         << "\"" << (iter->first) << "\": "
         << "\"" << (iter->second) << "\"";
      ++iter;
    }
  }
  ss << "}\n";
  return ss.str();
}


template<typename Array>
std::string arr2json(const Array& arr) {
  std::stringstream ss;
  ss << "[\n";
  for (size_t i = 0; i < arr.size(); ++i) {
    if (i) 
      ss << ", ";
    ss << "\"" << (arr[i]) << "\"";
  }
  ss << "]\n";
  return ss.str();
}
