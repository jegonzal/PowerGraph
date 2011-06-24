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


#ifndef GRAPHLAB_JUDY_UTIL_HPP
#define GRAPHLAB_JUDY_UTIL_HPP
#include <judyhash/judy_map_m.h>
#include <judyhash/judy_map_kdcell.h>

#include <vector>

#include <graphlab/macros_def.hpp>
namespace graphlab {
  /**
   * Get the set of keys in a map as a vector
   */
  template <typename Key, typename T, typename H, typename L, typename E>
  std::vector<Key> keys_as_vector(const judy_map_m<Key, T, H, L, E>& map) {
    std::vector<Key> output(map.size());
    typedef std::pair<Key, T> pair_type;
    size_t i = 0;
    foreach(const pair_type& pair, map) {
      output[i++] = pair.first;
    }
    return output;
  }


  template <typename Key, typename T, typename H, typename L, typename E>
  const T& safe_get(const judy_map_m<Key, T, H, L, E>& map, Key k) {
    typedef judy_map_m<Key, T, H, L, E> map_type;
    typename map_type::const_iterator i = map.find(k);
    assert(i != map.end());
    return i->second;
  }



  template <typename Key, typename T>
  std::vector<Key> keys_as_vector(const judy_map_kdcell<Key, T>& map) {
    std::vector<Key> output(map.size());
    typedef std::pair<Key, T> pair_type;
    size_t i = 0;
    foreach(const pair_type& pair, map) {
      output[i++] = pair.first;
    }
    return output;
  }


  template <typename Key, typename T>
  const T safe_get(const judy_map_kdcell<Key, T>& map, Key k) {
    typedef judy_map_kdcell<Key, T> map_type;
    typename map_type::const_iterator i = map.find(k);
    assert(i != map.end());
    return (*i).second;
  }
}
#include <graphlab/macros_undef.hpp>
#endif

