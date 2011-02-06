// Probabilistic Reasoning Library (PRL)
// Copyright 2009 (see AUTHORS.txt for a list of contributors)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef GRAPHLAB_STL_UTIL_HPP
#define GRAPHLAB_STL_UTIL_HPP


#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>

#include <graphlab/serialization/serialize.hpp>
#include <graphlab/serialization/set.hpp>
#include <graphlab/serialization/map.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {

  // Functions on sets
  //============================================================================

  /**
   * computes the union of two sets.
   */
  template <typename T>
  std::set<T> set_union(const std::set<T>& a, const std::set<T>& b) {
    std::set<T> output;
    std::set_union(a.begin(), a.end(), 
                   b.begin(), b.end(),
                   std::inserter(output, output.begin()));
    return output;
  }
  
  template <typename T>
  std::set<T> set_union(const std::set<T>& a, const T& b) {
    std::set<T> output = a;
    output.insert(b);
    return output;
  }

  template <typename T>
  std::set<T> set_intersect(const std::set<T>& a, const std::set<T>& b) {
    std::set<T> output;
    std::set_intersection(a.begin(), a.end(), 
                          b.begin(), b.end(),
                          std::inserter(output, output.begin()));
    return output;
  }

  template <typename T>
  std::set<T> set_difference(const std::set<T>& a, const std::set<T>& b) {
    std::set<T> output;
    std::set_difference(a.begin(), a.end(), 
                        b.begin(), b.end(),
                        std::inserter(output, output.begin()));
    return output;
  }


  template <typename T>
  std::set<T> set_difference(const std::set<T>& a, const T& b) {
    std::set<T> output = a;
    output.erase(b);
    return output;
  }

  //! @return 2 sets: <s in partition, s not in partition>
  template <typename T>
  std::pair<std::set<T>,std::set<T> > 
  set_partition(const std::set<T>& s, const std::set<T>& partition) {
    std::set<T> a, b;
    a = set_intersect(s, partition);
    b = set_difference(s, partition);
    return std::make_pair(a, b);
  }

  template <typename T>
  bool set_disjoint(const std::set<T>& a, const std::set<T>& b) {
    return (intersection_size(a,b) == 0);
  }
  
  template <typename T>
  bool set_equal(const std::set<T>& a, const std::set<T>& b) {
    if (a.size() != b.size()) return false;
    return a == b; // defined in <set>
  }
  
  template <typename T>
  bool includes(const std::set<T>& a, const std::set<T>& b) {
    return std::includes(a.begin(), a.end(), b.begin(), b.end());
  }

  template <typename T>
  bool is_subset(const std::set<T>& a, const std::set<T>& b) {
    return includes(b, a);
  }

  template <typename T>
  bool is_superset(const std::set<T>& a,const std::set<T>& b) {
    return includes(a, b);
  }
  
  //! Writes a human representation of the set to the supplied stream.
  //! \relates set
  template <typename T>
  std::ostream& operator<<(std::ostream& out, const std::set<T>& s) {
    return print_range(out, s, '{', ' ', '}');
  }

  // Functions on maps
  //============================================================================

  /**
   * constant lookup in a map. assertion failure of key not found in map
   */
  template <typename Key, typename T>
  const T& safe_get(const std::map<Key, T>& map,
                    const Key& key) {
    typedef typename std::map<Key, T>::const_iterator iterator;
    iterator iter = map.find(key);
    assert(iter != map.end());
    return iter->second;
  } // end of safe_get

  /**
   * constant lookup in a map. If key is not found in map, 
   * 'default_value' is returned. Note that this can't return a reference
   * and must return a copy
   */
  template <typename Key, typename T>
  const T safe_get(const std::map<Key, T>& map,
                    const Key& key, const T default_value) {
    typedef typename std::map<Key, T>::const_iterator iterator;
    iterator iter = map.find(key);
    if (iter == map.end()) {
      return default_value;
    }
    else {
      return iter->second;
    }
  } // end of safe_get

  /**
   * Transform each key in the map using the key_map
   * transformation. The resulting map will have the form
   * output[key_map[i]] = map[i]
   */
  template <typename OldKey, typename NewKey, typename T>
  std::map<NewKey, T>
  rekey(const std::map<OldKey, T>& map,
        const std::map<OldKey, NewKey>& key_map) {
    std::map<NewKey, T> output;
    typedef std::pair<OldKey, T> pair_type;
    foreach(const pair_type& pair, map) {
      output[safe_get(key_map, pair.first)] = pair.second;
    }
    return output;
  }

  /**
   * Transform each key in the map using the key_map
   * transformation. The resulting map will have the form
   output[i] = remap[map[i]]
  */
  template <typename Key, typename OldT, typename NewT>
  std::map<Key, NewT>
  remap(const std::map<Key, OldT>& map,
        const std::map<OldT, NewT>& val_map) {
    std::map<Key, NewT> output;
    typedef std::pair<Key, OldT> pair_type;
    foreach(const pair_type& pair, map) {
      output[pair.first] = safe_get(val_map, pair.second);
    }
    return output;
  }

  /**
   * Inplace version of remap
   */
  template <typename Key, typename T>
  void remap(std::map<Key, T>& map,
             const std::map<T, T>& val_map) {
    typedef std::pair<Key, T> pair_type;
    foreach(pair_type& pair, map) {
      pair.second = safe_get(val_map, pair.second);
    }
  }

  /**
   * Computes the union of two maps
   */
  template <typename Key, typename T>
  std::map<Key, T> 
  map_union(const std::map<Key, T>& a,
            const std::map<Key, T>& b) {
    // Initialize the output map
    std::map<Key, T> output;
    std::set_union(a.begin(), a.end(),
                   b.begin(), b.end(),
                   std::inserter(output, output.begin()),
                   output.value_comp());
    return output;
  }
  
  /**
   * Computes the intersection of two maps
   */
  template <typename Key, typename T>
  std::map<Key, T> 
  map_intersect(const std::map<Key, T>& a,
                const std::map<Key, T>& b) {
    // Initialize the output map
    std::map<Key, T> output;
    // compute the intersection
    std::set_intersection(a.begin(), a.end(),
                          b.begin(), b.end(),
                          std::inserter(output, output.begin()),
                          output.value_comp());
    return output;
  }
  
  /**
   * Returns the entries of a map whose keys show up in the set keys
   */
  template <typename Key, typename T>
  std::map<Key, T> 
  map_intersect(const std::map<Key, T>& m,
                const std::set<Key>& keys) {
    std::map<Key, T> output;
    foreach(const Key& key, keys) {
      typename std::map<Key,T>::const_iterator it = m.find(key);
      if (it != m.end())
        output[key] = it->second;
    }
    return output;
  }

  /**
   * Computes the difference between two maps
   */
  template <typename Key, typename T>
  std::map<Key, T> 
  map_difference(const std::map<Key, T>& a,
                 const std::map<Key, T>& b) {
    // Initialize the output map
    std::map<Key, T> output;
    // compute the intersection
    std::set_difference(a.begin(), a.end(),
                        b.begin(), b.end(),
                        std::inserter(output, output.begin()),
                        output.value_comp());
    return output;
  }


  /**
   * Returns the set of keys in a map
   */
  template <typename Key, typename T>
  std::set<Key> keys(const std::map<Key, T>& map) {
    std::set<Key> output;
    typedef std::pair<Key, T> pair_type;
    foreach(const pair_type& pair, map) {
      output.insert(pair.first);
    }
    return output;
  }

  /**
   * Get teh set of keys in a map as a vector
   */
  template <typename Key, typename T>
  std::vector<Key> keys_as_vector(const std::map<Key, T>& map) {
    std::vector<Key> output(map.size());   
    typedef std::pair<Key, T> pair_type;
    size_t i = 0;
    foreach(const pair_type& pair, map) {
      output[i++] = pair.first;
    }
    return output;
  }


  /**
   * Gets the values from a map
   */
  template <typename Key, typename T>
  std::set<T> values(const std::map<Key, T>& map) {
    std::set<T> output;
    typedef std::pair<Key, T> pair_type;
    foreach(const pair_type& pair, map) {
      output.insert(pair.second);
    }
    return output;
  }
  
  template <typename Key, typename T>
  std::vector<T> values(const std::map<Key, T>& m, 
                        const std::set<Key>& keys) {
    std::vector<T> output;

    foreach(const Key &i, keys) {
      output.push_back(safe_get(m, i));
    }
    return output;
  }
  
  template <typename Key, typename T>
  std::vector<T> values(const std::map<Key, T>& m, 
                        const std::vector<Key>& keys) {
    std::vector<T> output;
    foreach(const Key &i, keys) {
      output.push_back(safe_get(m, i));
    }
    return output;
  }
  
  //! Creates an identity map (a map from elements to themselves)
  //! \relates map
  template <typename Key>
  std::map<Key, Key> make_identity_map(const std::set<Key>& keys) {
    std::map<Key, Key> m;
    foreach(const Key& key, keys) 
      m[key] = key;
    return m;
  }

  //! Writes a map to the supplied stream.
  //! \relates map
  template <typename Key, typename T>
  std::ostream& operator<<(std::ostream& out, const std::map<Key, T>& m) {
    out << "{";
    for (typename std::map<Key, T>::const_iterator it = m.begin(); 
         it != m.end();) {
      out << it->first << "-->" << it->second;
      if (++it != m.end()) out << " ";
    }
    out << "}";
    return out;
  }


  inline std::string trim(const std::string& str) {
    std::string::size_type pos1 = str.find_first_not_of(" \t");
    std::string::size_type pos2 = str.find_last_not_of(" \t");
    return str.substr(pos1 == std::string::npos ? 0 : pos1,
                      pos2 == std::string::npos ? str.size()-1 : pos2-pos1+1);
  }

  template<typename T>
  std::string tostr(const T &t) {
    std::stringstream strm;
    strm << t;
    return strm.str();
  }

  inline std::vector<std::string> strsplit(std::string s, std::string splitchars) {
    std::vector<std::string> ret;
    size_t pos = 0;
    while(1) {
      size_t nextpos = s.find_first_of(splitchars, pos);
      if (nextpos != std::string::npos) {
        ret.push_back(s.substr(pos, nextpos - pos));
        pos = nextpos + 1;
      }
      else {
        ret.push_back(s.substr(pos));
        break;
      }
    }
    return ret;
  }
}; // end of namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif
