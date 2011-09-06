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


/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */

#ifndef GRAPHLAB_CACHE_HPP
#define GRAPHLAB_CACHE_HPP

#include <algorithm>

#include <boost/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>


#include <graphlab/logger/assertions.hpp>


namespace graphlab {
  namespace cache { 

    // template<typename Cache, typename Source>
    // struct bind {
    //   typedef Cache cache_type;
    //   typedef typename cache_type::key_type key_type;
    //   typedef typename cache_type::value_type value_type;
    //   cache_type cache;
    //   Source& source;
    //   bind(Source& source, size_t capacity = 100) : 
    //     source(source), capacity(capacity) { }
    //   value_type get(const key_type& key) { return cache.get(key, source); }
    // }; // end of bind


    template<typename Key, typename Value>
    class lru {
    public:
      typedef Key key_type;
      typedef Value value_type;
      
      typedef boost::bimaps::bimap<
        boost::bimaps::unordered_set_of<key_type>, 
        boost::bimaps::list_of<value_type> > 
      cache_map_type;
      
      
    private:
      mutable cache_map_type cache_map;
      
      
    public:

      lru(size_t cache_reserve = 1024) {
        cache_map.left.rehash(cache_reserve); 
      }

      
      size_t size() { return cache_map.size(); }

      typename cache_map_type::iterator begin() { return cache_map.begin(); }
      typename cache_map_type::iterator end() { return cache_map.end(); }
      
      std::pair<key_type, value_type> evict() {
        ASSERT_FALSE(cache_map.empty());
        typedef typename cache_map_type::right_iterator iterator_type;
        iterator_type iter = cache_map.right.begin();
        const std::pair<key_type, value_type> 
          result(iter->get_left(), iter->get_right());
        cache_map.right.erase(iter);
        return result;
      } // end of evict

      std::pair<bool, value_type> evict(const key_type& key) {
        typedef typename cache_map_type::left_iterator iterator_type;
        iterator_type iter = cache_map.left.find(key);
        if(iter == cache_map.left.end()) 
          return std::make_pair(false, value_type());
        const value_type result = iter->get_right();
        cache_map.left.erase(iter);
        return std::make_pair(true, result);  
      } // end of evict(key)

      bool evict(const key_type& key, value_type& ret_value) {
        typedef typename cache_map_type::left_iterator iterator_type;
        iterator_type iter = cache_map.left.find(key);
        if(iter == cache_map.left.end()) return false;
        ret_value = iter->get_right();
        cache_map.left.erase(iter);
        return true;
      } // end of evict(key)
      

      bool contains(const key_type& key) const {
        typedef typename cache_map_type::left_const_iterator iterator_type;
        iterator_type iter = cache_map.left.find(key);
        return iter != cache_map.left.end();
      } // end of contains


      value_type& operator[](const key_type& key) {
        typedef typename cache_map_type::left_iterator iterator_type;
        iterator_type iter = cache_map.left.find(key);
        if(iter != cache_map.left.end()) { // already in cache
          // move it to the end
          cache_map.right.relocate(cache_map.right.end(), 
                                   cache_map.project_right(iter));
          return iter->get_right();
        } else {
          // add it to the cache
          // Get the true entry from the source
          typedef typename cache_map_type::value_type pair_type;
          cache_map.insert(pair_type(key, value_type()));
          return cache_map.left[key];
        }
      } // end of oeprator[]

      const value_type& operator[](const key_type& key) const {
        typedef typename cache_map_type::const_left_iterator iterator_type;
        iterator_type iter = cache_map.left.find(key);
        if(iter != cache_map.left.end()) { // already in cache
          // move it to the end
          cache_map.right.relocate(cache_map.right.end(), 
                                   cache_map.project_right(iter));
          return iter->get_right();
        }
        logstream(LOG_FATAL) << "Key not found!" << std::endl;
        return value_type();
      } // end of oeprator[]

      bool get(const key_type& key, value_type& ret_value) {
        typedef typename cache_map_type::left_iterator iterator_type;
        iterator_type iter = cache_map.left.find(key);
        if(iter != cache_map.left.end()) { // already in cache
          ret_value = iter->get_right();
          // move it to the end
          cache_map.right.relocate(cache_map.right.end(), 
                                   cache_map.project_right(iter));
          return true;
         } else return false;
      } // end of get

    }; // end of class lru




  }; // end of cache namespace 
}; // end of graphlab namespace

#endif



