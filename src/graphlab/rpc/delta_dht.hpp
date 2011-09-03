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


#ifndef GRAPHLAB_DELTA_DHT_HPP
#define GRAPHLAB_DELTA_DHT_HPP


#include <boost/unordered_map.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/functional/hash.hpp>


#include <graphlab/rpc/dc.hpp>
#include <graphlab/parallel/pthread_tools.hpp>






namespace graphlab {



  template<typename KeyType, typename ValueType>
  class delta_dht {
  public:
    typedef KeyType   key_type;
    typedef ValueType value_type;
    typedef size_t    size_type;    

    typedef boost::unordered_map<key_type, value_type> data_map_type;
    
    //! List hook required to use the intrusive containers
    typedef boost::instrusive::list_base_hook<
      boost::intrusive::link_mode<
        boost::intrusive::auto_unlink> > auto_unlink_hook;


    struct cache_entry : public auto_unlink_hook {
      value_type old;
      value_type current;
      size_t uses;
      cache_entry() : uses(0) { }
    };
    typedef boost::unordered_map<key_type, cache_entry> cache_map_type;

    //! LRU list used to queue cache entries
    typedef boost::intrusive::
    list<cache_entry, 
         boost::intrusive::constant_time_size<false> >  lru_list_type;

  private:

    //! The remote procedure call manager 
    mutable dc_dist_object<delta_dht> rpc;

    //! The data stored locally on this machine
    data_map_type  data_map;

    //! The cache of remote data currently present on this machine
    mutable cache_map_type cache_map;

    //! The LRU list which is used to order cache evictions
    mutable lru_list_type  lru_list;

    //! The maximum cache size
    size_t max_cache_size;

    //! The locks
    mutex data_lock;

    mutable size_t cache_requests, cache_misses;

  public:

    delta_dht(distributed_control& dc, 
              size_t max_cache_size = 1600) : 
      rpc(dc, this), 
      cache_map(max_cache_size),
      max_cache_size(max_cache_size), 
      cache_requests(0), cache_misses(0) {  }



    value_type& operator[](const key_type& key) {
      typedef cache_map_type::iterator iterator_type;
      // test for the key in the cache
      iterator_type iter = cache_map.find(key);
      if(key != iter.end()) { // found entry in cache
        cache_entry& entry = iter->second;
        entry.uses++; // count an additional read and write
        // update the LRU list
        lru_list.erase(cache_entry);
        lru_list.push_front(cache_entry);
        // return the value
        return entry.current;
      } else { // the data is not in the cache
        // Determine which machine owns the data

      }

    } // end of operator [] 
    


  }; // end of delta_dht


}; // end of namespace graphlab




#endif


