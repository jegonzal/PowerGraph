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
#include <boost/functional/hash.hpp>


#include <graphlab/rpc/dc.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/cache.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {



  template<typename KeyType, typename ValueType>
  class delta_dht {
  public:
    typedef KeyType   key_type;
    typedef ValueType value_type;
    typedef size_t    size_type;    

    typedef boost::unordered_map<key_type, value_type> data_map_type;
    

    struct cache_entry {
      value_type old;
      value_type current;
      size_t uses;
      cache_entry() : uses(0) { }
    };
    typedef cache::lru<key_type, cache_entry> cache_type;

  private:

    //! The remote procedure call manager 
    mutable dc_dist_object<delta_dht> rpc;

    //! The data stored locally on this machine
    data_map_type  data_map;

    //! The cache of remote data currently present on this machine
    mutable cache_type cache;

    //! The maximum cache size
    size_t max_cache_size;

    //! The maximum number of uses
    size_t max_uses;

    //! The locks
    mutex data_lock;

    mutable size_t cache_hits, cache_misses;

    boost::hash<key_type> hasher;


  public:

    delta_dht(distributed_control& dc, 
              size_t max_cache_size = 2056) : 
      rpc(dc, this), 
      max_cache_size(max_cache_size), max_uses(100), 
      cache_hits(0), cache_misses(0) {  
      rpc.full_barrier();
    }

    ~delta_dht() {
      rpc.full_barrier();
      std::cout << "cache_hits: " << cache_hits << std::endl;
      std::cout << "cache_misses: " << cache_misses << std::endl;
    }

    value_type& operator[](const key_type& key) {
      // test for the key in the cache
      if(!cache.contains(key)) {
        cache_misses++;
        // make room for the new entry
        while(cache.size() + 1 > max_cache_size) {
          const std::pair<key_type, cache_entry> pair = cache.evict();
          const key_type& key = pair.first;
          const cache_entry& entry = pair.second;
          send_delta(key, entry.current - entry.old);
        }
        // get the new entry from the server
        cache_entry& entry = cache[key];
        entry.old = get_master(key);
        entry.current = entry.old;
        return entry.current;
      } else {
        cache_hits++;
      }

      cache_entry& entry = cache[key];
      if(entry.uses > max_uses) {
        synchronize(key, entry);
        entry.uses = 0;
      }
      return entry.current;
    } // end of operator []


    //! empty the local cache
    void flush() {
      while(cache.size() > 0) {
        const std::pair<key_type, cache_entry> pair = cache.evict();
        const key_type& key = pair.first;
        const cache_entry& entry = pair.second;
        send_delta(key, entry.current - entry.old);
      }
    }

    void synchronize(const key_type& key) {
      if(cache.contains(key)) synchronize(key, cache[key]);
    }


    size_t owning_cpu(const key_type& key) const {
      const size_t hash_value = hasher(key);
      const size_t cpuid = hash_value % rpc.numprocs();
      return cpuid;
    }
      


    bool is_local(const key_type& key) const {
      return owning_cpu(key) == rpc.procid();
    } // end of is local

   

    value_type delta(const key_type& key) const { 
      if(cache.contains(key)) { 
        cache_entry& entry = cache[key];
        return entry.current - entry.old;
      } else {
        return value_type();
      }
    }

    void send_delta(const key_type& key) {
      if(cache.contains(key)) send_delta(key, delta(key));
    }


    void send_delta(const key_type& key, const value_type& delta) {
      send_delta_rpc(key, delta);
    } // end of send_delta


    size_t local_size() const {
      data_lock.lock();
      const size_t result = data_map.size(); 
      data_lock.unlock();
      return result;
    }

    size_t size() const {
      std::vector<size_t> sizes(rpc.numprocs());
      sizes[rpc.procid()] = local_size();
      rpc.all_gather(sizes);
      size_t sum = 0;
      foreach(const size_t& val, sizes) sum += val;
      return sum;
    }

    size_t numprocs() const { return rpc.num_procs(); }
    size_t procid() const { return rpc.procid(); }


    value_type get_master(const key_type& key) {
      // If the data is stored locally just read and return
      if(is_local(key)) {
        data_lock.lock();
        const value_type ret_value = data_map[key];
        data_lock.unlock();
        return ret_value;
      } else {
        return rpc.remote_request(owning_cpu(key), 
                                  &delta_dht::get_master, key);
      }
    } // end of direct get

  private:
    void send_delta_rpc(const key_type& key, const value_type& delta)  {
      // If the data is stored locally just read and return
      if(is_local(key)) {
        data_lock.lock();
        typename data_map_type::iterator iter = data_map.find(key);
        ASSERT_TRUE(iter != data_map.end());
        iter->second += delta;
        data_lock.unlock();
      } else {
        rpc.remote_call(owning_cpu(key), 
                        &delta_dht::send_delta_rpc, 
                        key, delta);
      }
    } // end of send_delta_rpc
    
    void synchronize(const key_type& key, cache_entry& entry)  {
      const value_type delta = entry.current - entry.old;
      entry.old = synchronize_rpc(key, delta);
      entry.current = entry.old;
    } // end of synchronize

    value_type synchronize_rpc(const key_type& key, const value_type& delta) {
      if(is_local(key)) {
        data_lock.lock();
        typename data_map_type::iterator iter = data_map.find(key);
        ASSERT_TRUE(iter != data_map.end());
        const value_type ret_value = (iter->second += delta);
        data_lock.unlock();
        return ret_value;
      } else {
        return rpc.remote_request(owning_cpu(key), 
                                  &delta_dht::synchronize_rpc, 
                                  key, delta);
      }
    } // end of synchronize_rpc





  }; // end of delta_dht


}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>



#endif


