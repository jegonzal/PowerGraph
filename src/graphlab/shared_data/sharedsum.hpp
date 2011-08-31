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


#ifndef GRAPHLAB_SHAREDSUM_HPP
#define GRAPHLAB_SHAREDSUM_HPP
// #include <boost/shared_ptr.hpp>
// #include <boost/function.hpp>
// #include <boost/type_traits/function_traits.hpp>
// #include <boost/type_traits/remove_reference.hpp>

//#include <graphlab/parallel/atomic.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/logger/assertions.hpp>



namespace graphlab {

  struct isharedsum { 
    struct icache_entry { 
      uint16_t reads, writes;
      icache_entry(uint16_t reads = 0, uint16_t writes = 0) : 
        reads(reads), writes(writes) { }
      virtual ~icache_entry() { }
    };
    virtual void flush(icache_entry* entry) = 0; 
  };


  namespace sharedsum_impl {
    typedef isharedsum::icache_entry icache_entry;

    /**
     * Get the cache entry if it is available.  Returns null if its
     * not available.
     */
    icache_entry* get_cache_entry(isharedsum* id);
    
    /**
     * add a local cache entry
     */
    void add_cache_entry(isharedsum* id, icache_entry* entry_ptr);

    /**
     * Evict a particular key
     */
    bool evict(isharedsum* key);
    
    /**
     * evict an entry from the cache
     */
    void evict();

    size_t size();

  };

  template <typename T>
  class sharedsum : public isharedsum {
  public:
    typedef T value_type;
    typedef isharedsum::icache_entry icache_entry;

    struct cache_entry : public icache_entry {      
      value_type current, old;
      cache_entry(const value_type& value) : 
      current(value), old(value) { }
    };    

  private:

    cache_entry* get_cache_entry() const {
      return static_cast<cache_entry*>
        (sharedsum_impl::get_cache_entry(const_cast<sharedsum*>(this)));
    }

    cache_entry* create_cache_entry() const {
      cache_entry* entry_ptr = new cache_entry(value);
      ASSERT_NE(entry_ptr, NULL);
      sharedsum_impl::add_cache_entry(const_cast<sharedsum*>(this),
                                      entry_ptr);
      return entry_ptr;
    }


    const value_type& get_read() const {     
      cache_entry* entry_ptr = get_cache_entry();
      const bool is_cached = entry_ptr != NULL;
      if(is_cached) {
        ASSERT_NE(entry_ptr, NULL);
        if(entry_ptr->reads == lag) { 
          rwlock.readlock();
          entry_ptr->current += value - entry_ptr->old;
          entry_ptr->old = value;
          rwlock.unlock();
          entry_ptr->reads = 0;
        }
      } else {
        ASSERT_EQ(entry_ptr, NULL);
        // if it is not cached we go ahead an force the creation of
        // a cache entry
        rwlock.readlock();
        entry_ptr = create_cache_entry();
        rwlock.unlock();
      }
      ASSERT_NE(entry_ptr, NULL);
      entry_ptr->reads++;
      return entry_ptr->current; 
    } // end of get

    value_type& get_write() {    
      cache_entry* entry_ptr = get_cache_entry();
      const bool is_cached = entry_ptr != NULL;
      if(is_cached) {
        ASSERT_NE(entry_ptr, NULL);
        if(entry_ptr->writes == lag) { 
          rwlock.writelock();
          value += (entry_ptr->current - entry_ptr->old);
          entry_ptr->old = value;
          rwlock.unlock();
          entry_ptr->current = entry_ptr->old;
          entry_ptr->writes = 0;
          entry_ptr->reads = 0;
        }
      } else {
        ASSERT_EQ(entry_ptr, NULL);
        // if it is not cached we go ahead an force the creation of
        // a cache entry
        rwlock.readlock();
        entry_ptr = create_cache_entry();
        rwlock.unlock();
      }
      ASSERT_NE(entry_ptr, NULL);
      entry_ptr->writes++;
      return entry_ptr->current; 
    } // end of get 


    graphlab::rwlock rwlock;
    value_type value;
    uint16_t lag;    

  public:

    //! The eviction interface
    void flush(icache_entry* ientry_ptr) {
      ASSERT_NE(ientry_ptr, NULL);      
      cache_entry* entry_ptr = static_cast<cache_entry*>(ientry_ptr);
      rwlock.writelock();
      value += (entry_ptr->current - entry_ptr->old);
      rwlock.unlock();
      delete entry_ptr;
    }

    sharedsum(const T& value = T(), uint16_t lag = 10) : 
      value(value), lag(lag){ }
  

    const value_type& val() const {     
      return get_read();
    } // end of get 

    value_type& val() {     
      return get_write();
    } // end of get 

    operator T() const { return val(); }


    //! Assign a new value 
    void operator=(const T& new_val) {
      val() = val;
    }
    
    void operator+=(const T& other) {
      val() += other;
    }

    void operator-=(const T& other) { 
      val() -= other;
    }


    void flush() { sharedsum_impl::evict(this); }


  };


 
} 
#endif

