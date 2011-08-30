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
    virtual void flush(icache_entry* entry); 
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

    cache_entry* get_cache_entry() {
      return static_cast<cache_entry*>
        (sharedsum_impl::get_cache_entry(this));
    }

    cache_entry* create_cache_entry() {
      cache_entry* entry_ptr = new cache_entry(value);
      ASSERT_NE(entry_ptr, NULL);
      sharedsum_impl::add_cache_entry(this, entry_ptr);
      return entry_ptr;
    }


    const value_type& get_read() const {     
      cache_entry* entry_ptr = get_cache_entry();
      const bool is_cached = entry_ptr != NULL;
      const bool grabbed_readlock = rwlock.try_readlock();
      // Base case
      if(!is_cached && grabbed_readlock) {
        return value;
      } else if(is_cached && grabbed_readlock) {
        entry_ptr->current += value - entry_ptr->old;
        entry_ptr->old = value;
        rwlock.unlock();
        entry_ptr->reads = 1;
        return entry_ptr->current;
      } else if(is_cached && !grabbed_readlock) {
        if(entry_ptr.reads == 10) { 
          rwlock.readlock();
          entry_ptr->current += value - entry_ptr->old;
          entry_ptr->old = value;
          rwlock.unlock();
          entry_ptr->reads = 0;
        }
        entry_ptr->reads++;
        return entry_ptr->current; 
      } else {
        ASSERT_FALSE(is_cached);
        ASSERT_FALSE(grabbed_readlock);
        // if it is not cached we go ahead an force the creation of
        // a cache entry
        rwlock.readlock();
        entry_ptr = create_cache_entry();
        rwlock.unlock();
        entry_ptr->reads++;
        return entry_ptr->current; 
      }
    } // end of get 

    value_type& get_write() const {     
      cache_entry* entry_ptr = get_cache_entry();
      const bool is_cached = entry_ptr != NULL;
      const bool grabbed_writelock = rwlock.try_writelock();
      // Base case
      if(!is_cached && grabbed_writelock) {
        return value;
      } else if(is_cached && grabbed_writelock) {
        value += (entry_ptr->current - entry_ptr->old);
        entry_ptr->old = value;
        rwlock.unlock();
        entry_ptr->current = entry_ptr->old;
        entry_ptr->writes = 1;
        entry_ptr->reads = 0;
        return entry_ptr->current;
      } else if(is_cached && !grabbed_writelock) {
        if(entry_ptr.writes == 10) { 
          rwlock.writelock();
          value += (entry_ptr->current - entry_ptr->old);
          entry_ptr->old = value;
          rwlock.unlock();
          entry_ptr->current = entry_ptr->old;
          entry_ptr->writes = 0;
          entry_ptr->reads = 0;
        }
        entry_ptr->writes++;
        return entry_ptr->current; 
      } else {
        ASSERT_FALSE(is_cached);
        ASSERT_FALSE(grabbed_writelock);
        // if it is not cached we go ahead an force the creation of
        // a cache entry
        rwlock.readlock();
        entry_ptr = create_cache_entry();
        rwlock.unlock();
        entry_ptr->writes++;
        return entry_ptr->current; 
      }
    } // end of get 


    void release() {
      if(get_cache_entry() == NULL) rwlock.unlock();
    }

    graphlab::rwlock rwlock;
    value_type value;
    

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

    sharedsum(const T& value = T()) : value(value) { }
  

    value_type get() const {     
      const value_type ret_val = get_read();
      release();
      return ret_val;
    } // end of get 

    operator T() { return get(); }


    //! Assign a new value 
    void operator=(const T& val) {
      value_type& cache_val = get_write();
      cache_val = val;
      release();
    }
    
    void operator+=(const T& val) { 
      value_type& cache_val = get_write();
      cache_val += val;
      release();
    }

    void operator-=(const T& val) { 
      value_type& cache_val = get_write();
      cache_val -= val;
      release();
    }


    void flush() { sharedsum_impl::evict(this); }


  };


 
} 
#endif

