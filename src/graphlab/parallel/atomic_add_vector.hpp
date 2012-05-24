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



#ifndef GRAPHLAB_ATOMIC_ADD_VECTOR_HPP
#define GRAPHLAB_ATOMIC_ADD_VECTOR_HPP


#include <vector>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/lock_free_pool.hpp>



namespace graphlab {

  /**
   * \TODO DOCUMENT THIS CLASS
   */ 
  
  template<typename ValueType>
  class atomic_add_vector {
  public:
    typedef ValueType value_type;

    // We needed a second "NULL" pointer value to indicate a value
    // that is being swapped.
#define VALUE_PENDING (value_type*)(size_t)(-1)
    
  private:
    lock_free_pool<value_type> pool;
    atomic<size_t> joincounter;
    
    class atomic_box_type {
    private:
      value_type* volatile value_ptr;
    public:
      atomic_box_type() : value_ptr(NULL) {  }
      
      // void assign_unsafe(const atomic_box_type &other) {
      //   value_ptr = other.value_ptr;
      // }
      
      // bool peek_unsafe(const value_type& other) {
      //   if (value_ptr != NULL) {
      //     other = (*value_ptr);
      //     return true;
      //   } else { return false; }
      // }
      
      // bool get_reference_unsafe(value_type*& ret) {
      //   ret = value_ptr;
      //   return ret != NULL;
      // }

      // /** returns true if set for the first time */
      // inline bool set_unsafe(lock_free_pool<value_type>& pool,
      //                        const value_type& other,
      //                        atomic<size_t> &joincounter) {
      //   if (value_ptr == NULL) {
      //     value_ptr = pool.alloc();
      //     (*value_ptr) = other;
      //     return true;
      //   } else {
      //     (*value_ptr) += other;
      //     joincounter.inc();
      //     return false;
      //   }
      // } // end of set_unsafe
      
      /** returns true if set for the first time */
      inline bool set(lock_free_pool<value_type>& pool,
                      const value_type& other,
                      value_type& new_value,
                      atomic<size_t>& joincounter) {
        bool ret = false;
        value_type toinsert = other;
        while(1) {
          value_type* vptr = VALUE_PENDING;
          // pull it out to process it
          atomic_exchange(value_ptr, vptr);
          // if there is nothing in there, set it
          // otherwise add it
          if (vptr == NULL) {
            vptr = pool.alloc();
            (*vptr) = toinsert;
            ret = true;
          } else if (vptr == VALUE_PENDING) {
            // a pending is in here. it is not ready for reading. try again.
            continue;
          } else { (*vptr) += toinsert; joincounter.inc(); }
          // swap it back in
          ASSERT_TRUE(vptr != VALUE_PENDING);
          new_value = *value_ptr;
          atomic_exchange(value_ptr, vptr);
          //aargh! I swapped something else out. Now we have to
          //try to put it back in
          if (__unlikely__(vptr != NULL && vptr != VALUE_PENDING)) {
            toinsert = (*vptr);
          } else { break; }
        }
        return ret;
      }
                
      void clear(lock_free_pool<value_type>& pool) {
        value_type val; test_and_get(val);
      }
      
      bool empty() { return value_ptr == NULL; }
      
      inline bool test_and_get(lock_free_pool<value_type>& pool,
                               value_type& r) {
        value_type* ret;
        while (1) {
          ret = value_ptr;
          if (ret == NULL) return false;
          else if (__likely__(ret != VALUE_PENDING)) {
            if (__likely__(atomic_compare_and_swap(value_ptr, ret, 
                                                   (value_type*)NULL))) {
              r = *ret;
              pool.free(ret);
              return true;
            }
          }
        }
        return false;
      } // end of test_and_get
    }; // end of atomic_box_type;


   
    typedef std::vector< atomic_box_type > atomic_box_vec_type; 
    atomic_box_vec_type atomic_box_vec;


    /** Not assignable */
    void operator=(const atomic_add_vector& other) { }


  public:
    /** Initialize the per vertex task set */
    atomic_add_vector(size_t num_vertices = 0) :
      pool(num_vertices + 256), atomic_box_vec(num_vertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      atomic_box_vec.resize(num_vertices);
      pool.reset_pool(num_vertices + 256);
    }

    /** Add a task to the set returning false if the task was already
        present. */
    bool add(const size_t& idx, 
             const value_type& val) {
      ASSERT_LT(idx, atomic_box_vec.size());
      value_type new_value;
      return atomic_box_vec[idx].set(pool, val, new_value, joincounter);
    } // end of add task to set 


    // /** Add a task to the set returning false if the task was already
    //     present. */
    // bool add_unsafe(const size_t& idx,
    //                 const value_type& val) {
    //   ASSERT_LT(idx, atomic_box_vec.size());
    //   return atomic_box_vec[idx].set_unsafe(pool, val, joincounter);
    // } // end of add task to set


    bool add(const size_t& idx, 
             const value_type& val,
             value_type& new_value) {
      ASSERT_LT(idx, atomic_box_vec.size());
      return atomic_box_vec[idx].set(pool, val, new_value, joincounter);
    } // end of add task to set 


    
    // /** Add a task to the set returning false if the task was already
    //     present. Returns the priority of the task before and after
    //     insertion. If the task did not exist prior to the add, 
    //     prev_priority = 0 */
    // bool add(const size_t& idx, 
    //          const value_type& val,
    //          double& prev_priority,
    //          double& new_priority) {
    //   ASSERT_LT(idx, atomic_box_vec.size());
    //   return atomic_box_vec[idx].set(pool, val, prev_priority, new_priority, 
    //                            joincounter);
    // } // end of add task to set 

    // bool get_nondestructive_unsafe(const size_t& idx,
    //                                value_type& ret_val) {
    //   return atomic_box_vec[idx].get_nondestructive_unsafe(ret_val);
    // }

    // bool get_reference_unsafe(const size_t& idx,
    //                           value_type*& ret_val) {
    //   return atomic_box_vec[idx].get_reference_unsafe(ret_val);
    // }


    bool test_and_get(const size_t& idx,
                      value_type& ret_val) {
      ASSERT_LT(idx, atomic_box_vec.size());
      return atomic_box_vec[idx].test_and_get(pool, ret_val);
    }
    
    bool empty(const size_t& idx) {
      return atomic_box_vec[idx].empty();
    }

    size_t size() const { 
      return atomic_box_vec.size(); 
    }
    
    size_t num_joins() const { 
      return joincounter.value;
    }
    
    
    void clear() {
      for (size_t i = 0; i < atomic_box_vec.size(); ++i) clear(i);
    }

    void clear(size_t i) { atomic_box_vec[i].clear(pool); }
    
  }; // end of vertex map

}; // end of namespace graphlab

#undef VALUE_PENDING

#endif

