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



#ifndef GRAPHLAB_ATOMIC_ADD_VECTOR2_HPP
#define GRAPHLAB_ATOMIC_ADD_VECTOR2_HPP


#include <vector>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/lock_free_pool.hpp>



namespace graphlab {

  /**
   * \TODO DOCUMENT THIS CLASS
   */ 
  
  template<typename ValueType>
  class atomic_add_vector2 {
  public:
    typedef ValueType value_type;

    // We needed a second "NULL" pointer value to indicate a value
    // that is being swapped.
#define VALUE_PENDING (value_type*)(size_t)(-1)
    
  private:
    atomic<size_t> joincounter;
    
    class atomic_box_type {
    private:
      simple_spinlock lock;
      bool _empty;
      value_type value;
    public:
      atomic_box_type() : _empty(true) { }
      /** returns true if set for the first time */
      inline bool set(const value_type& other,
                      value_type& new_value,
                      atomic<size_t>& joincounter) {
        bool first_set = false;
        lock.lock();
        if(!_empty) value += other; 
        else { value = other; first_set = true; }
        new_value = value;
        _empty = false;
        lock.unlock();
        return first_set;
      }
                
      void clear() {
        value_type val; test_and_get(val);
      }
      
      bool empty() { return _empty; }
      
      inline bool test_and_get(value_type& ret_val) {
        bool success = false;
        lock.lock();
        if(!_empty) {
          success = true;
          ret_val = value;
          _empty = true;
        }
        lock.unlock();
        return success;
      } // end of test_and_get

    }; // end of atomic_box_type;


   
    typedef std::vector< atomic_box_type > atomic_box_vec_type; 
    atomic_box_vec_type atomic_box_vec;


    /** Not assignable */
    void operator=(const atomic_add_vector2& other) { }


  public:
    /** Initialize the per vertex task set */
    atomic_add_vector2(size_t num_vertices = 0) :
      atomic_box_vec(num_vertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      atomic_box_vec.resize(num_vertices);
    }

    /** Add a task to the set returning false if the task was already
        present. */
    bool add(const size_t& idx, 
             const value_type& val) {
      ASSERT_LT(idx, atomic_box_vec.size());
      value_type new_value;
      return atomic_box_vec[idx].set( val, new_value, joincounter);
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
      return atomic_box_vec[idx].set(val, new_value, joincounter);
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
    //   return atomic_box_vec[idx].set( val, prev_priority, new_priority, 
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
      return atomic_box_vec[idx].test_and_get( ret_val);
    }

    // bool test_peek(const size_t& idx,
    //                value_type& ret_val) {
    //   ASSERT_LT(idx, atomic_box_vec.size());
    //   return atomic_box_vec[idx].test_and_peek(pool, ret_val);
    // }
    
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

    void clear(size_t i) { atomic_box_vec[i].clear(); }
    
  }; // end of vertex map

}; // end of namespace graphlab

#undef VALUE_PENDING

#endif

