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



#ifndef GRAPHLAB_VERTEX_MAP_HPP
#define GRAPHLAB_VERTEX_MAP_HPP


#include <vector>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/lock_free_pool.hpp>


#define VALUE_PENDING (value_type*)(size_t)(-1)

namespace graphlab {

  /**
   * \TODO DOCUMENT THIS CLASS
   */ 
  
  template<typename ValueType>
  class vertex_map {
  public:
    typedef ValueType value_type;

    
  private:
    lock_free_pool<value_type> pool;
    atomic<size_t> joincounter;
    
    class vfun_type {
    private:
      value_type* volatile value_ptr;
    public:
      vfun_type() : value_ptr(NULL) {  }
      
      void assign_unsync(const vfun_type &other) {
        value_ptr = other.value_ptr;
      }
      
      bool get_nondestructive_unsync(const value_type& other) {
        if (value_ptr != NULL) {
          other = (*value_ptr);
          return true;
        } else { return false; }
      } // end of get_nondestructive_unsync
      
      bool get_reference_unsync(value_type*& ret) {
        ret = value_ptr;
        return ret != NULL;
      }

      /** returns true if set for the first time */
      inline bool set_unsafe(lock_free_pool<value_type>& pool,
                             const value_type& other,
                             atomic<size_t> &joincounter) {
        if (value_ptr == NULL) {
          value_ptr = pool.alloc();
          (*value_ptr) = other;
          return true;
        } else {
          (*value_ptr) += other;
          joincounter.inc();
          return false;
        }
      } // end of set_unsafe
      
      /** returns true if set for the first time */
      inline bool set(lock_free_pool<value_type>& pool,
                      const value_type& other,
                      atomic<size_t> &joincounter) {
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
          atomic_exchange(value_ptr, vptr);
          //aargh! I swapped something else out. Now we have to
          //try to put it back in
          if (__unlikely__(vptr != NULL && vptr != VALUE_PENDING)) {
            toinsert = (*vptr);
          } else { break; }
        }
        return ret;
      }


      /** returns true if set for the first time */
      inline bool merge(lock_free_pool<value_type>& pool,
                        const value_type& other) {
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
          } else { (*vptr).merge(toinsert); }
          // swap it back in
          ASSERT_TRUE(vptr != VALUE_PENDING);
          atomic_exchange(value_ptr, vptr);
          //aargh! I swapped something else out. Now we have to
          //try to put it back in
          if (__unlikely__(vptr != NULL && vptr != VALUE_PENDING)) {
            toinsert = (*vptr);
          }
          else {
            break;
          }
        }
        return ret;
      }
      
      /** returns true if set for the first time */
      inline bool set(lock_free_pool<value_type>& pool,
                      const value_type& other, 
                      double& prev_priority, 
                      double& new_priority,
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
            prev_priority = 0;
            vptr = pool.alloc();
            (*vptr) = toinsert;
            ret = true;
          } else if (vptr == VALUE_PENDING) {
            // a pending is in here. it is not ready for reading. try again.
            continue;
          } else {
            prev_priority = vptr->priority();
            (*vptr) += toinsert;
            joincounter.inc();
          }
          new_priority = vptr->priority();
          // swap it back in
          ASSERT_TRUE(vptr != VALUE_PENDING);
          atomic_exchange(value_ptr, vptr);
          //aargh! I swapped something else out. Now we have to
          //try to put it back in
          if (__unlikely__(vptr != NULL && vptr != VALUE_PENDING)) {
            toinsert = (*vptr);
          } else { break; }
        }
        return ret;
      }

      inline value_type get(lock_free_pool<value_type>& pool) {
        value_type* ret;
        while (1) {
          ret = value_ptr;
          if (ret == NULL) return value_type();
          else if (ret != VALUE_PENDING) {
            if (__likely__(atomic_compare_and_swap(value_ptr, ret, 
                                                   (value_type*)NULL))) {
              value_type r = *ret;
              pool.free(ret);
              return r;
            }
          }
        }
        return value_type();
      }
            
      void reset_unsync() {
        value_ptr = NULL;
      }
      
      bool has_task() {
        return value_ptr != NULL;
      }
      
      inline bool test_and_get(lock_free_pool<value_type>& pool,
                               value_type& r) {
        value_type* ret;
        while (1) {
          ret = value_ptr;
          if (__unlikely__(ret == NULL)) return false;
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

    }; // end of vfun_type;


   
    typedef std::vector< vfun_type > vfun_set_type; 
    vfun_set_type vfun_set;


  public:
    /** Initialize the per vertex task set */
    vertex_map(size_t num_vertices = 0) :
      pool(num_vertices + 256), vfun_set(num_vertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      vfun_set.resize(num_vertices);
      pool.reset_pool(num_vertices + 256);
    }

    void operator=(const vertex_map& other) {
      resize(other.vfun_set.size());
      for (size_t i = 0;i < vfun_set.size(); ++i) {
        vfun_set[i].assign_unsync(other.vfun_set[i]);
      }
    }

    
    /** Add a task to the set returning false if the task was already
        present. */
    bool add(const vertex_id_type& vid, 
             const value_type& fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(pool, fun, joincounter);
    } // end of add task to set 

    /** Add a task to the set returning false if the task was already
        present. */
    bool add_unsafe(const vertex_id_type& vid,
                    const value_type& fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set_unsafe(pool, fun, joincounter);
    } // end of add task to set

    /** Add a task to the set returning false if the task was already
        present. */
    bool merge(const vertex_id_type& vid, 
               const value_type& fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].merge(pool, fun);
    } // end of add task to set 


    /** Add a task to the set returning false if the task was already
        present. Also returns the combined priority of the task. */
    bool add(const vertex_id_type& vid, 
             const value_type& fun,
             double& new_priority) {
      ASSERT_LT(vid, vfun_set.size());
      double unused = 0;
      return vfun_set[vid].set(pool, fun, unused, new_priority, joincounter);
    } // end of add task to set 


    
    /** Add a task to the set returning false if the task was already
        present. Returns the priority of the task before and after
        insertion. If the task did not exist prior to the add, 
        prev_priority = 0 */
    bool add(const vertex_id_type& vid, 
             const value_type& fun,
             double& prev_priority,
             double& new_priority) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(pool, fun, prev_priority, new_priority, 
                               joincounter);
    } // end of add task to set 

    bool get_nondestructive_unsync(const vertex_id_type& vid,
                                   value_type& ret_fun) {
      return vfun_set[vid].get_nondestructive_unsync(ret_fun);
    }

    bool get_reference_unsync(const vertex_id_type& vid,
                              value_type*& ret_fun) {
      return vfun_set[vid].get_reference_unsync(ret_fun);
    }


    bool test_and_get(const vertex_id_type& vid,
                      value_type& ret_fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].test_and_get(pool, ret_fun);
    }
    
    bool has_task(const vertex_id_type& vid) {
      return vfun_set[vid].has_task();
    }

    size_t size() const { 
      return vfun_set.size(); 
    }
    
    size_t num_joins() const { 
      return joincounter.value;
    }
    
    void clear_unsync() {
      for (size_t i = 0; i < vfun_set.size(); ++i) vfun_set[i].reset_unsync();
    }
    
  }; // end of vertex map

}; // end of namespace graphlab

#undef VALUE_PENDING

#endif

