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



#ifndef GRAPHLAB_VERTEX_FUNCTOR_SET_HPP
#define GRAPHLAB_VERTEX_FUNCTOR_SET_HPP


#include <vector>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/lock_free_pool.hpp>


#define UPDATE_FUNCTOR_PENDING (update_functor_type*)(size_t)(-1)
namespace graphlab {

  template<typename UpdateFunctor>
  class vertex_functor_set {
  public:
    typedef UpdateFunctor update_functor_type;

    
  private:
    lock_free_pool<update_functor_type> pool;
    class vfun_type {
    private:
      update_functor_type* functor;

    public:
      vfun_type() : functor(NULL) { }
      
      void assign_unsync(const vfun_type &other) {
        functor = other.functor;
      }
      /** returns true if set for the first time */
      inline bool set(lock_free_pool<update_functor_type>& pool,
                      const update_functor_type& other) {
        bool ret = false;
        update_functor_type toinsert = other;
        while(1) {
          update_functor_type* uf = UPDATE_FUNCTOR_PENDING;
          // pull it out to process it
          atomic_exchange(uf, functor);
          // if there is nothing in there, set it
          // otherwise add it
          if (uf == NULL) {
            uf = pool.alloc();
            (*uf) = toinsert;
            ret = true;
          }
          else if (uf == UPDATE_FUNCTOR_PENDING) {
            // a pending is in here. it is not ready for reading. try again.
            continue;
          }
          else {
            (*uf) += toinsert;
          }
          // swap it back in
          atomic_exchange(uf, functor);
          //aargh! I swapped something else out. Now we have to
          //try to put it back in
          if (uf != NULL && uf != UPDATE_FUNCTOR_PENDING) {
            toinsert = (*uf);
          }
          else {
            break;
          }
        }
        return ret;
      }

      /** returns true if set for the first time */
      inline bool set(lock_free_pool<update_functor_type>& pool,
                      const update_functor_type& other, 
                      double& ret_priority) {
        bool ret = false;
        update_functor_type toinsert = other;
        while(1) {
          update_functor_type* uf = UPDATE_FUNCTOR_PENDING;
          // pull it out to process it
          atomic_exchange(uf, functor);
          // if there is nothing in there, set it
          // otherwise add it
          if (uf == NULL) {
            uf = pool.alloc();
            (*uf) = toinsert;
            ret = true;
          }
          else if (uf == UPDATE_FUNCTOR_PENDING) {
            // a pending is in here. it is not ready for reading. try again.
            continue;
          }
          else {
            (*uf) += toinsert;
          }
          ret_priority = uf->priority();
          // swap it back in
          atomic_exchange(uf, functor);
          //aargh! I swapped something else out. Now we have to
          //try to put it back in
          if (uf != NULL && uf != UPDATE_FUNCTOR_PENDING) {
            toinsert = (*uf);
          }
        }
        return ret;
      }

      inline update_functor_type get(lock_free_pool<update_functor_type>& pool) {
        update_functor_type* ret;
        while (1) {
          ret = functor;
          if (ret == NULL) return update_functor_type();
          else if (ret != UPDATE_FUNCTOR_PENDING) {
            if (atomic_compare_and_swap(functor, ret, (update_functor_type*)NULL)) {
              update_functor_type r = *ret;
              pool.free(ret);
              return r;
            }
          }
        }
        return update_functor_type();
      }
            
      void reset_unsync() {
        functor = NULL;
      }
      
      bool has_task() {
        return functor != NULL;
      }
      
      inline bool test_and_get(lock_free_pool<update_functor_type>& pool,
                               update_functor_type& r) {
        update_functor_type* ret;
        while (1) {
          ret = functor;
          if (ret == NULL) return false;
          else if (ret != UPDATE_FUNCTOR_PENDING) {
            if (atomic_compare_and_swap(functor, ret, (update_functor_type*)NULL)) {
              r = *ret;
              pool.free(ret);
              return true;
            }
          }
        }
        return false;
      }
    }; // end of vfun_type;

   
    typedef std::vector< vfun_type > vfun_set_type; 
    vfun_set_type vfun_set;


  public:
    /** Initialize the per vertex task set */
    vertex_functor_set(size_t num_vertices = 0) :
      pool(num_vertices + 256), vfun_set(num_vertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      vfun_set.resize(num_vertices);
    }

    void operator=(const vertex_functor_set& other) {
      resize(other.vfun_set.size());
      for (size_t i = 0;i < vfun_set.size(); ++i) {
        vfun_set[i].assign_unsync(other.vfun_set[i]);
      }
    }

    
    /** Add a task to the set returning false if the task was already
        present. Promote task to max(old priority, new priority) */
    bool add(const vertex_id_type& vid, 
             const update_functor_type& fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(pool, fun);
    } // end of add task to set 

    
    /** Add a task to the set returning false if the task was already
        present. Promote task to max(old priority, new priority) */
    bool add(const vertex_id_type& vid, 
             const update_functor_type& fun,
             double& ret_priority) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(pool, fun, ret_priority);
    } // end of add task to set 



    bool test_and_get(const vertex_id_type& vid,
                      update_functor_type& ret_fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].test_and_get(pool, ret_fun);
    }
    
    bool has_task(const vertex_id_type& vid) {
      return vfun_set[vid].has_task();
    }

    size_t size() const { 
      return vfun_set.size(); 
    }
    
    void clear_unsync() {
      for (size_t i = 0; i < vfun_set.size(); ++i) vfun_set[i].reset_unsync();
    }
    
  }; // end of vertex functor set

}; // end of namespace graphlab

#undef UPDATE_FUNCTOR_PENDING

#endif

