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

#include <graphlab/parallel/pthread_tools.hpp>



namespace graphlab {

  template<typename Engine>
  class vertex_functor_set {
  public:
    typedef Engine engine_type;
    typedef typename engine_type::vertex_id_type vertex_id_type;
    typedef typename engine_type::update_functor_type update_functor_type;


  private:
    
    class vfun_type {
    private:
      mutex lock;
      bool is_set;
      update_functor_type functor;
    public:
      vfun_type() : is_set(false) { }
      inline bool set(const update_functor_type& other) {
        lock.lock();
        const bool already_set(is_set);
        if(is_set) functor += other;
        else functor = other;
        lock.unlock();
        return already_set;
      }
      inline update_functor_type get() {
        update_functor_type ret;
        lock.lock();
        ASSERT_TRUE(is_set);
        ret = functor;
        is_set = false;
        lock.unlock();
        return functor;
      }
      inline bool test_and_get(update_functor_type& ret) {
        lock.lock();
        const bool success(is_set);
        if(success) {
          ret = functor;
          is_set = false;
        }
        lock.unlock();
        return success;
      }
    }; // end of vfun_type;

   
    typedef std::vector< vfun_type > vfun_set_type; 
    vfun_set_type vfun_set;


  public:
    /** Initialize the per vertex task set */
    vertex_functor_set(size_t num_vertices) :
      vfun_set(num_vertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      vfun_set.resize(num_vertices);
    }

    
    /** Add a task to the set returning false if the task was already
        present. Promote task to max(old priority, new priority) */
    bool add(const vertex_id_type& vid, 
             const update_functor_type& fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(fun);
    } // end of add task to set 


    bool test_and_get(const vertex_id_type& vid,
                      update_functor_type& ret_fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].test_and_get(ret_fun);
    }
    
    size_t size() const { 
      return vfun_set.size(); 
    }
    
  }; // end of vertex functor set

}; // end of namespace graphlab


#endif

