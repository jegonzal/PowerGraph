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


/** \file 
 *
 * This file describes the icontext interface as well as the the
 * context_range_enum.
 *
 */
#ifndef GRAPHLAB_IGLOBAL_CONTEXT_HPP
#define GRAPHLAB_IGLOBAL_CONTEXT_HPP

#include <string>

#include <graphlab/util/generics/any_vector.hpp>
#include <graphlab/parallel/pthread_tools.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {



  class iglobal_context {
  public:
    virtual ~iglobal_context() { }

    /**
     * Get the number of vertices in the graph.
     */
    virtual size_t num_vertices() const = 0;

    /**
     * Get the number of edges in the graph
     */
    virtual size_t num_edges() const = 0;

    /**
     * Get an estimate of the number of update functions executed up
     * to this point.
     */
    virtual size_t num_updates() const = 0;
    
    /**
     * Atomically set a global value
     */
    template<typename T>
    void set_global(const std::string& key, const T& value, size_t index = 0) {
      graphlab::any_vector* any_vec_ptr(NULL);
      bool is_const = true;
      get_global(key, any_vec_ptr, is_const);
      if(any_vec_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
        return;
      }
      if(is_const) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key 
          << "\" refers to a constant variable and cannot be set!" << std::endl;
        return;
      }      
      ASSERT_LT(index, any_vec_ptr->size());
      acquire_lock(key, index);
      // Get the actual value [ This could generate a dynamic cast error]
      any_vec_ptr->as<T>(index) = value;
      commit_change(key, index);
      release_lock(key, index);
    } // end of set global

    /**
     * Atomically increment a global value
     */
    template<typename T>
    void increment_global(const std::string& key, const T& delta, 
                          size_t index = 0) {
      graphlab::any_vector* any_vec_ptr(NULL);
      bool is_const = true;
      get_global(key, any_vec_ptr, is_const);
      if(any_vec_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
        return;
      }
      if(is_const) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key 
          << "\" refers to a constant variable and cannot be set!" << std::endl;
        return;
      }      
      ASSERT_LT(index, any_vec_ptr->size());
      acquire_lock(key, index);
      // Get the actual value [ This could generate a dynamic cast error]
      any_vec_ptr->as<T>(index) += delta;
      commit_change(key,index);
      release_lock(key, index);
    } // end of increment global

    /**
     * Atomically get a global value
     */
    template<typename T>
    void get_global(const std::string& key, T& ret_value, size_t index = 0) {
      graphlab::any_vector* any_vec_ptr(NULL);
      bool is_const = true;
      get_global(key, any_vec_ptr, is_const);
      if(any_vec_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
        return;
      }
      ASSERT_LT(index, any_vec_ptr->size());
      acquire_lock(key, index);
      // Get the actual value [ This could generate a dynamic cast error]
      ret_value = any_vec_ptr->as<T>(index);
      release_lock(key, index);
    } // end of get global


    /**
     * Atomically get a global const
     */
    template<typename T>
    const T& get_global_const(const std::string& key, size_t index = 0) {
      graphlab::any_vector* any_vec_ptr(NULL);
      bool is_const = true;
      get_global(key, any_vec_ptr, is_const);
      if(any_vec_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
        return;
      }
      if(!is_const) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" is not a constant!" 
          << std::endl;
      }      
      ASSERT_LT(index, any_vec_ptr->size());
      // Get the actual value [ This could generate a dynamic cast error]
      return any_vec_ptr->as<T>(index);
    } // end of get global const


    /**
     * Force the engine to stop executing additional update functions.
     */
    virtual void terminate() = 0;


  protected:

    /**
     * Acquire the lock on a particular entry in a global table. 
     */
    virtual void acquire_lock(const std::string& key, size_t index = 0) = 0;

    /**
     * Get the global table for a particular key.  This will be called
     * before acquiring the lock but no entires will be touched until
     * the lock has been acquired.  
     *
     * If the entry is not present in the table then the ret_vec_ptr
     * should be set to NULL. 
     *
     * Some tables are strictly constant an therefore ret_is_const
     * should return true in these cases.  If an entry is strictly
     * const then it cannot be set and users can take references to
     * that entry.
     */
    virtual void get_global(const std::string& key, 
                            any_vector*& ret_vec_ptr,
                            bool& ret_is_const) = 0;

    /**
     * Commit change is called if the entry at the location index has
     * been modified.  This is called while holding the lock on that
     * location.
     */
    virtual void commit_change(const std::string& key, size_t index = 0) = 0;

    /**
     * Release the lock on a particular entry in a global table.
     */
    virtual void release_lock(const std::string& key, size_t index = 0) = 0;



  }; // end of iglobal_context

  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

