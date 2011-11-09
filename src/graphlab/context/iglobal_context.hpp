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

#include <graphlab/util/generics/any.hpp>
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
      bool is_const = true;
      std::vector<spinlock>* locks_ptr;
      graphlab::any* values_ptr;
      get_global(key, is_const, locks_ptr, values_ptr);
      if(locks_ptr == NULL || values_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
      }
      if(is_const) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key 
          << "\" refers to a constant variable and cannot be set!" << std::endl;
      }      
      std::vector<spinlock>& locks = *locks_ptr;
      // Get the actual value [ This could generate a dynamic cast error]
      std::vector<T>& values = values_ptr->as< std::vector<T> > ();
      ASSERT_EQ(locks.size(), values.size());
      ASSERT_LT(index, values.size());
      // Update the value
      locks[index].lock(); values[index] = value; locks[index].unlock();
    } // end of set global

    /**
     * Atomically increment a global value
     */
    template<typename T>
    void increment_global(const std::string& key, const T& delta, 
                          size_t index = 0) {
      bool is_const = true;
      std::vector<spinlock>* locks_ptr;
      graphlab::any* values_ptr;
      get_global(key, is_const, locks_ptr, values_ptr);
      if(locks_ptr == NULL || values_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
      }
      if(is_const) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key 
          << "\" refers to a constant variable and cannot be set!" << std::endl;
      }      
      std::vector<spinlock>& locks = *locks_ptr;
      // Get the actual value [ This could generate a dynamic cast error]
      std::vector<T>& values = values_ptr->as< std::vector<T> > ();
      ASSERT_EQ(locks.size(), values.size());
      ASSERT_LT(index, values.size());
      // Update the value
      locks[index].lock(); values[index] += delta; locks[index].unlock();
    } // end of increment global

    /**
     * Atomically get a global value
     */
    template<typename T>
    void get_global(const std::string& key, T& ret_value, size_t index = 0) {
      bool is_const = true;
      std::vector<spinlock>* locks_ptr;
      graphlab::any* values_ptr;
      get_global(key, is_const, locks_ptr, values_ptr);
      if(locks_ptr == NULL || values_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
      }
      std::vector<spinlock>& locks = *locks_ptr;
      // Get the actual value [ This could generate a dynamic cast error]
      std::vector<T>& values = values_ptr->as< std::vector<T> > ();
      ASSERT_EQ(locks.size(), values.size());
      ASSERT_LT(index, values.size());
      // Get the value.  If it is constant we don't need the locks
      if(!is_const) locks[index].lock(); 
      ret_value = values[index]; 
      if(!is_const) locks[index].unlock();
    } // end of get global


    /**
     * Atomically get a global const
     */
    template<typename T>
    const T& get_global_const(const std::string& key, size_t index = 0) {
      bool is_const = true;
      std::vector<spinlock>* locks_ptr;
      graphlab::any* values_ptr;
      get_global(key, is_const, locks_ptr, values_ptr);
      if(locks_ptr == NULL || values_ptr == NULL) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" does not exist!" 
          << std::endl;
      }
      if(!is_const) {
        logstream(LOG_FATAL) 
          << "The global variable \"" << key << "\" is not a constant!" 
          << std::endl;
      }      
      std::vector<spinlock>& locks = *locks_ptr;
      // Get the actual value [ This could generate a dynamic cast error]
      std::vector<T>& values = values_ptr->as< std::vector<T> > ();
      ASSERT_LT(index, values.size());
      return values[index]; 
    } // end of get global const


    /**
     * Force the engine to stop executing additional update functions.
     */
    virtual void terminate() = 0;


  protected:

    /**
     * Get the internal information needed to access a global
     * variable.
     */
    virtual void get_global(const std::string& key,      
                            bool& ret_is_const,
                            std::vector<spinlock>*& ret_locks_ptr,
                            graphlab::any*& ret_values_ptr) = 0;



  }; // end of iglobal_context

  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

