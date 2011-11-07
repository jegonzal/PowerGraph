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

    virtual size_t num_vertices() const = 0;
    virtual size_t num_edges() const = 0;
    
    template<typename T>
    void set_global(const std::string& key, const T& value, size_t index = 0) {
      std::pair<std::vector<spinlock>*, any*> pair = get_global_pair(key);
      std::vector<spinlock>& locks = *pair.first;
      // Get the actual value [ This could generate a dynamic cast error]
      std::vector<T>& values = pair.second->as< std::vector<T> > ();
      ASSERT_EQ(locks.size(), values.size());
      ASSERT_LT(index, values.size());
      // Update the value
      locks[index].lock(); values[index] = value; locks[index].unlock();
    }

    template<typename T>
    void get_global(const std::string& key, T& ret_value, size_t index = 0) {
      std::pair<std::vector<spinlock>*, any*> pair = get_global_pair(key);
      const std::vector<spinlock>& locks = *pair.first;
      // Get the actual value [ This could generate a dynamic cast error]
      const std::vector<T>& values = pair.second->as< std::vector<T> > ();
      ASSERT_EQ(locks.size(), values.size());
      ASSERT_LT(index, values.size());
      // Update the value
      locks[index].lock(); ret_value = values[index]; locks[index].unlock();
    }

    /**
     * Calling this function will force the engine to abort
     * immediately
     */
    virtual void terminate() = 0;


  protected:

    virtual std::pair<std::vector<spinlock>*, any*> get_global_pair(const std::string& key) = 0;
  }; // end of iglobal_context

  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

