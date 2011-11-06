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
    void set_global(const std::string& key, const T& value) {
      std::pair<graphlab::mutex*, graphlab::any*> pair = get_any_pair(key);
      pair->first->lock();
      pair->second->as<T>() = value;
      pair->first->unlock();
    }

    template<typename T>
    void get_global(const std::string& key, T& ret_value) const {
      std::pair<graphlab::mutex*, graphlab::any*> pair = get_any_pair(key);
      pair->first->lock();
      ret_value = pair->second->as<T>(); 
      pair->first->unlock();
    }

    /**
     * Calling this function will force the engine to abort
     * immediately
     */
    virtual void terminate() = 0;


  protected:

    virtual std::pair<mutex*, any*>  get_any_pair(const std::string& key) = 0;
  }; // end of iglobal_context

  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

