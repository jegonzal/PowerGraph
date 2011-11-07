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


#ifndef GRAPHLAB_GLOBAL_CONTEXT_HPP
#define GRAPHLAB_GLOBAL_CONTEXT_HPP

#include <string>

#include <graphlab/util/generics/any.hpp>
#include <graphlab/parallel/pthread_tools.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename Engine>
  class global_context : public iglobal_context {
  public:
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
  protected:
    engine_type* engine_ptr;
    graph_type* graph_ptr;
  public:
    global_context(engine_type* engine_ptr = NULL, 
                   graph_type* graph_ptr = NULL) : 
      engine_ptr(engine_ptr), graph_ptr(graph_ptr) { }
    size_t num_vertices() const { return graph_ptr->num_vertices(); }
    size_t num_edges() const { return graph_ptr->num_edges(); }
    void terminate() { engine_ptr->stop(); }

  protected:
    std::pair<mutex*, any*> get_any_pair(const std::string& key) {
      return engine_ptr->get_global_pair(key);
    } // end of get_any_pair


  }; // end of global_context

  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

