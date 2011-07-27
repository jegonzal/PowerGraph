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
 * \todo: Evaluate if we still really need this interface.
 */

#ifndef GRAPHLAB_ISCOPE_MANAGER_HPP
#define GRAPHLAB_ISCOPE_MANAGER_HPP

#include <graphlab/scope/consistency_model.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {



  template<typename Graph>
  class iscope_manager {
  public:


    
    typedef Graph graph_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef iscope<Graph> iscope_type;

    /**  \note This constructor here does not actually do anything. It
         just exists to force the derived class constructors to look
         like this  */
    iscope_manager(Graph& graph, size_t ncpus) {}
 
    virtual ~iscope_manager() {}

    //!  Returns a scope around a particular vertex
    virtual iscope_type& get_scope(size_t cpuid,
                                   vertex_id_type vertex,
                                   consistency_model::model_enum s = 
                                   consistency_model::USE_DEFAULT) = 0;

    //! Set the default scope type
    virtual void set_default_scope(consistency_model::model_enum
                                   default_scope_range) = 0;
  
    //! Destroys a scope
    virtual void release_scope(iscope<graph_type>& scope) = 0;

    //! Get the number of vertices
    virtual size_t num_vertices() const = 0;

  };

}

#include <graphlab/macros_undef.hpp>
#endif

