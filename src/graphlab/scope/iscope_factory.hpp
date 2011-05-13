/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
   \author Yucheng Low (ylow), Joseph Gonzalez
   Interface for Scope Factories
*/
#ifndef GRAPHLAB_ISCOPE_FACTORY_HPP
#define GRAPHLAB_ISCOPE_FACTORY_HPP

#include <graphlab/scope/iscope.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {



  template<typename Graph>
  class iscope_factory {
  public:


    
    typedef Graph graph_type;
    typedef iscope<Graph> iscope_type;

    /**  \note This constructor here does not actually do anything. It just exists
         to force the derived class constructors to look like this     */
    iscope_factory(Graph& graph, size_t ncpus) {}
 
    virtual ~iscope_factory() {}

    //!  Returns a scope around a particular vertex
    virtual iscope_type* get_scope(size_t cpuid,
                                   vertex_id_t vertex,
                                   scope_range::scope_range_enum s = scope_range::USE_DEFAULT) = 0;

    //! Set the default scope type
    virtual void set_default_scope(scope_range::scope_range_enum default_scope_range) = 0;
  
    //! Destroys a scope
    virtual void release_scope(iscope<Graph>* scope) = 0;

    //! Get the number of vertices
    virtual size_t num_vertices() const = 0;

  };

}

#include <graphlab/macros_undef.hpp>
#endif

