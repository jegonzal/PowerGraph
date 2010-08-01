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
