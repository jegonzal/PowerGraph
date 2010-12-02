/**
 * \file 
 *
 * This file contains the shared data interface which is used to
 * represent the global shared data within a parallel execution of the
 * graphlab program.
 *
 */

#ifndef GRAPHLAB_ISHARED_DATA_HPP
#define GRAPHLAB_ISHARED_DATA_HPP

#include <cmath>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/scope/iscope_factory.hpp>
#include <graphlab/shared_data/glshared.hpp>

namespace graphlab {



  /**
   * \brief The shared data interface
   *
   * Often in a graphlab program there is data (constant and mutable)
   * that cannot be directly expressed in the graph.  Typical forms of
   * shared data may be counters or global statistics that evolve over
   * time.  The ishared_data class describes the interface given to
   * each update function that can be used to read and modify the
   * global shared data.  
   * 
   */
  template<typename Graph>
  class ishared_data {
  public:
    typedef Graph         graph_type;
    typedef typename Graph::vertex_data_type vertex_data_type;
    typedef iscope<Graph> iscope_type;
    typedef iscope_factory<Graph> iscope_factory_type;

    typedef void(*sync_function_type)(size_t index,
                                      const ishared_data& shared_data,
                                      iscope_type& scope,
                                      any& accumulator);
    
    typedef void(*apply_function_type)(size_t index,
                                       const ishared_data& shared_data,
                                       any& current_data,
                                       const any& new_data);

    typedef void(*merge_function_type)(size_t index,
                                       const ishared_data& shared_data,
                                       any& merge_dest,
                                       const any& merge_src);

    
  public:

    // Virtual Functions
    // ======================================================>    

    /** This is used to manage constants */
    virtual const any& get_constant(size_t index) const = 0;
        
    /** Return the shared data at the specific location.
        In the distributed setting. This value could be slightly outdated*/
    virtual any get(size_t index) const = 0;


    /** Return the shared data at the specific location */
    virtual any atomic_get(size_t index) const = 0;
    
    /** Atomically sets the shared data at the specific location */
    virtual void atomic_set(size_t index, const any& data) = 0;
   
    /** Atomically sets the shared data at the specific location */
    virtual any atomic_exchange(size_t index, const any& data) = 0;


    /** Atomically applies the function to the shared data with the
        provided closure.  The application is done eventually */
    virtual void atomic_apply(size_t index,
                              apply_function_type fun,     
                              const any& closure) = 0;

    virtual void trigger_sync(size_t index) = 0;
    virtual void trigger_sync_all() = 0;
    
    /// TODO: backward compatibility hack
    virtual void trigger_sync(glshared_base &var) { }
  };



}; // end of namespace graphlab


#endif
