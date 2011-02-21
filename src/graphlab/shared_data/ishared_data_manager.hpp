/**
  * \file Interface for shared data manager implementations.
  **/

#ifndef GRAPHLAB_ISHARED_DATA_MANAGER
#define GRAPHLAB_ISHARED_DATA_MANAGER

#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/scope/iscope_factory.hpp>


namespace graphlab {

   /** \brief \deprecated Use glshared */
  template<typename Graph>
  class ishared_data_manager : 
    public ishared_data<Graph> {

  public:
    typedef Graph graph_type;
    typedef ishared_data<Graph> base;

    typedef typename base::iscope_type iscope_type;
    typedef typename base::iscope_factory_type iscope_factory_type;
    
    typedef typename base::sync_function_type sync_function_type;
    typedef typename base::apply_function_type apply_function_type;


    /**
     * \brief run sync on a particular field
     *
     * User code should not call this function and instead use the
     * sync which takes a graph argument.
     *
     * \todo this should be optimized
     */
    virtual void sync(size_t index) = 0;
    
    /**
     * \brief Run the particular sync task using the graph data.
     *
     * This function should not be called by user code outside of the
     * graphlab engine.
     *
     * \todo this should be optimized
     *
     */
    virtual void sync_all() = 0;

    
    /**
     * \brief Run the particular sync task using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine.
     *
     */
    virtual void sync(Graph& graph, size_t index) = 0;


    /**
     * \brief Run all sync tasks using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine.
     *
     */

    virtual void sync_all(Graph& graph) = 0;
    
    

    /** Notify that a sync may be needed on index */
    virtual void signal(size_t index) = 0;
    virtual void signal_all() = 0;

    
    /** syncs as soon as possible */
    virtual void trigger_sync(size_t index) = 0;
    virtual void trigger_sync_all() = 0;

    /* Note: only used by pushy engine (Added by Aapo 5/31/10) */
    virtual void progress(size_t cpuid, iscope_type * scope) {}
    
    /**
     * \brief register a sync.
     *
     * If the sync_interval is 0 then the sync is never run in the
     * background.
     * vertices in the range [rangelow, rangehigh) will be executed.
     * (including rangelow, not including rangehigh)
     */
    virtual void set_sync(size_t index,
                          sync_function_type sync_fun,
                          apply_function_type apply_fun,
                          const any& zero,
                          size_t sync_interval,
                          size_t rangelow = 0,
                          size_t rangehigh = -1) = 0;
    
    virtual void create_atomic(size_t index, const any& initial_value) = 0;


    /** Set the scpe factory, called by the engine */
    virtual void set_scope_factory(iscope_factory_type* factory) = 0;
    
    /** Creates a constant value*/
    virtual void set_constant(size_t index, const any& new_value) = 0;


  }; // end of class ishared data manager


}; // end of graphlab namespace


#endif 
