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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#ifndef GRAPHLAB_SHARED_MEMORY_ENGINE_HPP
#define GRAPHLAB_SHARED_MEMORY_ENGINE_HPP

#include <cmath>
#include <cassert>
#include <algorithm>
#include <boost/bind.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/mutable_queue.hpp>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/shared_data/glshared.hpp>
#include <graphlab/engine/scope_manager_and_scheduler_wrapper.hpp>
#include <graphlab/metrics/metrics.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  
  /**
   * This class defines a basic shared_memory engine
   */
  template<typename Graph>
  class shared_memory_engine : public iengine<Graph> {
    
  public:

    // Include parent types
    typedef iengine<Graph> iengine_base;
    typedef typename iengine_base::vertex_id_type vertex_id_type;
    typedef typename iengine_base::update_task_type update_task_type;
    typedef typename iengine_base::update_function_type update_function_type;
    typedef typename iengine_base::ischeduler_type ischeduler_type;
    typedef typename iengine_base::imonitor_type imonitor_type;
    typedef typename iengine_base::termination_function_type termination_function_type;
    typedef typename iengine_base::iscope_type iscope_type;
    typedef typename iengine_base::sync_function_type sync_function_type;
    typedef typename iengine_base::merge_function_type merge_function_type;

    
    typedef general_scope_factor<Graph> scope_factory_type;
    
  private:

    //! The Graph that the engine is operating on
    graph_type* graph;
   
    /** 
     * The number vertices in the graph.  This number is set when
     * internal data structures are initialized and is used to track
     * whether the graph has changed externally since the engine was
     * last executed.
     */
    size_t nverts;
    
    //! The active scope factory
    scope_factory_type* scope_factory_ptr;

    //! The active scheduler
    ischeduler_type* scheduler_ptr;
    
    //! The local engine options
    engine_options eopts;
   
  public:
   
    //! Create an engine for the given graph
    shared_memory_engine(Graph& graph);
    
    //! Clear internal members
    ~shared_memory_engine();
    
    //! Start the engine
    void start();
    
    //! Stop the engine
    void stop();
    
    //! \brief Describe the reason for termination.
    exec_status last_exec_status();

    //! \brief Get the number of updates executed by the engine.
    size_t last_update_count();

        
    //! \brief Register a monitor with an engine. 
    void register_monitor(imonitor_type* listener);
    
    //! \brief Adds an update task with a particular priority.
    void add_task(update_task_type task, double priority);

    //! \brief Add a vector of tasks

    void add_task(const std::vector<vertex_id_type>& vertices,
                  update_function_type func, double priority);

 
    //! \brief Apply update function to all the vertices in the graph
    void add_task_to_all(update_function_type func,
                         double priority);

    //! \brief associate a termination function with this engine.
    void add_terminator(termination_function_type term) = 0;

    //!  remove all associated termination functions
    void clear_terminators();
    
    //! \brief The timeout is the total
    void set_timeout(size_t timeout_secs);

    
    //! \brief set a limit on the number of tasks that may be executed.
    void set_task_budget(size_t max_tasks);

    //! \brief Update the engine options.  
    void set_engine_options(const engine_options& opts);


    //! \brief Registers a sync with the engine.
    void set_sync(glshared_base& shared,
                  sync_function_type sync,
                  glshared_base::apply_function_type apply,
                  const any& zero,
                  size_t sync_interval = 0,
                  merge_function_type merge = NULL,
                  vertex_id_type rangelow = 0,
                  vertex_id_type rangehigh = -1);


    //! Performs a sync immediately.
    virtual void sync_now(glshared_base& shared) = 0;
    



    ///////////////////////////////////////////////////////////////
    /// New Functions


    /** 
     * This function clears an internal engine state other than the
     * graph and the engine options.
     */
    void reset();

        
  }; // end of shared_memory engine




  /////////////////////////////////////////////////////////////////////////
  /// Implementation

  template<typename Graph> 
  shared_memory_engine<Graph>::
  shared_memory_engine(Graph& graph) : 
    graph(&graph), 
    nvert(graph.num_vertices()),
    scope_factory_ptr(NULL),
    scheduler_ptr(NULL) { } // end of constructor


  template<typename Graph> 
  shared_memory_engine<Graph>::
  ~shared_memory_engine(Graph& graph) {
    clear_members();
  } // end of destructor


  template<typename Graph>
  void
  shared_memory_engine<Graph>::
  set_engine_options(const engine_options& new_opts) {
    clear_members();
    eopts = new_opts;
  } // end of set_engine_options

  template<typename Graph>
  void
  shared_memory_engine<Graph>::reset() {
    ASSERT_TRUE((scope_factory_ptr != NULL && scheduler_ptr != NULL) ||
                (scope_factor_ptr == NULL && scheduler_ptr == NULL));
    if(scope_factory_ptr != NULL) {
      delete scope_factory_ptr;
      scope_factory_ptr = NULL;
    }
    if(scheduler_ptr != NULL) {
      delete scheduler_ptr; 
      scheduler_ptr = NULL;
    }
  } // end of clear_members



  

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

