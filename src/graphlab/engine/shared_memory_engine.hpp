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


#include <graphlab/logger/logger.hpp>

#include <graphlab/options/graphlab_options.hpp>

#include <graphlab/graph/graph.hpp>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/scope/general_scope_factory.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>

#include <graphlab/shared_data/glshared.hpp>



#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/engine/callback/direct_callback.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>
#include <graphlab/engine/terminator/task_count_terminator.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  
  /**
   * This class defines a basic shared_memory engine
   */
  template<typename Graph, typename UpdateFunctor>
  class shared_memory_engine : 
    public iengine<Graph, UpdateFunctor> {
    
  public:

    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef typename iengine_base::graph_type graph_type;
    typedef typename iengine_base::update_functor_type update_functor_type;

    typedef typename iengine_base::vertex_id_type vertex_id_type;
    typedef typename iengine_base::ischeduler_type ischeduler_type;
   
   
    typedef typename iengine_base::termination_function_type 
    termination_function_type;
    typedef typename iengine_base::iscope_type iscope_type;
    typedef typename iengine_base::sync_function_type sync_function_type;
    typedef typename iengine_base::merge_function_type merge_function_type;

    
    typedef task_count_terminator terminator_type;
    typedef direct_callback<shared_memory_engine> callback_type;
    
    typedef general_scope_factory<graph_type> scope_factory_type;
    
  private:

    //! The Graph that the engine is operating on
    graph_type& graph;
   
    //! The local engine options
    graphlab_options opts; 

    /** 
     * The number vertices in the graph.  This number is set when
     * internal data structures are initialized and is used to track
     * whether the graph has changed externally since the engine was
     * last executed.
     */
    size_t previous_nverts;
    
    //! The active scope factory
    scope_factory_type* scope_factory_ptr;
    
    //! The active scheduler
    ischeduler_type* scheduler_ptr;
    
    //! the active terminator
    iterminator* terminator_ptr;    
   

    //! The reason for last termination
    execution_status::status_enum exec_status;

    /**
     * Local state for each thread
     */
    struct thread_state_type {
      size_t update_count;
      callback_type callback;
      char FALSE_CACHE_SHARING_PAD[128];
      thread_state_type(shared_memory_engine* ptr = NULL) : 
        update_count(0), callback(ptr) { }
    }; //end of thread_state
    std::vector<thread_state_type> tls_array; 
    

  public:
   
    //! Create an engine for the given graph
    shared_memory_engine(graph_type& graph);
    
    //! Clear internal members
    ~shared_memory_engine() { clear(); }
    
    //! Start the engine
    void start() { }
    
    //! Stop the engine
    void stop() { }
    
    //! \brief Describe the reason for termination.
    execution_status::status_enum last_exec_status() const { return exec_status; }

    //! \brief Get the number of updates executed by the engine.
    size_t last_update_count() const { return 0; }

        
    
    //! \brief Adds an update task with a particular priority.
    void schedule(vertex_id_type vid,
                  const update_functor_type& update_functor) { }

 
    //! \brief Apply update function to all the vertices in the graph
    void schedule_all(const update_functor_type& update_functor) { }


    //! \brief associate a termination function with this engine.
    void add_termination_condition(termination_function_type term) { }

    //!  remove all associated termination functions
    void clear_termination_conditions() { };
    
    //! \brief The timeout is the total
    void set_timeout(size_t timeout_secs) { }

    //! \brief set a limit on the number of tasks that may be executed.
    void set_task_budget(size_t max_tasks) { }

    //! \brief Update the engine options.  
    void set_options(const graphlab_options& newopts);

    //! \brief Get the current engine options for this engine
    const graphlab_options& get_options() { return opts; } 


    //! \brief Registers a sync with the engine.
    void set_sync(iglshared& shared,
                  sync_function_type sync,
                  iglshared::apply_function_type apply,
                  const any& zero,
                  size_t sync_interval = 0,
                  merge_function_type merge = NULL,
                  vertex_id_type rangelow = 0,
                  vertex_id_type rangehigh = -1) { };


    //! Performs a sync immediately.
    void sync_now(iglshared& shared) { };

    //! reset the engine
    void clear();

  private:
    void initialize_members();


    void run_once(size_t cpuid);
    void evaluate_sync_queue(size_t cpuid) { };
    bool evaluate_termination_conditions(size_t cpuid) { return false; };
        
  }; // end of shared_memory engine




  /////////////////////////////////////////////////////////////////////////
  /// Implementation

  template<typename Graph, typename UpdateFunctor> 
  shared_memory_engine<Graph, UpdateFunctor>::
  shared_memory_engine(graph_type& graph) : 
    graph(graph), 
    previous_nverts(graph.num_vertices()),
    scope_factory_ptr(NULL),
    scheduler_ptr(NULL),
    terminator_ptr(NULL), 
    exec_status(execution_status::EXEC_UNSET) {  
  } // end of constructor




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  set_options(const graphlab_options& new_opts) {
    clear();
    opts = new_opts;
  } // end of set_options



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  clear() {
    if(scope_factory_ptr != NULL) {
      delete scope_factory_ptr;
      scope_factory_ptr = NULL;
    }
    if(scheduler_ptr != NULL) {
      delete scheduler_ptr; 
      scheduler_ptr = NULL;
    }
    if(terminator_ptr != NULL) {
      delete terminator_ptr; 
      terminator_ptr = NULL;
    }
    tls_array.clear();
    previous_nverts = 0;
    exec_status = execution_status::EXEC_UNSET;
  } // end of clear_members




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  initialize_members() {
    // If everything is already properly initialized then we don't
    // need to do anything.
    if(previous_nverts == graph.nverts &&
       scope_factory_ptr != NULL &&
       terminator_ptr != NULL &&
       scheduler_ptr != NULL) {
      return;
    } else {
      // Force a reset and start over
      clear();
      // construct the scope factory
      ASSERT_TRUE(scope_factory_ptr == NULL);
      const consistency_model::model_enum scope_range =
        consistency_model::from_string(opts.scope_type);
      scope_factory_ptr = 
        new scope_factory_type(graph, 
                               opts.get_ncpus(),
                               scope_range);
      ASSERT_TRUE(terminator_ptr != NULL);

      // construct the terminator
      ASSERT_TRUE(terminator_ptr == NULL);                               
      terminator_ptr = new terminator_type();
      ASSERT_TRUE(terminator_ptr != NULL);

      // construct the scheduler
      ASSERT_TRUE(scheduler_ptr == NULL);
      scheduler_ptr = scheduler_factory<shared_memory_engine>::
        get_scheduler(opts.scheduler_type,
                      opts.scheduler_args,
                      graph,
                      *terminator_ptr,
                      opts.ncpus);
      ASSERT_TRUE(scheduler_ptr != NULL);

      // reset the execution status
      exec_status = execution_status::EXEC_UNSET;
      // reset the thread local state 
      const thread_state_type starting_state(this);
      tls_array.resize(opts.get_ncpus(), starting_state );
      
    }
      
  } // end of initialize_members




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  run_once(size_t cpuid) {
    // Evaluate pending sync operations
    evaluate_sync_queue(cpuid);
    // Evaluate the available termination conditions and if the
    // program is finished simply return
    const bool is_finished = evaluate_termination_conditions(cpuid);
    if(is_finished) return;

    // Get the next task from the scheduler
    vertex_id_type vid(-1);
    update_functor_type ufun;
    sched_status::status_enum stat = 
      scheduler_ptr->get_next(cpuid, vid, ufun);

    // If we failed to get a task enter the retry /termination loop
    while(stat == sched_status::EMPTY) {
      // Enter the critical section
      terminator_ptr->begin_critical_section(cpuid);
      // Try again in the critical section
      stat = scheduler_ptr->get_next(cpuid, vid, ufun);
      // If we fail again
      if (stat == sched_status::EMPTY) {
        // test if the scheduler is finished
        const bool scheduler_empty = 
          scheduler_ptr->get_terminator().end_critical_section(cpuid);
        if(scheduler_empty) {
          exec_status = execution_status::EXEC_TASK_DEPLETION;
          return;
        } else {
          // \todo use sched yield option
          // if(use_sched_yield) sched_yield();
        }
        // cancel the critical section
        terminator_ptr->cancel_critical_section(cpuid);
      }
    } // end of while loop
    ASSERT_EQ(stat, sched_status::NEW_TASK);
    ASSERT_LT(vid, graph.num_vertices());
    // Get the scope
    iscope_type* scope_ptr = 
      scope_factory_ptr->get_scope(cpuid, vid, ufun.consistency());
    ASSERT_TRUE(scope_ptr != NULL);
    // Apply the update functor
    ufun(*scope_ptr, tls_array[cpuid]);
    // Finish any pending transactions in the scope
    scope_ptr->commit();
    // Release the scope (and all corresponding locks)
    scope_factory_ptr->release_scope(scope_ptr);    
    // Mark scope as completed in the scheduler
    scheduler_ptr->completed(cpuid, vid, ufun);
    // Record an increase in the update counts
    tls_array[cpuid].update_count++;
    
   


    

    
      
  } // end of run_once

  




}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

