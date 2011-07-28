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
#include <graphlab/parallel/thread_pool.hpp>


#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>


#include <graphlab/logger/logger.hpp>

#include <graphlab/options/graphlab_options.hpp>

#include <graphlab/graph/graph.hpp>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/scope/scope_manager.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>

#include <graphlab/shared_data/glshared.hpp>
#include <graphlab/aggregation/fold_combine_aggregator.hpp>




#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/execution_status.hpp>
#include <graphlab/engine/callback/direct_callback.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>

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

    typedef ischeduler<shared_memory_engine> ischeduler_type;
    typedef typename iengine_base::iscope_type iscope_type;
   
    typedef typename iengine_base::termination_function_type 
    termination_function_type;

    // typedef typename iengine_base::sync_function_type sync_function_type;
    // typedef typename iengine_base::merge_function_type merge_function_type;

    
    typedef direct_callback<shared_memory_engine> callback_type;
    
    typedef scope_manager<graph_type> scope_manager_type;
    
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
    size_t nverts;
    
    //! The scope factory
    scope_manager_type* scope_manager_ptr;
    
    //! The scheduler
    ischeduler_type* scheduler_ptr;
    
    //! A boolean indicating that tasks have been added to the
    //! scheduler since the last run
    bool new_tasks_added;


    /**
     * The thread pool contain the active threads.  The pool is stored
     * as a pointer so that it does not need to be constructed in the
     * event that the engine is being simulated.  This helps reduce
     * unecessary overhead.
     */
    thread_pool* threads_ptr;   


    /**
     * The execution status describes the current state of the engine.
     * If unset the engine is not currently running.  
     *
     * The main thread loops check the execution status while running
     * and continue to run if and only iff the execution status is
     * currently set to RUNNING.  
     */
    execution_status::status_enum exec_status;

    /** The cause of the last termination condition */
    const char* exception_message;
    
    /** The time in millis that the engine was started */
    size_t start_time_millis;


    /**
     * Local state for each thread
     */
    struct thread_state_type {
      size_t update_count;
      callback_type callback;
      char FALSE_CACHE_SHARING_PAD[64];
      thread_state_type(shared_memory_engine* engine_ptr = NULL,
                        ischeduler_type* ischeduler_ptr = NULL) : 
        update_count(0), callback(engine_ptr, ischeduler_ptr) { }
    }; //end of thread_state
    std::vector<thread_state_type> tls_array; 



    //! Termination related variables
    struct termination_members {
      size_t last_check_time_in_millis;
      size_t task_budget;
      size_t timeout_millis;
      std::vector<termination_function_type> functions;
      termination_members() { clear(); }
      void clear() { 
        last_check_time_in_millis = 0;
        task_budget = 0;
        timeout_millis = 0;
        functions.clear();
      }
    };
    termination_members termination;


    

  public:
   
    //! Create an engine for the given graph
    shared_memory_engine(graph_type& graph);
    
    //! Clear internal members
    ~shared_memory_engine() { clear(); }
    
    //! Start the engine
    void start();
    
    //! Stop the engine
    void stop();
    
    //! \brief Describe the reason for termination.
    execution_status::status_enum last_exec_status() const { 
      return exec_status; }

    //! \brief Get the number of updates executed by the engine.
    size_t last_update_count() const;

        
    
    //! \brief Adds an update task with a particular priority.
    void schedule(vertex_id_type vid,
                  const update_functor_type& update_functor);

 
    //! \brief Apply update function to all the vertices in the graph
    void schedule_all(const update_functor_type& update_functor);


    //! \brief associate a termination function with this engine.
    void add_termination_condition(termination_function_type term);

    //!  remove all associated termination functions
    void clear_termination_conditions();
    
    //! \brief The timeout is the total
    void set_timeout(size_t timeout_in_seconds = 0);

    //! \brief set a limit on the number of tasks that may be executed.
    void set_task_budget(size_t max_tasks = 0);

    //! \brief Update the engine options.  
    void set_options(const graphlab_options& newopts);

    //! \brief Get the current engine options for this engine
    const graphlab_options& get_options() { return opts; } 


    //! \brief Registers a sync with the engine.
    template<typename Target, typename Accum>
    void set_sync(Target& target,
                  const Accum zero,
                  size_t sync_interval,
                  typename fold_combine_aggregator<Graph,Target,Accum>::
                  fold_function_type fold_fun,
                  typename fold_combine_aggregator<Graph,Target,Accum>::
                  apply_function_type apply_fun,
                  typename fold_combine_aggregator<Graph,Target,Accum>::
                  combine_function_type combine_fun,
                  size_t rangelow = 0,
                  size_t rangehigh = 0) { 
      
    }



    //! Performs a sync immediately.
    template<typename Target>
    void sync_now(Target& t) { };

    //! reset the engine
    void clear();

  private:
    void initialize_members();

    
    void launch_threads();
    void join_threads();
    void thread_mainloop(size_t cpuid);
    void run_simulated();


    void run_once(size_t cpuid);
    void evaluate_sync_queue(size_t cpuid);
    void evaluate_termination_conditions(size_t cpuid);
        
  }; // end of shared_memory engine




  /////////////////////////////////////////////////////////////////////////
  /// Implementation

  template<typename Graph, typename UpdateFunctor> 
  shared_memory_engine<Graph, UpdateFunctor>::
  shared_memory_engine(graph_type& graph) : 
    graph(graph), 
    nverts(graph.num_vertices()),
    scope_manager_ptr(NULL),
    scheduler_ptr(NULL),
    new_tasks_added(false),
    threads_ptr(NULL),
    exec_status(execution_status::UNSET),
    exception_message(NULL),   
    start_time_millis(0) {
    
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
  set_timeout(size_t timeout_in_seconds) {
    termination.timeout_millis = 
      timeout_in_seconds * 1000;
  } // end of set_timeout

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  set_task_budget(size_t max_tasks) {
     termination.task_budget = max_tasks; 
  } // end of set_timeout



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule(vertex_id_type vid,
           const update_functor_type& update_functor) { 
    initialize_members();
    ASSERT_TRUE(scheduler_ptr != NULL);
    scheduler_ptr->schedule(vid, update_functor);
  } // end of schedule

 
    //! \brief Apply update function to all the vertices in the graph
  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_all(const update_functor_type& update_functor) { 
    initialize_members();
    ASSERT_TRUE(scheduler_ptr != NULL);
    scheduler_ptr->schedule_all(update_functor);
  } // end of schedule_all


  template<typename Graph, typename UpdateFunctor> 
  size_t
  shared_memory_engine<Graph, UpdateFunctor>::
  last_update_count() const {
    size_t count = 0;
    for(size_t i = 0; i < tls_array.size(); ++i) 
      count += tls_array[i].update_count;
    return count;
  } // end of last update count



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  clear() {
    // Join any pending threads before proceeding to destroy local
    // data structures
    if(threads_ptr != NULL) {
      join_threads();
      delete threads_ptr;
      threads_ptr = NULL;
    }
    if(scope_manager_ptr != NULL) {
      delete scope_manager_ptr;
      scope_manager_ptr = NULL;
    }
    if(scheduler_ptr != NULL) {
      if(new_tasks_added) {
        logstream(LOG_WARNING) 
          << "Destroying scheduler without excuting new tasks!"
          << std::endl;
        new_tasks_added = false;
      }
      delete scheduler_ptr; 
      scheduler_ptr = NULL;
    }
    tls_array.clear();
    nverts = 0;
    exec_status = execution_status::UNSET;
  } // end of clear_members









  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  start() {      
    // ---------------------- Pre-execution Checks ------------------------- //
    graph.finalize();      // Do any last checks on the graph
    initialize_members();  // Initialize engine members
    
    // Check internal data-structures
    ASSERT_EQ(graph.num_vertices(), nverts);
    ASSERT_TRUE(scope_manager_ptr != NULL);
    ASSERT_TRUE(scheduler_ptr != NULL);
    ASSERT_EQ(tls_array.size(), opts.get_ncpus());
        
    // -------------------- Reset Internal Counters ------------------------ //
    for(size_t i = 0; i < tls_array.size(); ++i) { // Reset thread local state
      tls_array[i].update_count = 0;
    }
    exec_status = execution_status::RUNNING;  // Reset active flag
    exception_message = NULL;
    start_time_millis = lowres_time_millis(); 
    // Initialize the scheduler
    scheduler_ptr->start();
      
    // ------------------------ Start the engine --------------------------- //
    if(threads_ptr != NULL) launch_threads();
    else run_simulated();

    // --------------------- Finished running engine ----------------------- //
    // \todo: new_tasks_added should only be cleared when no tasks
    // remain in the scheduler
    new_tasks_added = false;  // The scheduler is "finished"
    
  } // end of start







  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  stop() {      
    // Setting the execution status to forced abort will cause any
    // thread loops to terminate after the current computaiton
    // completes
    exec_status = execution_status::FORCED_ABORT;
  } // end of stop







  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  clear_termination_conditions() {      
    termination.clear();
  } // end of clear_termination_conditions


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  add_termination_condition(termination_function_type fun) {
    ASSERT_TRUE(fun != NULL);
    termination.functions.push_back(fun);
  } // end of clear_termination_conditions




  /////////////////////////////////////////////////////////////////////////////
  ///////////////// Private Methods ///////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  initialize_members() {
    // If everything is already properly initialized then we don't
    // need to do anything.
    if(nverts == graph.num_vertices() &&
       scope_manager_ptr != NULL &&
       scheduler_ptr != NULL) {
      return;
    } else {
      if(nverts != graph.num_vertices() && new_tasks_added) {
        logstream(LOG_WARNING) 
          << "The was modified after adding tasks! All tasks will be removed."
          << std::endl;
      } 
      // Force a reset and start over
      clear();
      // construct the scope factory
      ASSERT_TRUE(scope_manager_ptr == NULL);
      const consistency_model::model_enum scope_range =
        consistency_model::from_string(opts.get_scope_type());
      scope_manager_ptr = 
        new scope_manager_type(graph, 
                               opts.get_ncpus(),
                               scope_range);
      ASSERT_TRUE(scope_manager_ptr != NULL);


      // construct the scheduler
      ASSERT_TRUE(scheduler_ptr == NULL);
      scheduler_ptr = scheduler_factory<shared_memory_engine>::
        new_scheduler(opts.scheduler_type,
                      opts.scheduler_args,
                      graph,
                      opts.get_ncpus());
      ASSERT_TRUE(scheduler_ptr != NULL);
      
      // construct the thread pool (IF NOT SIMUALTED)
      ASSERT_TRUE(threads_ptr == NULL);
      // Determine if the engine is threaded or simualted (default is simulated)
      std::string engine_type = "threaded"; // "simulated"
      opts.engine_args.get_string_option("type", engine_type);
      const bool is_simulated = engine_type == "simulated" ||
        opts.get_engine_type() == "async_sim";
      if(!is_simulated) {
        // Determine if the engine should use affinities
        std::string affinity = "false";
        opts.engine_args.get_string_option("affinity", affinity);
        const bool use_cpu_affinities = affinity == "true";
        // construct the thread pool
        threads_ptr = new thread_pool(opts.get_ncpus(), use_cpu_affinities);
        ASSERT_TRUE(threads_ptr != NULL);
      }     

      // reset the number of vertices
      nverts = graph.num_vertices();
      // reset the execution status
      exec_status = execution_status::UNSET;
      // reset the thread local state 
      const thread_state_type starting_state(this, scheduler_ptr);
      tls_array.resize(opts.get_ncpus(), starting_state );
    }
      
  } // end of initialize_members





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  launch_threads() {  
    ASSERT_TRUE(threads_ptr != NULL);
    const size_t ncpus = opts.get_ncpus();
    // launch the threads
    for(size_t i = 0; i < ncpus; ++i) {
      // Create the boost function which effectively calls:
      //    this->thread_mainloop(i);
      const boost::function<void (void)> run_function = 
        boost::bind(&shared_memory_engine::thread_mainloop, this, i);
      // Add the function to the thread pool
      threads_ptr->launch(run_function);
    }
    // Join the threads
    join_threads();
  } // end of run_threaded
  


  


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  join_threads() {
    ASSERT_TRUE(threads_ptr != NULL);
    // Join all the threads (looping while failed)
    bool join_successful = false;   
    while(!join_successful) {
      try { 
        // Join blocks until all exection threads return.  However if
        // a thread terminates early due to exception in user code
        // (e.g., update function) this join may fail before all
        // threads are joined.  If this happens we catch the exception
        // and then try to kill the rest of the threads and proceed to
        // join again.
        threads_ptr->join();
        join_successful = true; 
      } catch (const char* error_str) {
        logstream(LOG_ERROR) 
          << "Exception in execution thread caught: " << error_str 
          << std::endl;     
        // killall the running threads
        exception_message = error_str;
        exec_status = execution_status::EXCEPTION;
        join_successful = false;
      } 
    }
  } // end of join_threads




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  thread_mainloop(size_t cpuid) {  
    while(exec_status == execution_status::RUNNING) {
      run_once(cpuid);
    }
  } // end of thread_mainloop





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  run_simulated() {  
    ASSERT_TRUE(opts.get_ncpus() > 0);
    const size_t last_cpuid = opts.get_ncpus() - 1;
    while(exec_status == execution_status::RUNNING) {
      const size_t cpuid = random::fast_uniform<size_t>(0, last_cpuid);
      run_once(cpuid);
    }
  } // end of run_simulated






  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  run_once(size_t cpuid) {
    // -------------------- Execute Sync Operations ------------------------ //
    // Evaluate pending sync operations
    evaluate_sync_queue(cpuid);

    // --------------- Evaluate Termination Conditions --------------------- //
    // Evaluate the available termination conditions and if the
    // program is finished simply return
    evaluate_termination_conditions(cpuid);
    if(exec_status != execution_status::RUNNING);

    // -------------------- Get Next Update Functor ------------------------ //
    // Get the next task from the scheduler
    vertex_id_type vid(-1);
    update_functor_type ufun;
    sched_status::status_enum stat = 
      scheduler_ptr->get_next(cpuid, vid, ufun);
    // If we failed to get a task enter the retry /termination loop
    while(stat == sched_status::EMPTY) {
      // Enter the critical section
      scheduler_ptr->terminator().begin_critical_section(cpuid);
      // Try again in the critical section
      stat = scheduler_ptr->get_next(cpuid, vid, ufun);
      // If we fail again
      if (stat == sched_status::EMPTY) {
        // test if the scheduler is finished
        const bool scheduler_empty = 
          scheduler_ptr->terminator().end_critical_section(cpuid);
        if(scheduler_empty) {
          exec_status = execution_status::TASK_DEPLETION;
          return;
        } else {
          // \todo use sched yield option
          // if(use_sched_yield) sched_yield();
        }
        // cancel the critical section
        scheduler_ptr->terminator().cancel_critical_section(cpuid);
      }
    } // end of while loop

    // ------------------------- Update Code ------------------------------- //
    ASSERT_EQ(stat, sched_status::NEW_TASK);
    ASSERT_LT(vid, graph.num_vertices());
    // Get the scope
    iscope_type& scope = 
      scope_manager_ptr->get_scope(cpuid, vid, ufun.consistency());
    // Apply the update functor
    ufun(scope, tls_array[cpuid].callback);
    // ----------------------- Post Update Code ---------------------------- //
    // Finish any pending transactions in the scope
    scope.commit();
    // Release the scope (and all corresponding locks)
    scope_manager_ptr->release_scope(scope);    
    // Mark scope as completed in the scheduler
    scheduler_ptr->completed(cpuid, vid, ufun);
    // Record an increase in the update counts
    ASSERT_LT(cpuid, tls_array.size());
    tls_array[cpuid].update_count++;     
      
  } // end of run_once

  




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_sync_queue(size_t cpuid) {


  } // end of evaluate_sync_queue






  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_termination_conditions(size_t cpuid) {
    return;
    // return immediately if a check has already occured recently
    if(termination.last_check_time_in_millis == lowres_time_millis()) return;
    termination.last_check_time_in_millis = lowres_time_millis();    
    // ------------------------- Check timeout ----------------------------- //
    if(termination.timeout_millis > 0 &&
       start_time_millis + termination.timeout_millis < lowres_time_millis()) {
      exec_status = execution_status::TIMEOUT;
      return;
    }
    // ----------------------- Check task budget --------------------------- //
    if(termination.task_budget > 0 &&
       last_update_count() > termination.task_budget) {
      exec_status = execution_status::TASK_BUDGET_EXCEEDED;
      return;
    }
    // ------------------ Check termination functions ---------------------- //
    for (size_t i = 0; i < termination.functions.size(); ++i) {
      if (termination.functions[i]()) {
        exec_status = execution_status::TERM_FUNCTION;
        return;
      }
    }
  } // end of evaluate_termination_conditions

  




}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

