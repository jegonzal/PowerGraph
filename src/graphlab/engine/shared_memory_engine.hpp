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
#include <boost/type_traits.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/thread_pool.hpp>


#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/mutable_queue.hpp>

#include <graphlab/logger/logger.hpp>

#include <graphlab/options/graphlab_options.hpp>

#include <graphlab/graph/graph.hpp>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/scope/scope_manager.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>

#include <graphlab/shared_data/glshared.hpp>

#include <graphlab/sync/isync.hpp>
#include <graphlab/sync/sync.hpp>




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
    typedef typename iengine_base::iupdate_functor_type iupdate_functor_type;
    

    typedef typename iengine_base::vertex_id_type vertex_id_type;
    typedef typename iengine_base::edge_id_type   edge_id_type;
    typedef typename iengine_base::edge_list_type edge_list_type;

    typedef ischeduler<shared_memory_engine> ischeduler_type;
    typedef typename iengine_base::iscope_type iscope_type;
    typedef general_scope<Graph> general_scope_type;
   
    typedef typename iengine_base::termination_function_type 
    termination_function_type;

    typedef typename iengine_base::isync_type isync_type;

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
    thread_pool threads;   


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
      mutex lock;
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



    //! All sync related variables
    struct sync_members {
      //! A record contains a single sync record
      struct record {
        mutex lock;
        isync_type* sync_ptr;
        size_t sync_interval;
        vertex_id_type begin_vid, end_vid;
        size_t subtasks_remaining;
        //! If true the syncs form a barrier durring execution
        bool barrier;
        record(isync_type* sync_ptr) :
          sync_ptr(sync_ptr), 
          sync_interval(0),
          begin_vid(0), end_vid(0),
          subtasks_remaining(0) { ASSERT_TRUE(sync_ptr != NULL); }
        ~record() { clear(); }
        void clear(const bool need_lock = true) {
          if(need_lock) lock.lock();
          ASSERT_TRUE(sync_ptr != NULL);
          delete sync_ptr; 
          sync_ptr = NULL; 
          if(need_lock) lock.unlock();
        }
      }; // end of sync record,
      struct subtask {
        record* record_ptr;
        isync_type* sync_ptr;
        vertex_id_type begin_vid, end_vid;
        subtask(record* record_ptr, isync_type* sync_ptr) : 
          record_ptr(record_ptr), sync_ptr(sync_ptr), 
          begin_vid(0), end_vid(0) { 
          ASSERT_TRUE(record_ptr != NULL);
          ASSERT_TRUE(sync_ptr != NULL);
        }
        ~subtask() { 
          ASSERT_TRUE(sync_ptr != NULL);
          delete sync_ptr;
          sync_ptr = NULL;
          record_ptr = NULL;
        }
        
      }; // end of subtask

      //! This lock is held while a sync is proceeding
      mutex lock;
      //! barrier vertex locks
      std::vector<mutex> vlocks;
      //! map from iglshared2sync
      typedef std::map<iglshared*, record*> map_type;   
      map_type iglshared2sync;
      typedef mutable_queue<record*, long> queue_type;   
      queue_type queue;
      ~sync_members() {
        lock.lock();
        typedef typename map_type::value_type pair_type;
        foreach(pair_type& pair, iglshared2sync) {
          record*& record_ptr = pair.second;
          ASSERT_TRUE(record_ptr != NULL);
          delete record_ptr;
       
          record_ptr = NULL;
        
        }
        lock.unlock();
      }
      thread_pool threads;
    }; // end of sync members
    sync_members syncs;


    

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
    template< typename Accum, typename T>
    void set_sync(glshared<T>& shared,
                  size_t sync_interval = 100,
                  bool barrier = false,
                  const Accum& zero = Accum(),                
                  void(*apply_function)(T& lvalue, const Accum& rvalue) =
                  (sync_defaults::apply<T, Accum >),
                  void(*fold_function)(iscope_type& scope, Accum& acc) = 
                  (sync_defaults::fold<iscope_type, Accum>),
                  vertex_id_type begin_vid = 0,
                  vertex_id_type end_vid = 
                  std::numeric_limits<vertex_id_type>::max());


    //! Performs a sync immediately.
    void sync_now(iglshared& shared); 

    //! reset the engine
    void clear();

  private:
    void initialize_members();

    
    void launch_threads();
    void join_threads(thread_pool& threads);
    void join_all_threads();

    void thread_mainloop(size_t cpuid);
    void run_simulated();

    void run_once(size_t cpuid);

    void evaluate_update_functor(vertex_id_type vid,
                                 update_functor_type& ufun, 
                                 size_t cpuid); 
    
    void evaluate_factorized_update_functor(vertex_id_type vid,
                                              update_functor_type& ufun,
                                              size_t cpuid);
    

    
    void evaluate_sync_queue();

    /**
     * Add a sync to the scheduling queue.  This should be called
     * while holding the sync.lock.
     */
    void schedule_sync(typename sync_members::record* record_ptr);
    void background_sync(typename sync_members::record* record_ptr);
    void run_sync_subtask(typename sync_members::subtask* subtask_ptr);



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
    threads(opts.get_ncpus()),
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
    join_all_threads();
    // Delete any local datastructures if necessary
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
    ASSERT_EQ(syncs.vlocks.size(), nverts);
        
    // -------------------- Reset Internal Counters ------------------------ //
    for(size_t i = 0; i < tls_array.size(); ++i) { // Reset thread local state
      tls_array[i].update_count = 0;
    }
    exec_status = execution_status::RUNNING;  // Reset active flag
    exception_message = NULL;
    start_time_millis = lowres_time_millis(); 
    // Initialize the scheduler
    scheduler_ptr->start();

    // -------------------------- Initialize Syncs ------------------------- //
    {
      typedef typename sync_members::map_type::value_type pair_type;
      typedef typename sync_members::record record_type;
      syncs.lock.lock();
      syncs.queue.clear();
      foreach(pair_type& pair, syncs.iglshared2sync) {
        schedule_sync(pair.second);
      }
      syncs.lock.unlock();
    } // end of initialize syncs
        
    // ------------------------ Start the engine --------------------------- //
    launch_threads();
    // \todo: Decide if we really still want to support run_simulated();

    // --------------------- Finished running engine ----------------------- //
    // Join all the active threads
    join_all_threads();

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




  template<typename Graph, typename UpdateFunctor> 
  template<typename Accum, typename T>
  void 
  shared_memory_engine<Graph, UpdateFunctor>::
  set_sync(glshared<T>& shared,
           size_t sync_interval,
           bool barrier,
           const Accum& zero,
           void(*apply_function)(T& lvalue, const Accum& rvalue),
           void(*fold_function)(iscope_type& scope, Accum& acc),
           vertex_id_type begin_vid,
           vertex_id_type end_vid) {
    // Update the syncs map
    syncs.lock.lock();
    typedef typename sync_members::record record_type;
    // test if the glshared had already been registered
    if(syncs.iglshared2sync.find(&shared) != 
       syncs.iglshared2sync.end()) {
      record_type*& record_ptr = syncs.iglshared2sync[&shared];
      ASSERT_TRUE(record_ptr != NULL);
      delete record_ptr;
      record_ptr = NULL;
    }

    record_type*& record_ptr = syncs.iglshared2sync[&shared];
    ASSERT_TRUE(record_ptr == NULL);
    // Construct a local isync
    typedef fold_sync<Graph, T, Accum> fold_type;    
    isync_type* sync_ptr = new fold_type(shared,
                                         zero,
                                         fold_function,
                                         apply_function);
    record_ptr = new record_type(sync_ptr);
    record_type& rec = *record_ptr;
    // Initialize the rest of the record
    rec.lock.lock();
    rec.sync_interval = sync_interval;
    rec.begin_vid = begin_vid;
    rec.end_vid = end_vid;
    rec.barrier = barrier;
    rec.lock.unlock(); 
    syncs.lock.unlock();
  } // end of set_sync


  
  template<typename Graph, typename UpdateFunctor> 
  void 
  shared_memory_engine<Graph, UpdateFunctor>::
  sync_now(iglshared& shared) {
    typedef typename sync_members::record record_type;
    syncs.lock.lock();
    // \todo: barrier should be an engine parameter
    syncs.vlocks.resize(graph.num_vertices());
    syncs.threads.set_nthreads(opts.get_ncpus());    
    // lookup the glshared 
    ASSERT_TRUE(syncs.iglshared2sync.find(&shared) !=
                syncs.iglshared2sync.end());
    record_type* record_ptr = syncs.iglshared2sync[&shared];
    ASSERT_TRUE(record_ptr != NULL);
    background_sync(record_ptr);
    // block until sync is free
    join_threads(syncs.threads);        
  } // end of sync_now




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
      
      // Determine if the engine should use affinities
      std::string affinity = "false";
      opts.engine_args.get_string_option("affinity", affinity);
      const bool use_cpu_affinities = affinity == "true";
      // construct the thread pool
      threads.set_nthreads(opts.get_ncpus());
      threads.set_cpu_affinity(use_cpu_affinities);

      // reset the number of vertices
      nverts = graph.num_vertices();
      // reset the execution status
      exec_status = execution_status::UNSET;
      // reset the thread local state 
      const thread_state_type starting_state(this, scheduler_ptr);
      tls_array.resize(opts.get_ncpus(), starting_state );
      
      // Initialize the syncs data structure
      syncs.lock.lock();     
      syncs.vlocks.resize(graph.num_vertices());
      syncs.threads.set_nthreads(opts.get_ncpus());
      syncs.lock.unlock();
    }
      
  } // end of initialize_members





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  launch_threads() {  
    const size_t ncpus = opts.get_ncpus();
    // launch the threads
    for(size_t i = 0; i < ncpus; ++i) {
      logstream(LOG_INFO) << "Scheduling thread: " << i << std::endl;
      // Create the boost function which effectively calls:
      //    this->thread_mainloop(i);
      const boost::function<void (void)> run_function = 
        boost::bind(&shared_memory_engine::thread_mainloop, this, i);
      // Add the function to the thread pool
      threads.launch(run_function);
    }
  } // end of run_threaded
  


  

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  join_all_threads() {
    join_threads(threads);
    join_threads(syncs.threads);
  } // end of join_all_threads


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  join_threads(thread_pool& threads) {
    // Join all the threads (looping while failed)
    bool join_successful = false;   
    while(!join_successful) {
      try { 
        // Join blocks until all exection threads return.  However if
        // a thread terminates early due to exception in user code
        // (e.g., update function throws an exception) this join may
        // fail before all threads are joined.  If this happens we
        // catch the exception and then try to kill the rest of the
        // threads and proceed to join again.
        threads.join();
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
    logstream(LOG_INFO) 
      << "Starting thread " << cpuid << " started." << std::endl;
    while(exec_status == execution_status::RUNNING) {
      run_once(cpuid);
    }
    // Flush the thread local cache associated with vertex data
    scope_manager_ptr->flush_cache(cpuid);
   
    logstream(LOG_INFO) 
      << "Thread " << cpuid << " finished." << std::endl;
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
    // std::cout << "Run once on " << cpuid << std::endl;
    // -------------------- Execute Sync Operations ------------------------ //
    // Evaluate pending sync operations
    evaluate_sync_queue();

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
          continue;
        }
        // cancel the critical section
        scheduler_ptr->terminator().cancel_critical_section(cpuid);
      }
    } // end of while loop

    // ------------------- Run The Update Functor -------------------------- //
    
    ASSERT_EQ(stat, sched_status::NEW_TASK);
    ASSERT_LT(vid, graph.num_vertices());
    // Grab the sync lock
    syncs.vlocks[vid].lock();
    // Call the correct update functor
    if(ufun.is_factorizable()) {
      evaluate_factorized_update_functor(vid, ufun, cpuid);
    } else { 
      evaluate_update_functor(vid, ufun, cpuid);
    }
    // release the lock
    syncs.vlocks[vid].unlock();

    // ----------------------- Post Update Code ---------------------------- //   
    // Mark scope as completed in the scheduler
    scheduler_ptr->completed(cpuid, vid, ufun);
    // Record an increase in the update counts
    ASSERT_LT(cpuid, tls_array.size());
    tls_array[cpuid].update_count++;     
      
  } // end of run_once

  // unfactorized version
  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_update_functor(vertex_id_type vid,
                          update_functor_type& ufun, 
                          size_t cpuid) {              
    // Get the scope
    general_scope_type& scope = 
      scope_manager_ptr->get_scope(cpuid, vid, ufun.consistency());
    // get the callback
    callback_type& callback = tls_array[cpuid].callback;
    // Apply the update functor
    ufun(scope, callback);
    // Finish any pending transactions in the scope
    scope.commit();
    // Release the scope (and all corresponding locks)
    scope_manager_ptr->release_scope(cpuid, scope); 
  }


  // Factorized version
  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_factorized_update_functor(vertex_id_type vid, 
                                     update_functor_type& ufun,
                                     size_t cpuid) {
    //    std::cout << "Running vid " << vid << " on " << cpuid << std::endl;
    // get the callback
    callback_type& callback = tls_array[cpuid].callback;
    // Gather phase -----------------------------------------------------------
    if(ufun.gather_edges() == iupdate_functor_type::IN_EDGES ||
       ufun.gather_edges() == iupdate_functor_type::ALL_EDGES) {
      const edge_list_type edges = graph.in_edge_ids(vid);
      foreach(const edge_id_type eid, edges) {
        general_scope_type& scope = 
          scope_manager_ptr->get_single_edge_scope(cpuid, vid, eid, 
                                                   ufun.writable_gather());
        ufun.gather(scope, callback, eid);
        scope.commit();
        scope_manager_ptr->release_scope(cpuid, scope);
      }
    }
    if(ufun.gather_edges() == iupdate_functor_type::OUT_EDGES ||
       ufun.gather_edges() == iupdate_functor_type::ALL_EDGES) {
      const edge_list_type edges = graph.out_edge_ids(vid);
      foreach(const edge_id_type eid, edges) {
        general_scope_type& scope = 
          scope_manager_ptr->get_single_edge_scope(cpuid, vid, eid, 
                                                   ufun.writable_gather());
        ufun.gather(scope, callback, eid);
        scope.commit();
        scope_manager_ptr->release_scope(cpuid, scope);
      }
    }

    // Apply phase ------------------------------------------------------------
    general_scope_type& scope = 
      scope_manager_ptr->get_vertex_scope(cpuid, vid);
    ufun.apply(scope, callback);
    scope.commit();
    scope_manager_ptr->release_scope(cpuid, scope);

    // Scatter phase ----------------------------------------------------------
    if(ufun.scatter_edges() == iupdate_functor_type::IN_EDGES ||
       ufun.scatter_edges() == iupdate_functor_type::ALL_EDGES) {
      const edge_list_type edges = graph.in_edge_ids(vid);
      foreach(const edge_id_type eid, edges) {
        general_scope_type& scope = 
          scope_manager_ptr->get_single_edge_scope(cpuid, vid, eid,
                                                   ufun.writable_scatter());
        ufun.scatter(scope, callback, eid);
        scope.commit();
        scope_manager_ptr->release_scope(cpuid, scope);
      }
    }
    if(ufun.scatter_edges() == iupdate_functor_type::OUT_EDGES ||
       ufun.scatter_edges() == iupdate_functor_type::ALL_EDGES) {
      const edge_list_type edges = graph.out_edge_ids(vid);
      foreach(const edge_id_type eid, edges) {
        general_scope_type& scope = 
          scope_manager_ptr->get_single_edge_scope(cpuid, vid, eid,
                                                   ufun.writable_scatter());
        ufun.scatter(scope, callback, eid);
        scope.commit();
        scope_manager_ptr->release_scope(cpuid, scope);
      }
    }
  } // end of evaluate_update_functor

  




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_sync_queue() {
    typedef typename sync_members::record record_type;
    // if the engine is no longer running or there is nothing in the
    // sync queue then we terminate early
    if(exec_status != execution_status::RUNNING ||
       syncs.queue.empty()) return;
    // Try to grab the lock if we fail just return
    if(!syncs.lock.try_lock()) return;
    // ASSERT: the lock has been aquired
    // Test for a task at the top of the queue
    const long negated_next_ucount = syncs.queue.top().second;
    ASSERT_LE(negated_next_ucount, 0);
    const size_t next_ucount = size_t(-negated_next_ucount);
    const size_t last_ucount = last_update_count();
    // if we have more updates than the next update count for this
    // task then run it
    if(next_ucount < last_ucount) {
      // Get the key and then remove the task from the queue
      record_type* record_ptr = syncs.queue.pop().first;
      ASSERT_TRUE(record_ptr != NULL);
      background_sync(record_ptr);
    } else { 
      // Otherwsie the next thing to do is not ready yet so just
      // release the lock and let the thread continue.
      syncs.lock.unlock();
    }
  } // end of evaluate_sync_queue



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_sync(typename sync_members::record* record_ptr) {
    ASSERT_TRUE(record_ptr != NULL);
    const size_t ucount = last_update_count();
    const long negated_next_ucount = 
      -long(ucount + record_ptr->sync_interval); 
    syncs.queue.push(record_ptr, negated_next_ucount);
  }; // end of schedule_sync





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  background_sync(typename sync_members::record* record_ptr) {
    const size_t MIN_SPAN(1);
    typedef typename sync_members::record record_type;
    ASSERT_TRUE(record_ptr != NULL);
    // Get the record
    record_type& rec = *record_ptr;
    // Start the sync by grabbing the lock
    rec.lock.lock();
    // Zero the current master accumulator
    rec.sync_ptr->clear();
    // Compute the true begin and end.  Note that if we find that
    // the user set begin == end then we will assume the user just
    // wants to run on all vertices in the graph.
    vertex_id_type true_begin = 
      std::min(graph.num_vertices(), size_t(rec.begin_vid));;
    vertex_id_type true_end = 
      std::min(graph.num_vertices(), size_t(rec.end_vid));
    // if(true_begin == true_end) {
    //   true_begin = 0; true_end = graph.num_vertices();
    // } 
    ASSERT_LT(true_begin, true_end);
    ASSERT_LE(true_end, graph.num_vertices());
    // Compute the span of each subtask.  The span should not be
    // less than some minimal span.
    const vertex_id_type span = 
      std::max((true_end - true_begin)/opts.get_ncpus(), MIN_SPAN);
    ASSERT_GT(span, 0);
    // Launch threads to compute the fold on each span.
    for(vertex_id_type vid = true_begin;
        vid < true_end; vid += span) {
      // Create a subtask that will be passed to the thread
      typedef typename sync_members::subtask subtask_type;
      subtask_type* subtask_ptr =
        new subtask_type(record_ptr, rec.sync_ptr->clone());
      ASSERT_TRUE(subtask_ptr != NULL);
      subtask_ptr->begin_vid = vid;
      subtask_ptr->end_vid = std::min(vid + span, true_end);
      ASSERT_LT(subtask_ptr->begin_vid, subtask_ptr->end_vid);
      ASSERT_LE(subtask_ptr->end_vid, graph.num_vertices());
      // Increment the active subtasks
      rec.subtasks_remaining++;
      // Add the function to the thread pool
      const boost::function<void (void)> run_function = 
        boost::bind(&shared_memory_engine::run_sync_subtask, 
                    this, subtask_ptr);
      syncs.threads.launch(run_function);        
    } 
    // Release the locks to allow subtasks to commit their results.
    rec.lock.unlock();
  }; // end of background_sync
  





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  run_sync_subtask(typename sync_members::subtask* subtask_ptr) {
    ASSERT_TRUE(subtask_ptr != NULL);
    ASSERT_TRUE(subtask_ptr->record_ptr != NULL);
    ASSERT_TRUE(subtask_ptr->sync_ptr != NULL);
    ASSERT_LT(subtask_ptr->begin_vid, subtask_ptr->end_vid);
    ASSERT_LE(subtask_ptr->end_vid, graph.num_vertices());
    // Get a reference to the isync
    isync_type& local_sync = *(subtask_ptr->sync_ptr);
    // Get the master record
    typedef typename sync_members::record record_type;
    record_type& record = *(subtask_ptr->record_ptr);
    // Loop over all vertices applying the sync
    for(vertex_id_type vid = subtask_ptr->begin_vid;
        vid < subtask_ptr->end_vid; 
        ++vid) {
      // Get the scope;
      general_scope<Graph> scope(&graph, vid, 
                                 consistency_model::VERTEX_CONSISTENCY);
      syncs.vlocks[vid].lock();
      // Apply the sync operation
      local_sync += scope;
      // // Release the locks if we are not barriering
      if(!record.barrier) syncs.vlocks[vid].unlock();
    }
    // add the accumulator to the master
    record.lock.lock();
    *(record.sync_ptr) += local_sync;
    record.subtasks_remaining--;
    const bool last_thread = record.subtasks_remaining == 0;
    record.lock.unlock();
    /**
     * If this is the last thread than we must:
     *
     *   0) Apply the final result
     *  
     *   1) free the syncs.vlocks on all vertices if we are running in
     *   barrier mode.
     *
     *   2) Reschedule this sync task to be run again
     * 
     *   3) Trigger any other syncs that may have been pending
     *
     */
    if(last_thread) {
      record.lock.lock();
      // Step 0: 
      record.sync_ptr->apply();
      // Step 1:
      if(record.barrier) {
        // release all the locks
        for(vertex_id_type vid = 0; vid < syncs.vlocks.size(); ++vid)
          syncs.vlocks[vid].unlock();
      }
      // Step 2:
      if(exec_status == execution_status::RUNNING) {
        schedule_sync(subtask_ptr->record_ptr);
      }
      // Release the lock on the record
      record.lock.unlock();
      // Free the master lock to allow other syncs to proceed
      syncs.lock.unlock();
      // Step 3:
      if(exec_status == execution_status::RUNNING) {
        evaluate_sync_queue();
      }
    }
    // free the subtask aggregator and the subtask itself
    delete subtask_ptr;
    subtask_ptr = NULL;
  } // end of run_sync_subtask





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_termination_conditions(size_t cpuid) {
    // return immediately if a check has already occured recently
    if(termination.last_check_time_in_millis == lowres_time_millis()) return;
    // step the termination time forward
    termination.last_check_time_in_millis = lowres_time_millis();    
    // try to grab the termination check lock;
    if(!termination.lock.try_lock()) return;
    // ------------------------- Check timeout ----------------------------- //
    if(exec_status == execution_status::RUNNING &&
       termination.timeout_millis > 0 &&
       start_time_millis + termination.timeout_millis < lowres_time_millis()) {
      exec_status = execution_status::TIMEOUT;
    }
    // ----------------------- Check task budget --------------------------- //
    if(exec_status == execution_status::RUNNING &&
       termination.task_budget > 0 &&
       last_update_count() > termination.task_budget) {
      exec_status = execution_status::TASK_BUDGET_EXCEEDED;
    }
    // ------------------ Check termination functions ---------------------- //
    for (size_t i = 0; i < termination.functions.size() && 
           exec_status == execution_status::RUNNING; ++i) {
      if (termination.functions[i]()) {
        exec_status = execution_status::TERM_FUNCTION;
      }
    }
    // step the termination time forward
    termination.last_check_time_in_millis = lowres_time_millis();
    termination.lock.unlock();
  } // end of evaluate_termination_conditions

  




}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

