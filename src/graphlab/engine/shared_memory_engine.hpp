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
#include <graphlab/util/tracepoint.hpp>

#include <graphlab/logger/logger.hpp>

#include <graphlab/util/generics/any_vector.hpp>
#include <graphlab/options/graphlab_options.hpp>

#include <graphlab/graph/graph.hpp>

#include <graphlab/context/icontext.hpp>
#include <graphlab/context/iglobal_context.hpp>
#include <graphlab/context/context_manager.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>

#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/execution_status.hpp>

#include <graphlab/scheduler/terminator/iterminator.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  
  /**
   * This class defines a basic shared_memory engine
   */
  template<typename Graph, typename UpdateFunctor>
  class shared_memory_engine : public iengine<Graph, UpdateFunctor> {
    
  public:

    // Include parent types
    typedef iengine<Graph, UpdateFunctor> iengine_base;
    typedef typename iengine_base::graph_type graph_type;
    typedef typename iengine_base::update_functor_type update_functor_type;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename graph_type::edge_id_type   edge_id_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::edge_type edge_type;

    typedef ischeduler<shared_memory_engine>      ischeduler_type;
    
    typedef typename iengine_base::icontext_type  icontext_type;
    typedef context<shared_memory_engine>         context_type;
    typedef context_manager<shared_memory_engine> context_manager_type;
   
   
    typedef typename iengine_base::termination_function_type 
    termination_function_type;

    
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
    
    //! The context factory
    context_manager_type* context_manager_ptr;
    
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
      char FALSE_CACHE_SHARING_PAD[64];
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


    // Global Values ----------------------------------------------------------
    struct global_record { 
      std::vector<spinlock> locks;
      graphlab::any_vector values;
      bool is_const;
    }; // end of global_record
    typedef std::map<std::string, global_record> global_map_type;
    global_map_type global_records;

    // Global Aggregates ------------------------------------------------------
    struct isync {
      size_t interval;
      vertex_id_type begin_vid, end_vid;
      bool use_barrier;
      virtual ~isync() { }
      virtual void run_aggregator(const std::string key,
                                  const graphlab::barrier* barrier_ptr,
                                  const std::vector<mutex>* sync_vlocks_ptr,
                                  context_type context, 
                                  size_t ncpus, size_t cpuid) = 0;
    }; // end of isync
    
    template<typename Aggregator >
    struct sync : public isync {
      typedef Aggregator       aggregator_type;
      using isync::begin_vid;
      using isync::end_vid;
      using isync::use_barrier;
      const aggregator_type zero;
      aggregator_type shared_aggregator;
      mutex lock;
      sync(const aggregator_type& zero) : zero(zero), 
                                          shared_aggregator(zero) { }
      void run_aggregator(const std::string key,
                          const graphlab::barrier* barrier_ptr,
                          const std::vector<mutex>* sync_vlocks_ptr,
                          context_type context, size_t ncpus, size_t cpuid);
    }; // end of sync


    std::vector<mutex> sync_vlocks;
    mutex sync_master_lock;
    typedef std::map<std::string, isync*> sync_map_type;
    sync_map_type sync_map;
    typedef mutable_queue<std::string, long> sync_queue_type;   
    sync_queue_type sync_queue;
    thread_pool sync_threads;


    

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
        
    
    /**
     * Schedule the exeuction of an update functor on a particular
     * vertex.  
     */
    void schedule(vertex_id_type vid,
                  const update_functor_type& update_functor);

    /**
     * Schedule the execution of an update functor on a collection of
     * vertices
     */
    void schedule(const std::vector<vertex_id_type>& vid,
                  const update_functor_type& update_functor);


 
    /**
     * Schedule the execution of an update functor on all the vertices
     * in the graph.
     */
    void schedule_all(const update_functor_type& update_functor);

    /**
     * \brief associate a termination function with this engine.
     */
    void add_termination_condition(termination_function_type term);

    //!  remove all associated termination functions
    void clear_termination_conditions();
    
    //! \brief The timeout is the total
    void set_timeout(size_t timeout_in_seconds = 0);

    //! \brief set a limit on the number of tasks that may be executed.
    void set_task_budget(size_t max_tasks = 0);

    /**
     * \brief Update the engine options.  
     *
     * Setting the engine options will cause all existing state,
     * including scheduled update functors, to be cleared.
     */
    void set_options(const graphlab_options& newopts);

    //! \brief Get the current engine options for this engine
    const graphlab_options& get_options() { return opts; } 

    /**
     * Define a global mutable variable (or vector of variables).
     *
     * \param key the name of the variable (vector)
     * \param value the initial value for the variable (vector)
     * \param size the initial size of the global vector (default = 1)
     * 
     */
    template< typename T >
    void add_global(const std::string& key, const T& value, 
                    size_t size = 1);

    /**
     * Define a global constant.
     */
    template< typename T >
    void add_global_const(const std::string& key, const T& value, 
                          size_t size = 1);


    //! Change the value of a global entry
    template< typename T >
    void set_global(const std::string& key, const T& value, 
                    size_t index = 0);

    //! Get a copy of the value of a global entry
    template< typename T >
    T get_global(const std::string& key, size_t index = 0) const;


    //! \brief Registers an aggregator with the engine
    template<typename Aggregate>
    void add_aggregator(const std::string& key,            
                        const Aggregate& zero,                 
                        size_t interval,
                        bool use_barrier = false,
                        vertex_id_type begin_vid = 0,
                        vertex_id_type end_vid = 
                        std::numeric_limits<vertex_id_type>::max());
    

    //! Performs a sync immediately.
    void aggregate_now(const std::string& key);

    //! reset the engine
    void clear();

  private:
    friend class context<shared_memory_engine>;
    //! Get the global data and lock
    void get_global(const std::string& key,      
                    graphlab::any_vector*& ret_values_ptr,
                    bool& ret_is_const); 

    //! Get the global data and lock
    void acquire_global_lock(const std::string& key,      
                             size_t index = 0);
    //! Release the global data lock
    void release_global_lock(const std::string& key,      
                             size_t index = 0);


    /**
     * Initialize all the engine members.  This is called before
     * running the engine or populating the schedule.  Repeated calls
     * will only initialize once (unless the graph or options change).
     */
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
    void launch_sync_prelocked(const std::string& str, 
                               isync* sync);
    void schedule_sync_prelocked(const std::string& key, size_t sync_interval);


    void evaluate_termination_conditions(size_t cpuid);

        
  }; // end of shared_memory engine




  /////////////////////////////////////////////////////////////////////////
  /// Implementation

  template<typename Graph, typename UpdateFunctor> 
  shared_memory_engine<Graph, UpdateFunctor>::
  shared_memory_engine(graph_type& graph) : 
    graph(graph), 
    nverts(graph.num_vertices()),
    context_manager_ptr(NULL),
    scheduler_ptr(NULL),
    new_tasks_added(false),
    threads(opts.get_ncpus()),
    exec_status(execution_status::UNSET),
    exception_message(NULL),   
    start_time_millis(0) { 
  
    REGISTER_TRACEPOINT(eng_syncqueue, 
                        "shared_memory_engine: Evaluating Sync Queue");
    REGISTER_TRACEPOINT(eng_schednext, 
                        "shared_memory_engine: Reading Task from Scheduler");
    REGISTER_TRACEPOINT(eng_schedcrit, 
                        "shared_memory_engine: Time in Engine Termination Critical Section");
    REGISTER_TRACEPOINT(eng_basicupdate, 
                        "shared_memory_engine: Time in Basic Update user code");
    REGISTER_TRACEPOINT(eng_factorized, 
                        "shared_memory_engine: Time in Factorized Update user code");
    REGISTER_TRACEPOINT(eng_locktime, 
                        "shared_memory_engine: Time Acquiring Locks");
    REGISTER_TRACEPOINT(eng_evalfac, 
                        "shared_memory_engine: Total time in factorized update");
    REGISTER_TRACEPOINT(eng_evalbasic, 
                        "shared_memory_engine: Total time in factorized update");


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
    termination.timeout_millis = timeout_in_seconds * 1000;
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
    const size_t cpuid = random::fast_uniform<size_t>(0, threads.size() - 1);
    scheduler_ptr->schedule(cpuid, vid, update_functor);
  } // end of schedule

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule(const std::vector<vertex_id_type>& vids,
           const update_functor_type& update_functor) { 
    initialize_members();
    ASSERT_TRUE(scheduler_ptr != NULL);
    foreach(const vertex_id_type& vid, vids) {
      const size_t cpuid = random::fast_uniform<size_t>(0, threads.size() - 1);
      scheduler_ptr->schedule(cpuid, vid, update_functor);
    }
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
    if(context_manager_ptr != NULL) {
      delete context_manager_ptr;
      context_manager_ptr = NULL;
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
    ASSERT_TRUE(context_manager_ptr != NULL);
    ASSERT_TRUE(scheduler_ptr != NULL);
    ASSERT_EQ(tls_array.size(), opts.get_ncpus());
    ASSERT_EQ(sync_vlocks.size(), nverts);
        
    // -------------------- Reset Internal Counters ------------------------ //
    for(size_t i = 0; i < tls_array.size(); ++i) { // Reset thread local state
      tls_array[i].update_count = 0;
    }
    exec_status = execution_status::RUNNING;  // Reset active flag
    exception_message = NULL;
    start_time_millis = lowres_time_millis(); 
    // Initialize the scheduler
    scheduler_ptr->start();
    // Initialize the context manager
    context_manager_ptr->start();

    // -------------------------- Initialize Syncs ------------------------- //
    {
      typedef typename sync_map_type::value_type pair_type;
      sync_master_lock.lock();
      sync_queue.clear();
      foreach(const pair_type& pair, sync_map) {        
        // If their is a sync associated with the global record and
        // the sync has a non-zero interval than schedule it.
        if(pair.second != NULL && pair.second->interval > 0) 
          schedule_sync_prelocked(pair.first, pair.second->interval);
      }
      sync_master_lock.unlock();
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
  template<typename T>
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  add_global(const std::string& key, const T& value, size_t size) {
    global_record& record = global_records[key];
    // Set the initial value (this can change the type)
    typedef std::vector<T> vector_type;
    record.values = vector_type(size, value);
    record.is_const = false;
    record.locks.resize(size);
  } //end of set_global

  template<typename Graph, typename UpdateFunctor> 
  template<typename T>
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  add_global_const(const std::string& key, const T& value, size_t size) {
    global_record& record = global_records[key];
    // Set the initial value (this can change the type)
    typedef std::vector<T> vector_type;
    record.values = vector_type(size, value);
    record.is_const = true;
    record.locks.resize(size);
  } //end of set_global


  template<typename Graph, typename UpdateFunctor> 
  template<typename T>
  void 
  shared_memory_engine<Graph, UpdateFunctor>::
  set_global(const std::string& key, const T& value, size_t index) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;
      return;
    }
    global_record& record = iter->second;
    typedef std::vector<T> vector_type;    
    vector_type& values = record.values.template as<T>();
    ASSERT_EQ(values.size(), record.locks.size());
    ASSERT_LT(index, values.size());
    record.locks[index].lock();
    values[index] = value; 
    record.locks[index].unlock();
  } // end of set_global


  template<typename Graph, typename UpdateFunctor> 
  template<typename T>
  T
  shared_memory_engine<Graph, UpdateFunctor>::
  get_global(const std::string& key, size_t index) const {
    typename global_map_type::const_iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;      
    }
    const global_record& record = iter->second;
    typedef std::vector<T> vector_type;
    const vector_type& values = record.values.template as<T>();
    ASSERT_EQ(values.size(), record.locks.size());
    ASSERT_LT(index, values.size());
    record.locks[index].lock();
    T ret_value = values[index];
    record.locks[index].unlock();
    return ret_value;
  } //end of get_global


  template<typename Graph, typename UpdateFunctor> 
  template<typename Aggregator>
  void 
  shared_memory_engine<Graph, UpdateFunctor>::
  add_aggregator(const std::string& key,
                const Aggregator& zero,                 
                size_t interval,
                bool use_barrier,
                vertex_id_type begin_vid,
                vertex_id_type end_vid) {
    isync*& sync_ptr = sync_map[key];
    // Clear the old sync and remove from scheduling queue
    if(sync_ptr != NULL) { delete sync_ptr; sync_ptr = NULL; }
    sync_queue.remove(key);
    ASSERT_TRUE(sync_ptr == NULL);
    // Attach a new sync type
    typedef sync<Aggregator> sync_type;
    sync_ptr = new sync_type(zero);
    sync_ptr->interval    = interval;
    sync_ptr->use_barrier = use_barrier;
    sync_ptr->begin_vid   = begin_vid;
    sync_ptr->end_vid     = end_vid;
  }// end of add_sync


  
  
  template<typename Graph, typename UpdateFunctor> 
  void 
  shared_memory_engine<Graph, UpdateFunctor>::
  aggregate_now(const std::string& key) {
    initialize_members();    
    typename sync_map_type::iterator iter = sync_map.find(key);
    if(iter == sync_map.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in sync map!"
        << std::endl;
      return;
    }
    isync* sync = iter->second;
    ASSERT_NE(sync, NULL);
    // The current implementation will lead to a deadlock if called
    // from within an update function
    sync_master_lock.lock();
    launch_sync_prelocked(key, sync);
    sync_master_lock.unlock();
  } // end of sync_now




  /////////////////////////////////////////////////////////////////////////////
  ///////////////// Private Methods ///////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  get_global(const std::string& key,                         
             graphlab::any_vector*& ret_values_ptr,
             bool& ret_is_const) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) ret_values_ptr = NULL;
    else {
      // Get the global record
      global_record& rec = iter->second;
      ret_is_const = rec.is_const;      
      ret_values_ptr = &rec.values;
    }
  } // end of get_global


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  acquire_global_lock(const std::string& key, size_t index) {
    typename global_map_type::iterator iter = global_records.find(key);
    ASSERT_TRUE(iter != global_records.end());
    // Get the global record
    global_record& rec = iter->second;
    ASSERT_LT(index, rec.locks.size());
    rec.locks[index].lock();
  }

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  release_global_lock(const std::string& key, size_t index) {
    typename global_map_type::iterator iter = global_records.find(key);
    ASSERT_TRUE(iter != global_records.end());
    // Get the global record
    global_record& rec = iter->second;
    ASSERT_LT(index, rec.locks.size());
    rec.locks[index].unlock();
  } // end of release global lock


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  initialize_members() {
    // If everything is already properly initialized then we don't
    // need to do anything.
    if(nverts == graph.num_vertices() &&
       context_manager_ptr != NULL &&
       scheduler_ptr != NULL &&
       !tls_array.empty()) {
      return;
    } else {
      if(nverts != graph.num_vertices() && new_tasks_added) {
        logstream(LOG_WARNING) 
          << "The graph was modified after adding tasks! "
          << "All previously scheduled tasks will be removed!"
          << std::endl;
      } 
      // Force a reset and start over
      clear();
      
      // construct the scheduler
      ASSERT_TRUE(scheduler_ptr == NULL);
      scheduler_ptr = scheduler_factory<shared_memory_engine>::
        new_scheduler(opts.scheduler_type,
                      opts.scheduler_args,
                      graph,
                      opts.get_ncpus());
      ASSERT_TRUE(scheduler_ptr != NULL);

      // construct the context manager
      ASSERT_TRUE(context_manager_ptr == NULL);
      const consistency_model context_range =
        string_to_consistency_model(opts.get_scope_type());
      context_manager_ptr = 
        new context_manager_type(this,
                                 scheduler_ptr,
                                 &graph,
                                 opts.get_ncpus(),
                                 context_range);
      ASSERT_TRUE(context_manager_ptr != NULL);
      
      // Determine if the engine should use affinities
      std::string affinity = "false";
      opts.engine_args.get_option("affinity", affinity);
      const bool use_cpu_affinities = affinity == "true";
      if(use_cpu_affinities) 
        logstream(LOG_INFO) << "Using cpu affinities." << std::endl;
      // construct the thread pool
      threads.resize(opts.get_ncpus());
      threads.set_cpu_affinity(use_cpu_affinities);

      // reset the number of vertices
      nverts = graph.num_vertices();
      // reset the execution status
      exec_status = execution_status::UNSET;
      // reset the thread local state
      tls_array.resize(opts.get_ncpus());
      for(size_t i = 0; i < tls_array.size(); ++i) 
        tls_array[i].update_count = 0;
      // Initialize the sync data structures
      sync_vlocks.resize(graph.num_vertices());
      sync_threads.resize(opts.get_ncpus());
      sync_queue.clear();
    }      
  } // end of initialize_members





  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  launch_threads() {  
    const size_t ncpus = opts.get_ncpus();
    // launch the threads
    for(size_t i = 0; i < ncpus; ++i) {
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
    join_threads(sync_threads);
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
      } 
    }
  } // end of join_threads


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  thread_mainloop(size_t cpuid) {  
    while(exec_status == execution_status::RUNNING) { run_once(cpuid); }
    context_manager_ptr->flush_cache(cpuid);
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
    // flush the final cache
    for(size_t i = 0; i < opts.get_ncpus(); ++i) 
      context_manager_ptr->flush_cache(i);
  } // end of run_simulated


  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  run_once(size_t cpuid) {
     
    // std::cout << "Run once on " << cpuid << std::endl;
    // -------------------- Execute Sync Operations ------------------------ //
    // Evaluate pending sync operations
    BEGIN_TRACEPOINT(eng_syncqueue);
    evaluate_sync_queue();
    END_TRACEPOINT(eng_syncqueue);
    // --------------- Evaluate Termination Conditions --------------------- //
    // Evaluate the available termination conditions and if the
    // program is finished simply return
    evaluate_termination_conditions(cpuid);
    if(exec_status != execution_status::RUNNING);
    // -------------------- Get Next Update Functor ------------------------ //
    // Get the next task from the scheduler
    vertex_id_type vid(-1);
    update_functor_type ufun;
    
    BEGIN_TRACEPOINT(eng_schednext);
    sched_status::status_enum stat = 
      scheduler_ptr->get_next(cpuid, vid, ufun);
    END_TRACEPOINT(eng_schednext);
    // If we failed to get a task enter the retry /termination loop
    while(stat == sched_status::EMPTY) {
      // Enter the critical section
      BEGIN_TRACEPOINT(eng_schedcrit);
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
      }
      // cancel the critical section
      scheduler_ptr->terminator().cancel_critical_section(cpuid);
      END_TRACEPOINT(eng_schedcrit);
    } // end of while loop

    // ------------------- Run The Update Functor -------------------------- //    
    ASSERT_EQ(stat, sched_status::NEW_TASK);
    ASSERT_LT(vid, graph.num_vertices());
    // Grab the sync lock
    sync_vlocks[vid].lock();
    // Call the correct update functor
    if(ufun.is_factorizable()) {
      BEGIN_TRACEPOINT(eng_evalfac);
      evaluate_factorized_update_functor(vid, ufun, cpuid);
      END_TRACEPOINT(eng_evalfac);
    } else { 
      BEGIN_TRACEPOINT(eng_evalbasic);
      evaluate_update_functor(vid, ufun, cpuid);
      END_TRACEPOINT(eng_evalbasic);
    }
    // release the lock
    sync_vlocks[vid].unlock();
    // ----------------------- Post Update Code ---------------------------- //   
    // Mark context as completed in the scheduler
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
    BEGIN_TRACEPOINT(eng_locktime);
    // Get the context
    context_type& context = 
      context_manager_ptr->get_context(cpuid, vid, ufun.consistency());
    END_TRACEPOINT(eng_locktime);
    BEGIN_TRACEPOINT(eng_basicupdate);
    // Apply the update functor
    ufun(context);
    END_TRACEPOINT(eng_basicupdate);
    // Finish any pending transactions in the context
    context.commit();
    // Release the context (and all corresponding locks)
    context_manager_ptr->release_context(cpuid, context); 
  }

  // Factorized version
  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_factorized_update_functor(vertex_id_type vid, 
                                     update_functor_type& ufun,
                                     size_t cpuid) {

    //    std::cout << "Running vid " << vid << " on " << cpuid << std::endl;
    // Gather phase -----------------------------------------------------------
    {
      iglobal_context& context = context_manager_ptr->get_global_context(cpuid);
      ufun.init_gather(context);
    }
    const consistency_model gather_consistency = ufun.gather_consistency();
    if(ufun.gather_edges() == graphlab::IN_EDGES || 
       ufun.gather_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.in_edges(vid);

      foreach(const edge_type& edge, edges) {
        BEGIN_TRACEPOINT(eng_locktime);
        context_type& context = 
          context_manager_ptr->get_single_edge_context(cpuid, vid, edge, 
                                                       gather_consistency);
        END_TRACEPOINT(eng_locktime);
        BEGIN_TRACEPOINT(eng_factorized);
        ufun.gather(context, edge);
        END_TRACEPOINT(eng_factorized);
        context.commit();
        context_manager_ptr->release_single_edge_context(cpuid, context, edge);
      }
    }

    if(ufun.gather_edges() == graphlab::OUT_EDGES ||
       ufun.gather_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.out_edges(vid);
      foreach(const edge_type& edge, edges) {
        BEGIN_TRACEPOINT(eng_locktime);
        context_type& context = 
          context_manager_ptr->get_single_edge_context(cpuid, vid, edge, 
                                                       gather_consistency);
        END_TRACEPOINT(eng_locktime);
        BEGIN_TRACEPOINT(eng_factorized);
        ufun.gather(context, edge);
        END_TRACEPOINT(eng_factorized);
        context.commit();
        context_manager_ptr->release_single_edge_context(cpuid, context, edge);
      }
    }
    // Apply phase ------------------------------------------------------------
    BEGIN_TRACEPOINT(eng_locktime);
    context_type& context = 
      context_manager_ptr->get_vertex_context(cpuid, vid);
    END_TRACEPOINT(eng_locktime);
    BEGIN_TRACEPOINT(eng_factorized);
    ufun.apply(context);
    END_TRACEPOINT(eng_factorized);
    context.commit();
    context_manager_ptr->release_context(cpuid, context);

    // Scatter phase ----------------------------------------------------------
    const consistency_model scatter_consistency = ufun.gather_consistency();
    if(ufun.scatter_edges() == graphlab::IN_EDGES ||
       ufun.scatter_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.in_edges(vid);
      foreach(const edge_type& edge, edges) {
        BEGIN_TRACEPOINT(eng_locktime);
        context_type& context = 
          context_manager_ptr->get_single_edge_context(cpuid, vid, edge,
                                                       scatter_consistency);
        END_TRACEPOINT(eng_locktime);
        BEGIN_TRACEPOINT(eng_factorized);
        ufun.scatter(context, edge);
        END_TRACEPOINT(eng_factorized);
        context.commit();
        context_manager_ptr->release_single_edge_context(cpuid, context, edge);
      }
    }
    if(ufun.scatter_edges() == graphlab::OUT_EDGES ||
       ufun.scatter_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.out_edges(vid);
      foreach(const edge_type& edge, edges) {
        BEGIN_TRACEPOINT(eng_locktime);
        context_type& context = 
          context_manager_ptr->get_single_edge_context(cpuid, vid, edge,
                                                       scatter_consistency);
        END_TRACEPOINT(eng_locktime);
        BEGIN_TRACEPOINT(eng_factorized);
        ufun.scatter(context, edge);
        END_TRACEPOINT(eng_factorized);
        context.commit();
        context_manager_ptr->release_single_edge_context(cpuid, context, edge);
      }
    }
  } // end of evaluate_update_functor

  

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  launch_sync_prelocked(const std::string& key,
                        isync* sync) { 
    ASSERT_NE(sync, NULL);
    graphlab::barrier barrier(sync_threads.size());
    for(size_t i = 0; i < sync_threads.size(); ++i) {
      const boost::function<void (void)> sync_function = 
        boost::bind(&(isync::run_aggregator), sync, 
                    key, &barrier, &sync_vlocks, 
                    context_type(this, &graph, scheduler_ptr, i),
                    sync_threads.size(), i);
      sync_threads.launch(sync_function);
    }
    join_threads(sync_threads);
  } // end of launch sync prelocked



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_sync_queue() {
    // if the engine is no longer running or there is nothing in the
    // sync queue then we terminate early
    if(exec_status != execution_status::RUNNING || sync_queue.empty()) return;
    // Try to grab the lock if we fail just return
    if(!sync_master_lock.try_lock()) return;
    // ASSERT: the lock has been aquired. Test for a task at the top
    // of the queue
    const long negated_next_ucount = sync_queue.top().second;
    ASSERT_LE(negated_next_ucount, 0);
    const size_t next_ucount = size_t(-negated_next_ucount);
    const size_t ucount = last_update_count();
    // if we have more updates than the next update count for this
    // task then run it
    if(next_ucount < ucount) { // Run the actual sync
      const std::string key = sync_queue.top().first;
      sync_queue.pop();
      isync* sync = sync_map[key];
      ASSERT_NE(sync, NULL);
      launch_sync_prelocked(key, sync);
      // Reschedule the sync record
      schedule_sync_prelocked(key, sync->interval);
    }    
    sync_master_lock.unlock();    
  } // end of evaluate_sync_queue

  
  template<typename Graph, typename UpdateFunctor> 
  template<typename Aggregator>
  void
  shared_memory_engine<Graph, UpdateFunctor>::sync<Aggregator>::
  run_aggregator(const std::string key,
                 const graphlab::barrier* barrier_ptr,
                 const std::vector<mutex>* sync_vlocks_ptr,
                 context_type context, size_t ncpus, size_t cpuid) { 
    // Thread zero must initialize the the final shared accumulator
    const size_t nverts = context.num_vertices();
    const std::vector<mutex>& sync_vlocks = *sync_vlocks_ptr;
    const graphlab::barrier& barrier = *barrier_ptr;  
    // Compute partitioning of vertices over threads
    // Compute the true begin and end.
    const size_t global_begin = std::min(nverts, size_t(begin_vid));;
    const size_t global_end = std::min(nverts, size_t(end_vid));
    ASSERT_LE(global_begin, global_end);
    ASSERT_LE(global_end, nverts);
    // Compute the span of each subtask.  The span should not be less
    // than some minimal span.
    const size_t MIN_SPAN(1);    
    const vertex_id_type span = (global_end - global_begin)/ncpus + MIN_SPAN;    
    // Shadow the global begin
    const size_t true_begin_vid = std::min(cpuid*span, nverts);
    const size_t true_end_vid = std::min((cpuid+1)*span, nverts);

    // If Barriers are in place go ahead and lock all update functions
    if(use_barrier){
      for(vertex_id_type vid = true_begin_vid; vid < true_end_vid; ++vid) 
        sync_vlocks[vid].lock();
      barrier.wait(); 
    }
    
    // construct the local (to this thread) accumulator and context
    aggregator_type local_accum(zero);
    // Do map computation;
    for(vertex_id_type vid = true_begin_vid; vid < true_end_vid; ++vid) {
      if(!use_barrier) sync_vlocks[vid].lock();
      context.init(vid, VERTEX_CONSISTENCY);
      local_accum(context);     
      if(!use_barrier) sync_vlocks[vid].unlock();
    }
    context.commit();
    // Merge with master
    lock.lock(); shared_aggregator += local_accum; lock.unlock();
    barrier.wait();  // Wait until all merges are complete
    
    if(cpuid == 0) {
      // Recast the context as a global context.  This ensures that
      // the user implements finalize correctly;
      iglobal_context& global_context = context;
      shared_aggregator.finalize(global_context);
      context.commit();
      // Zero out the shared accumulator for the next run
      shared_aggregator = zero;
    }
    barrier.wait();   
    // If Barriers are in place go ahead and lock all update functions
    if(use_barrier){
      for(vertex_id_type vid = true_begin_vid; vid < true_end_vid; ++vid) 
        sync_vlocks[vid].unlock();
      barrier.wait(); 
    }    
  } // end of run sync



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_sync_prelocked(const std::string& key, size_t sync_interval) {
    const size_t ucount = last_update_count();
    const long negated_next_ucount = -long(ucount + sync_interval); 
    sync_queue.push(key, negated_next_ucount);
  }; // end of schedule_sync


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

