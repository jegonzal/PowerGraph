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
#include <graphlab/util/event_log.hpp>

#include <graphlab/logger/logger.hpp>

#include <graphlab/util/generics/any_vector.hpp>
#include <graphlab/options/graphlab_options.hpp>

#include <graphlab/graph/graph.hpp>

#include <graphlab/context/icontext.hpp>
#include <graphlab/context/context.hpp>
#include <graphlab/context/iglobal_context.hpp>
#include <graphlab/engine/rw_lock_manager.hpp>

#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/scheduler_factory.hpp>

#include <graphlab/engine/iengine.hpp>
#include <graphlab/engine/execution_status.hpp>

#include <graphlab/aggregation/shared_memory_aggregator.hpp>

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
    typedef typename iengine_base::ischeduler_type ischeduler_type;
    typedef typename iengine_base::icontext_type  icontext_type;
    
    typedef scheduler_factory<graph_type, update_functor_type> 
    scheduler_factory_type;


    typedef shared_memory_aggregator<shared_memory_engine> aggregator_type;
    friend class shared_memory_aggregator<shared_memory_engine>;
    
    typedef typename graph_type::vertex_data_type vertex_data_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;
    typedef typename graph_type::edge_id_type   edge_id_type;
    typedef typename graph_type::edge_list_type edge_list_type;
    typedef typename graph_type::edge_type edge_type;

    typedef context<shared_memory_engine>         context_type;
    typedef rw_lock_manager<graph_type>           lock_manager_type;
   
   
    
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
    
    /**
     * The lock manager manages the locking of contexts.
     */
    lock_manager_type* lock_manager_ptr;

    /** The default consistency model to use when none is provided by
        the update function */
    consistency_model default_consistency_model; 

    //! The scheduler
    ischeduler_type* scheduler_ptr;

    //! The shared memory aggregator
    aggregator_type aggregator;

    
    /** 
     * A boolean indicating that tasks have been added to the
     * scheduler since the last run 
     */
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

    /** en estimate of the update count maintined while the engine is
        running*/
    atomic<size_t> update_count;


    //! Termination related variables
    struct termination_members {
      size_t last_check_time_in_millis;
      size_t task_budget;
      size_t timeout_millis;
      mutex lock;
      termination_members() { clear(); }
      void clear() { 
        last_check_time_in_millis = 0;
        task_budget = 0;
        timeout_millis = 0;   
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
    size_t last_update_count() const { return update_count; }


    context_type get_context(const vertex_id_type vid, 
                             const consistency_model consistency =
                             DEFAULT_CONSISTENCY) {
      context_type context(this, &graph);
      context.init(vid, consistency);
      return context;
    }
        
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
    void schedule_all(const update_functor_type& update_functor,
                      const std::string& order);

    /**
     * Schedule an update on all the neighbors of a particular vertex
     */
    void schedule_in_neighbors(const vertex_id_type& vertex, 
                               const update_functor_type& update_fun);

    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    void schedule_out_neighbors(const vertex_id_type& vertex, 
                                const update_functor_type& update_fun);
                                                  
    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    void schedule_neighbors(const vertex_id_type& vertex, 
                            const update_functor_type& update_fun);

    
    //! \brief The timeout is the total
    void set_timeout(size_t timeout_in_seconds = 0);

    //! Time since start in milliseconds
    size_t elapsed_time() const;

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
                        vertex_id_type begin_vid = 0,
                        vertex_id_type end_vid = 
                        std::numeric_limits<vertex_id_type>::max());
    

    //! Performs a sync immediately.
    void aggregate_now(const std::string& key);

    //! reset the engine
    void clear();


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


  private:
    /**
     * Initialize all the engine members.  This is called before
     * running the engine or populating the schedule.  Repeated calls
     * will only initialize once (unless the graph or options change).
     */
    void initialize_members();

    template<typename FunctorType>
    void evaluate_update(vertex_id_type vid, FunctorType& fun);

    template<typename FunctorType>
    void evaluate_classic_update(vertex_id_type vid, FunctorType& fun); 
    
    template<typename FunctorType>
    void evaluate_factorized_update(vertex_id_type vid, FunctorType& ufun);


    void launch_threads();
    void join_threads(thread_pool& threads);
    void join_all_threads();

    void thread_mainloop(size_t cpuid);

    void evaluate_termination_conditions(size_t cpuid);

    DECLARE_TRACER(eng_syncqueue);
    DECLARE_TRACER(eng_schednext);
    DECLARE_TRACER(eng_schedcrit);
    DECLARE_TRACER(eng_basicupdate);
    DECLARE_TRACER(eng_factorized);
    DECLARE_TRACER(eng_locktime);
    DECLARE_TRACER(eng_evalfac);
    DECLARE_TRACER(eng_evalbasic);
    DECLARE_TRACER(eng_lockrelease);
    DECLARE_EVENT_LOG(eventlog);
    enum { SCHEDULE_EVENT = 0, UPDATE_EVENT = 1 };
  }; // end of shared_memory engine




  /////////////////////////////////////////////////////////////////////////
  /// Implementation

  template<typename Graph, typename UpdateFunctor> 
  shared_memory_engine<Graph, UpdateFunctor>::
  shared_memory_engine(graph_type& graph) : 
    graph(graph), 
    nverts(graph.num_vertices()),
    lock_manager_ptr(NULL),
    scheduler_ptr(NULL),
    aggregator(*this),
    new_tasks_added(false),
    threads(opts.get_ncpus()),
    exec_status(execution_status::UNSET),
    exception_message(NULL),   
    start_time_millis(0),
    update_count(0) {
    INITIALIZE_TRACER(eng_syncqueue,
                      "shared_memory_engine: Evaluating Sync Queue");
    INITIALIZE_TRACER(eng_schednext,
                      "shared_memory_engine: Reading Task from Scheduler");
    INITIALIZE_TRACER(eng_schedcrit,
                      "shared_memory_engine: Time in Engine Termination Critical Section");
    INITIALIZE_TRACER(eng_basicupdate,
                      "shared_memory_engine: Time in Basic Update user code");
    INITIALIZE_TRACER(eng_factorized,
                      "shared_memory_engine: Time in Factorized Update user code");
    INITIALIZE_TRACER(eng_locktime,
                      "shared_memory_engine: Time Acquiring Locks");
    INITIALIZE_TRACER(eng_evalfac,
                      "shared_memory_engine: Total time in evaluate_factorized_update_functor");
    INITIALIZE_TRACER(eng_evalbasic,
                      "shared_memory_engine: Total time in evaluate_update_functor");
    INITIALIZE_TRACER(eng_lockrelease,
                      "shared_memory_engine: Total time releasing locks");
    INITIALIZE_EVENT_LOG(eventlog, std::cout, 100, event_log::RATE_BAR);
    ADD_EVENT_TYPE(eventlog, SCHEDULE_EVENT, "Schedule");
    ADD_EVENT_TYPE(eventlog, UPDATE_EVENT, "Updates");
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
  size_t
  shared_memory_engine<Graph, UpdateFunctor>::
  elapsed_time() const {
    return lowres_time_millis() - start_time_millis;
  } // end of elapsed time

  

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
    ACCUMULATE_EVENT(eventlog, SCHEDULE_EVENT, 1);
    if (exec_status == execution_status::RUNNING) 
      scheduler_ptr->schedule_from_execution_thread(thread::thread_id(),
                                                    vid, update_functor);    
    else scheduler_ptr->schedule(vid, update_functor);    
  } // end of schedule

  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule(const std::vector<vertex_id_type>& vids,
           const update_functor_type& update_functor) { 
    initialize_members();
    ASSERT_TRUE(scheduler_ptr != NULL);
    ACCUMULATE_EVENT(eventlog, SCHEDULE_EVENT, vids.size());
    if (exec_status == execution_status::RUNNING) {
      foreach(const vertex_id_type& vid, vids) 
        scheduler_ptr->schedule_from_execution_thread(thread::thread_id(),
                                                      vid, update_functor);      
    } else {
      foreach(const vertex_id_type& vid, vids) 
        scheduler_ptr->schedule(vid, update_functor);     
    }
  } // end of schedule


  
  //! \brief Apply update function to all the vertices in the graph
  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_all(const update_functor_type& update_functor,
               const std::string& order) { 
    initialize_members();
    scheduler_ptr->schedule_all(update_functor, order);
  } // end of schedule_all




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_in_neighbors(const vertex_id_type& vertex, 
                        const update_functor_type& update_fun) {
    const edge_list_type edges = graph.in_edges(vertex);
    ACCUMULATE_EVENT(eventlog, SCHEDULE_EVENT, edges.size());
    foreach(const edge_type& e, edges) schedule(e.source(), update_fun);
  } // end of schedule in neighbors




  
  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_out_neighbors(const vertex_id_type& vertex, 
                         const update_functor_type& update_fun) {
    const edge_list_type edges = graph.out_edges(vertex);
    ACCUMULATE_EVENT(eventlog, SCHEDULE_EVENT, edges.size());
    foreach(const edge_type& e, edges) schedule(e.target(), update_fun);
  } // end of schedule out neighbors
  



  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  schedule_neighbors(const vertex_id_type& vertex, 
                     const update_functor_type& update_fun) {
    schedule_in_neighbors(vertex, update_fun);
    schedule_out_neighbors(vertex, update_fun);
  } // end of schedule neighbors


  
  




  template<typename Graph, typename UpdateFunctor> 
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  clear() {
    join_all_threads();
    // Delete any local datastructures if necessary
    if(lock_manager_ptr != NULL) {
      delete lock_manager_ptr;
      lock_manager_ptr = NULL;
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
    ASSERT_TRUE(lock_manager_ptr != NULL);
    ASSERT_TRUE(scheduler_ptr != NULL);
    // -------------------- Reset Internal Counters ------------------------ //
    update_count = 0;
    exec_status = execution_status::RUNNING;  // Reset active flag
    exception_message = NULL;
    start_time_millis = lowres_time_millis(); 
    // Initialize the scheduler
    scheduler_ptr->start();
    // -------------------------- Initialize Syncs ------------------------- //
    aggregator.initialize_queue();
    // ------------------------ Start the engine --------------------------- //
    launch_threads();
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
  stop() { exec_status = execution_status::FORCED_ABORT; } // end of stop





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
                 vertex_id_type begin_vid,
                 vertex_id_type end_vid) {
    aggregator.add_aggregator(key, zero, interval,
                              begin_vid, end_vid);
  }// end of add_sync


  
  
  template<typename Graph, typename UpdateFunctor> 
  void 
  shared_memory_engine<Graph, UpdateFunctor>::
  aggregate_now(const std::string& key) {
    initialize_members();    
    aggregator.aggregate_now(key);
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
       lock_manager_ptr != NULL &&
       scheduler_ptr != NULL) {
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
      scheduler_ptr = scheduler_factory_type::
        new_scheduler(opts.scheduler_type,
                      opts.scheduler_args,
                      graph,
                      opts.get_ncpus());
      ASSERT_TRUE(scheduler_ptr != NULL);

      // construct the context manager
      ASSERT_TRUE(lock_manager_ptr == NULL);
      default_consistency_model = 
        string_to_consistency_model(opts.get_scope_type());
      
      lock_manager_ptr = new lock_manager_type(graph);
      ASSERT_TRUE(lock_manager_ptr != NULL);
      
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
      // Resize the threads in the aggregator
      aggregator.get_threads().resize(opts.get_ncpus());     
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
    join_threads(aggregator.get_threads());
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
    for(;exec_status == execution_status::RUNNING; ++update_count) { 
      // -------------------- Execute Sync Operations ------------------------ //
      // Evaluate pending sync operations
      BEGIN_TRACEPOINT(eng_syncqueue);
      aggregator.evaluate_queue(last_update_count());
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
      // Call the correct update functor
      ACCUMULATE_EVENT(eventlog, UPDATE_EVENT, 1);
      evaluate_update(vid, ufun);
      // ----------------------- Post Update Code ---------------------------- //   
      // Mark context as completed in the scheduler
      scheduler_ptr->completed(cpuid, vid, ufun);
    }
  } // end of thread_mainloop





  template<typename Graph, typename UpdateFunctor> 
  template<typename Functor>
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_update(vertex_id_type vid, Functor& ufun) {
    if(ufun.is_factorizable()) {
      BEGIN_TRACEPOINT(eng_evalfac);
      evaluate_factorized_update(vid, ufun);
      END_TRACEPOINT(eng_evalfac);
    } else { 
      BEGIN_TRACEPOINT(eng_evalbasic);
      evaluate_classic_update(vid, ufun);
      END_TRACEPOINT(eng_evalbasic);
    }
  } // end of evaluate update functor





  // unfactorized version
  template<typename Graph, typename UpdateFunctor> 
  template<typename Functor>
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_classic_update(vertex_id_type vid, Functor& ufun) {
    BEGIN_TRACEPOINT(eng_locktime);
    // Determin the lock consistency level
    const consistency_model ufun_consistency = ufun.consistency();
    const consistency_model consistency = (ufun_consistency == DEFAULT_CONSISTENCY)?
      default_consistency_model : ufun_consistency;
    // Get and initialize the context
    context_type context(this, &graph);
    context.init(vid, consistency);
    switch(consistency) {
    case EDGE_CONSISTENCY:
    case DEFAULT_CONSISTENCY: lock_manager_ptr->lock_edge_context(vid); break;
    case VERTEX_CONSISTENCY: lock_manager_ptr->writelock_vertex(vid); break;
    case FULL_CONSISTENCY: lock_manager_ptr->lock_full_context(vid); break;
    case NULL_CONSISTENCY: break;
    } // end of lock switch
    END_AND_BEGIN_TRACEPOINT(eng_locktime, eng_basicupdate);
    // Apply the update functor
    ufun(context);
    // Finish any pending transactions in the context
    context.commit();
    END_AND_BEGIN_TRACEPOINT(eng_basicupdate, eng_lockrelease);
    // release the locks
    switch(consistency) {
    case EDGE_CONSISTENCY:
    case FULL_CONSISTENCY:
    case DEFAULT_CONSISTENCY: lock_manager_ptr->release_context(vid); break;
    case VERTEX_CONSISTENCY: lock_manager_ptr->release_vertex(vid); break;
    case NULL_CONSISTENCY: break;
    } // end of lock switch
    END_TRACEPOINT(eng_lockrelease);
  } // end of evaluate_classic







  // Factorized version
  template<typename Graph, typename UpdateFunctor> 
  template<typename Functor>
  void
  shared_memory_engine<Graph, UpdateFunctor>::
  evaluate_factorized_update(vertex_id_type vid, Functor& ufun) {
    CREATE_ACCUMULATING_TRACEPOINT(eng_locktime); 
    CREATE_ACCUMULATING_TRACEPOINT(eng_factorized);
    CREATE_ACCUMULATING_TRACEPOINT(eng_lockrelease);
    // Create a context
    context_type context(this, &graph);
    // Gather phase -----------------------------------------------------------
    const consistency_model gather_consistency = ufun.gather_consistency();
    context.init(vid, gather_consistency);
    iglobal_context& global_context = context; 
    ufun.init_gather(global_context);

    if(ufun.gather_edges() == graphlab::IN_EDGES || 
       ufun.gather_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.in_edges(vid);
      for(size_t i = 0; i < edges.size(); ++i) {
        const vertex_id_type neighbor = edges[i].source();
        // Lock the edge
        BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime);
        if(gather_consistency != NULL_CONSISTENCY) {
          if(i > 0) { // if this is not the first lock
            const vertex_id_type old_neighbor = edges[i-1].source();
            lock_manager_ptr->swap_single_edge(vid, old_neighbor, neighbor);
          } else lock_manager_ptr->lock_single_edge(vid, neighbor);
        }
        // Execute the gather
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime, eng_factorized);
        ufun.gather(context, edges[i]);
        context.commit();
        // Release the lock if necessary
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_factorized, eng_lockrelease);
        if(gather_consistency != NULL_CONSISTENCY) {
          if(i + 1 == edges.size()) // if this is the last lock
            lock_manager_ptr->release_single_edge(vid, neighbor);
        }
        END_ACCUMULATING_TRACEPOINT(eng_lockrelease);
      }
    } // end of gather in edges
    
    if(ufun.gather_edges() == graphlab::OUT_EDGES ||
       ufun.gather_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.out_edges(vid);
      for(size_t i = 0; i < edges.size(); ++i) {
        const vertex_id_type neighbor = edges[i].target();
        // Lock the edge
        BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime);
        if(gather_consistency != NULL_CONSISTENCY) {      
          if(i > 0) { // if this is not the first lock
            const vertex_id_type old_neighbor = edges[i-1].target();
            lock_manager_ptr->swap_single_edge(vid, old_neighbor, neighbor);
          } else lock_manager_ptr->lock_single_edge(vid, neighbor);
        }
        // Execute the gather
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime, eng_factorized);
        ufun.gather(context, edges[i]);
        context.commit();
        // Release the lock if necessary
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_factorized, eng_lockrelease);
        if(gather_consistency != NULL_CONSISTENCY) {
          if(i + 1 == edges.size()) // if this is the last lock
            lock_manager_ptr->release_single_edge(vid, neighbor);
        }
        END_ACCUMULATING_TRACEPOINT(eng_lockrelease);
      }
    } // end of gather out edges

    // Apply phase ------------------------------------------------------------
    context.init(vid, VERTEX_CONSISTENCY);
    BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime);
    lock_manager_ptr->writelock_vertex(vid);
    END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime, eng_factorized);
    ufun.apply(context);
    context.commit(); 
    END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_factorized, eng_lockrelease);
    lock_manager_ptr->release_vertex(vid);
    END_ACCUMULATING_TRACEPOINT(eng_lockrelease);

    // Scatter phase ----------------------------------------------------------
    const consistency_model scatter_consistency = ufun.scatter_consistency();
    context.init(vid, scatter_consistency);
    if(ufun.scatter_edges() == graphlab::IN_EDGES || 
       ufun.scatter_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.in_edges(vid);
      for(size_t i = 0; i < edges.size(); ++i) {
        const vertex_id_type neighbor = edges[i].source();
        // Lock the edge
        BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime);
        if(scatter_consistency != NULL_CONSISTENCY) {
          if(i > 0) { // if this is not the first lock
            const vertex_id_type old_neighbor = edges[i-1].source();
            lock_manager_ptr->swap_single_edge(vid, old_neighbor, neighbor);
          } else lock_manager_ptr->lock_single_edge(vid, neighbor);
        }
        // Execute the gather
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime, eng_factorized);
        ufun.scatter(context, edges[i]);
        context.commit();
        // Release the lock if necessary
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_factorized, eng_lockrelease);
        if(scatter_consistency != NULL_CONSISTENCY) {
          if(i + 1 == edges.size()) // if this is the last lock
            lock_manager_ptr->release_single_edge(vid, neighbor);
        }
        END_ACCUMULATING_TRACEPOINT(eng_lockrelease);
      }
    } // end of scatter in edges
    
    if(ufun.scatter_edges() == graphlab::OUT_EDGES ||
       ufun.scatter_edges() == graphlab::ALL_EDGES) {
      const edge_list_type edges = graph.out_edges(vid);
      for(size_t i = 0; i < edges.size(); ++i) {
        const vertex_id_type neighbor = edges[i].target();
        // Lock the edge
        BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime);
        if(scatter_consistency != NULL_CONSISTENCY) {      
          if(i > 0) { // if this is not the first lock
            const vertex_id_type old_neighbor = edges[i-1].target();
            lock_manager_ptr->swap_single_edge(vid, old_neighbor, neighbor);
          } else lock_manager_ptr->lock_single_edge(vid, neighbor);
        }
        // Execute the scatter
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_locktime, eng_factorized);
        ufun.scatter(context, edges[i]);
        context.commit();
        // Release the lock if necessary
        END_AND_BEGIN_ACCUMULATING_TRACEPOINT(eng_factorized, eng_lockrelease);
        if(scatter_consistency != NULL_CONSISTENCY) {
          if(i + 1 == edges.size()) // if this is the last lock
            lock_manager_ptr->release_single_edge(vid, neighbor);
        }
        END_ACCUMULATING_TRACEPOINT(eng_lockrelease);
      }
    } // end of scatter out edges
    STORE_ACCUMULATING_TRACEPOINT(eng_locktime);
    STORE_ACCUMULATING_TRACEPOINT(eng_factorized);
    STORE_ACCUMULATING_TRACEPOINT(eng_lockrelease);
  } // end of evaluate_update_functor




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
    // update the update_count estimate
    if(exec_status == execution_status::RUNNING &&
       termination.task_budget > 0 && 
       last_update_count() > termination.task_budget) {
      exec_status = execution_status::TASK_BUDGET_EXCEEDED;
    }
    // step the termination time forward
    termination.last_check_time_in_millis = lowres_time_millis();
    termination.lock.unlock();
  } // end of evaluate_termination_conditions

  




}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

