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


#ifndef DISTRIBUTED2_DELTA_ENGINE_HPP
#define DISTRIBUTED2_DELTA_ENGINE_HPP
#include <functional>
#include <algorithm>
#include <ext/functional> // for select1st
#include <boost/bind.hpp>
#include <boost/unordered_map.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/options/graphlab_options.hpp>
#include <graphlab/engine/execution_status.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/async_consensus.hpp>
#include <graphlab/distributed2/graph/dgraph_context.hpp>
#include <graphlab/distributed2/distributed_shared_data.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {



/**
All processes must receive the same options at the same time.
i.e. if set_cpu_affinities is called, all processes mus call it at the same time.
This is true for all set_* functions.
*/
template <typename Graph, typename UpdateFunctor>
class distributed_delta_engine : public iengine<Graph, UpdateFunctor> {
 public:
  // Include parent types
  typedef iengine<Graph, UpdateFunctor> iengine_base;
  typedef distributed_delta_engine<Graph, UpdateFunctor> engine_type;
  typedef Graph graph_type;
  typedef typename iengine_base::update_functor_type update_functor_type;
  
  typedef typename graph_type::vertex_data_type vertex_data_type;
  typedef typename graph_type::vertex_id_type vertex_id_type;
  typedef typename graph_type::edge_id_type   edge_id_type;
  typedef typename graph_type::edge_list_type edge_list_type;
  typedef typename graph_type::edge_type edge_type;  
  
  typedef typename iengine_base::icontext_type  icontext_type;
  typedef dgraph_context<distributed_delta_engine>  context_type;
  typedef distributed_shared_data<context_type> shared_data_type;
  
  typedef typename iengine_base::termination_function_type termination_function_type;


 private:
  // the local rmi instance
  dc_dist_object<distributed_delta_engine<Graph, UpdateFunctor> > rmi;
  
  // the graph we are processing
  Graph &graph;
  
  /** Number of cpus to use */
  size_t ncpus; 

  /** Use processor affinities */
  bool use_cpu_affinity;

  /** Use schedule yielding when waiting on the scheduler*/
  bool use_sched_yield;
 
  
  /** Track the number of updates */
  std::vector<size_t> update_counts;

  atomic<size_t> numsyncs;
  
  shared_data_type shared_data;
  
  /** terminators */ 
  std::vector<termination_function_type> term_functions;

  /** The timeout time in millis */
  size_t timeout_millis;
  float start_time;
  timer ti;

  
  /// Used to identify when the engine is stopped through stop
  bool force_stop;
  
  /** The total number of tasks that should be executed */
  size_t task_budget;
  

  
  /** If dynamic scheduling is used, the number of scheduled tasks */
  atomic<size_t> num_pending_tasks; 
  
  /** The cause of the last termination condition */
  execution_status::status_enum termination_reason;

  consistency_model::model_enum default_scope_range;
 
  std::vector<std::vector<vertex_id_type> > color_block; // set of localvids in each color
  vertex_functor_set<UpdateFunctor> vfunset;
 
  graphlab_options opts;
  
  double barrier_time;
  size_t num_dist_barriers_called;
  async_consensus consensus;
  
  struct sync_task {
    std::string key;
    size_t sync_interval;
    size_t next_time;
    vertex_id_type rangelow;
    vertex_id_type rangehigh;
  };
  
  /// A list of all registered sync tasks
  std::vector<sync_task> sync_tasks;
  
  /// The list of tasks which are currently being evaluated
  std::vector<sync_task*> active_sync_tasks;

  std::vector<mutex> synclock;
  
 public:
  distributed_delta_engine(distributed_control &dc,
                                    Graph& graph,
                                    size_t ncpus = 1):
                            rmi(dc, this),
                            graph(graph),
                            ncpus( std::max(ncpus, size_t(1)) ),
                            use_cpu_affinity(false),
                            use_sched_yield(true),
                            update_counts(std::max(ncpus, size_t(1)), 0),
                            shared_data(dc, ncpus),
                            timeout_millis(0),
                            start_time(0),
                            force_stop(false),
                            task_budget(0),
                            termination_reason(execution_status::UNSET),
                            vfunset(graph.get_local_store().num_vertices()),
                            barrier_time(0.0),
                            consensus(dc, ncpus, &rmi),
                            synclock(ncpus),
                            thread_color_barrier(ncpus){ 
                              
    permuted_ordering.resize(graph.get_local_store().num_vertices());
    for (size_t i = 0;i < graph.get_local_store().num_vertices(); ++i) {
      permuted_ordering[i] = i;
    }
    std::random_shuffle(permuted_ordering.begin(), permuted_ordering.end());
    rmi.barrier();
  }
  
  ~distributed_delta_engine() {
    rmi.barrier();
  }
  
  
  //! Get the number of cpus
  size_t get_ncpus() const { return ncpus; }

  //! set sched yield
  void set_sched_yield(bool value) {
    use_sched_yield = value;
    rmi.barrier();
  }

  void set_cpu_affinities(bool value) {
    use_cpu_affinity = value;
    rmi.barrier();
  }


  /**
   * Set the default scope range.  The scope ranges are defined in
   * iscope.hpp
   */
  void set_default_scope(consistency_model::model_enum default_scope_range_) {
    default_scope_range = default_scope_range_;
    rmi.barrier();
  }
  
  

  /**
   * Return the reason why the engine last terminated
   */
  execution_status::status_enum last_exec_status() const {
    return termination_reason;
  }

  /**
   * This function computes the last update count by adding all the
   * update counts of the individual threads.  This is an underestimate
   * if the engine is currently running.
   */
  size_t thisproc_update_counts() const {
    size_t sum = 0;
    for(size_t i = 0; i < update_counts.size(); ++i)
      sum += update_counts[i];
    return sum;
  } 


  /**
   * Returns the total number of updates executed
   */
  size_t last_update_count() const {
    return total_update_count;
  } // end of last_update_count


  /**
   * Add a terminator to the engine.
   * Must be called by all machines simultaneously.
   */
  void add_termination_condition(termination_function_type term) {
    term_functions.push_back(term);
    rmi.barrier();
  }


  /**
   * Clear all terminators from the engine
   * Must be called by all machines simultaneously.
   */
   void clear_termination_conditions() {
    term_functions.clear();
    rmi.barrier();
  }



  /**
   * Set a timeout. Disabled if set to 0.
   * Must be called by all machines simultaneously.
   */
  void set_timeout(size_t timeout_seconds = 0) {
    timeout_millis = timeout_seconds * 1000;
    rmi.barrier();
  }
  
  /**
   * Sets a Task budget - max number of tasks to allow.
   * Disabled if set to 0.
   * Must be called by all machines simultaneously.
   */
  void set_task_budget(size_t max_tasks) {
    task_budget = max_tasks;
    rmi.barrier();
  }
  
  void schedule_from_context(size_t threadid, 
                             const vertex_id_type vid, 
                             const update_functor_type& fun) {
    if (vfunset.add(vid, fun)) {
      num_pending_tasks.inc();
    }
  }

  void schedule_impl(const vertex_id_type vid, 
                const update_functor_type& fun) {
    if (graph.is_owned(vid)) {
      if (vfunset.add(graph.globalvid_to_localvid(vid), fun)) {
        num_pending_tasks.inc();
      }
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(vid),
                &distributed_delta_engine<Graph, UpdateFunctor>::schedule_impl,
                vid, fun);
    }
  }


  /**
   * \brief Adds an update task with a particular priority.
   * add_task on any vertex can be called by any machine.
   * The call is asynchronous and may not be completed until
   * a full_barrier is issued.
   */
  void schedule(const vertex_id_type vid, 
                const update_functor_type& fun) {
    schedule_impl(vid, fun);
  }

  /**
   * \brief Creates a collection of tasks on all the vertices in
   * 'vertices', and all with the same update function and priority
   * This function is forwarded to the scheduler.
   */
  void schedule(const std::vector<vertex_id_type>& vertices,
                const update_functor_type& update_functor) {
    // not the most efficient way to do it...
    for (size_t i = 0;i < vertices.size(); ++i) {
      schedule(vertices[i], update_functor);
    }
  }

  void schedule_collection(const std::vector<std::pair<vertex_id_type, 
                                                     update_functor_type> >& tasks) {
    typename std::vector<std::pair<vertex_id_type, update_functor_type> >::const_iterator iter = tasks.begin();
    while(iter != tasks.end()) {
      schedule_impl(iter->first, iter->second);
      ++iter;
    }
    consensus.cancel();
  }
  
  /**
   * \brief Creates a collection of tasks on all the vertices in the graph,
   * with the same update function and priority
   * Must be called by all machines simultaneously
   */
  void schedule_all(const update_functor_type& update_functor) { 
    schedule(graph.owned_vertices(), update_functor);
    rmi.barrier();
  }

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
                  size_t size = 1) {
    shared_data.add_global(key, value, size);
  }

  /**
    * Define a global constant.
    */
  template< typename T >
  void add_global_const(const std::string& key, const T& value, 
                        size_t size = 1) {
    shared_data.add_global_const(key, value, size);
  }


  //! Change the value of a global entry
  template< typename T >
  void set_global(const std::string& key, const T& value, 
                  size_t index = 0) {
    shared_data.set_global(key, value, index);
    shared_data.synchronize_global(key, value, index);
  }

  //! Get a copy of the value of a global entry
  template< typename T >
  T get_global(const std::string& key, size_t index = 0) const {
    return shared_data.get_global<T>(key, index);
  }



  /**
    Registers a sync operation.
    Must be called by all machine simultaneously
  */
    template<typename Accum>
    void add_sync(const std::string& key,            
                  const Accum& zero,                 
                  size_t sync_interval,
                  bool use_barrier = false,
                  vertex_id_type begin_vid = 0,
                  vertex_id_type end_vid = 
                  std::numeric_limits<vertex_id_type>::max()) {
    sync_task st;
    st.key = key;
    st.sync_interval = sync_interval;
    st.next_time = 0;
    st.rangelow = begin_vid;
    st.rangehigh = end_vid;
    sync_tasks.push_back(st);
    shared_data.add_sync(key, zero);
    rmi.barrier();
  }

  
  
  /************  Actual Execution Engine ****************/
 private:

  atomic<size_t> curidx;
  barrier thread_color_barrier;
 public: 
  
  struct termination_evaluation{
    size_t pending_tasks;
    size_t executed_tasks;
    bool terminator;
    bool timeout;
    bool force_stop;
    termination_evaluation(): pending_tasks(0),
                              executed_tasks(0),
                              terminator(false),
                              timeout(false),
                              force_stop(false) { }
                              
    void save(oarchive &oarc) const {
      oarc << pending_tasks
           << executed_tasks
           << terminator
           << timeout
           << force_stop;
    }
    
    void load(iarchive &iarc) {
      iarc >> pending_tasks
           >> executed_tasks
           >> terminator
           >> timeout
           >> force_stop;
    }
  };

  /**
   * Initialize the sync tasks. Called by start()
   */
  void init_syncs() {
     active_sync_tasks.clear();
     shared_data.reset_all_syncs();
     for (size_t i = 0;i < sync_tasks.size(); ++i) {
       // everyone runs at the start even if scheduling interval is 0
       active_sync_tasks.push_back(&(sync_tasks[i]));
     }
  }


  
  /**
   * Called whenever a vertex is executed.
   * Accumulates the available syncs
   */
  void eval_syncs(vertex_id_type curvertex, 
                  context_type& context, 
                  size_t threadid) {
    synclock[threadid].lock();
     // go through all the active sync tasks
     foreach(sync_task* task, active_sync_tasks) {
       // if in range, sync!
       if (task->rangelow <= curvertex && curvertex <= task->rangehigh) {
         shared_data.accumulate(task->key, &context, threadid);
       }
     }
     synclock[threadid].unlock();
  }

  /** Called at the end of the iteration. Called by all threads after a barrier*/
  void sync_end_iteration(size_t threadid, context_type& globalcontext) {
    if (threadid == 0) {
      foreach(sync_task* task, active_sync_tasks) {
        shared_data.finalize(task->key, &globalcontext);
      }
      shared_data.wait_for_all_communication();
      shared_data.reset_all_syncs();
    }
    thread_color_barrier.wait();
  }

  /** clears the active sync tasks and figure out what syncs to run next.
      Called by one thread from each machine after sync_end_iteration */
  void compute_sync_schedule(size_t num_executed_tasks) {
    // update the next time variable
    for (size_t i = 0;i < active_sync_tasks.size(); ++i) {
      sync_tasks[i].next_time = num_executed_tasks + sync_tasks[i].sync_interval;
      //  if it is within the next iteration 
      if (sync_tasks[i].next_time < num_executed_tasks + graph.num_vertices()) sync_tasks[i].next_time = num_executed_tasks;
      // if sync interval of 0, this was the first iteration.
      // then I just set next time to infinity and it will never be run again
      if (sync_tasks[i].sync_interval == 0) {
        sync_tasks[i].next_time = size_t(-1);
      }
    }
    active_sync_tasks.clear();
    // figure out what to run next
    for (size_t i = 0;i < sync_tasks.size(); ++i) {
      if (sync_tasks[i].next_time <= num_executed_tasks) {
        active_sync_tasks.push_back(&(sync_tasks[i]));
      }
    }
  }



  /** Checks all machines for termination and sets the termination reason.
      Also returns the number of update tasks completed globally */
  size_t check_global_termination() {
    std::vector<termination_evaluation> termination_test;
    termination_test.resize(rmi.numprocs());
    
    termination_test[rmi.procid()].pending_tasks = num_pending_tasks.value;
    
    size_t numupdates = 0;
    for (size_t i = 0; i < update_counts.size(); ++i) numupdates += update_counts[i];
    termination_test[rmi.procid()].executed_tasks = numupdates;
  
    if (timeout_millis > 0 && ti.current_time_millis() > timeout_millis) {
      termination_test[rmi.procid()].timeout = true;
    }
    
    for (size_t i = rmi.procid(); i < term_functions.size(); i += rmi.numprocs()) {
      if (term_functions[i]()) {
        termination_test[rmi.procid()].terminator = true;
        break;
      }
    }
    termination_test[rmi.procid()].force_stop = force_stop;
    // gather all to 0.
    // machine 0 evaluates termiation
    rmi.gather(termination_test, 0);
    // used to globally evaluate termination
    termination_evaluation aggregate;
    if (rmi.procid() == 0) {
      for (size_t i = 0;i < termination_test.size(); ++i) {
        aggregate.pending_tasks += termination_test[i].pending_tasks;
        aggregate.executed_tasks += termination_test[i].executed_tasks;
        aggregate.terminator |= termination_test[i].terminator;
        aggregate.timeout |= termination_test[i].timeout;
        aggregate.force_stop |= termination_test[i].force_stop;
      }
      
      if (aggregate.pending_tasks == 0) {
        termination_reason = execution_status::TASK_DEPLETION;
      }
      else if (task_budget > 0 && aggregate.executed_tasks >= task_budget) {
        termination_reason = execution_status::TASK_BUDGET_EXCEEDED;
      }
      else if (timeout_millis > 0 && aggregate.timeout) {
        termination_reason = execution_status::TIMEOUT;
      }
      else if (aggregate.terminator) {
        termination_reason = execution_status::TERM_FUNCTION;
      }
      else if (aggregate.force_stop) {
        termination_reason = execution_status::FORCED_ABORT;
      }
    }
    size_t treason = termination_reason;
    // note this is OK because only machine 0 will have the right value for
    // executed_tasks. And everyone is receiving from machine 0
    std::pair<size_t, size_t> reason_and_task(treason, aggregate.executed_tasks);
    rmi.broadcast(reason_and_task, rmi.procid() == 0);
    termination_reason = (execution_status::status_enum)(reason_and_task.first);
    return reason_and_task.second;
  }



  std::vector<vertex_id_type> permuted_ordering;
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_context<engine_type> context(this, &graph, &shared_data);
    timer ti;
    while(1) {
      if (rmi.procid() == 0 && threadid == 0) {
        std::cout << ".";
        std::cout.flush();
      }
      for (size_t i = threadid; i < permuted_ordering.size(); i += ncpus) {
        // otherwise, get the local and globalvid
        vertex_id_type localvid = i;
        vertex_id_type globalvid = graph.localvid_to_globalvid(localvid);

        update_functor_type functor_to_run;
        bool has_functor_to_run = vfunset.test_and_get(localvid, functor_to_run);
        // this is to remote machine
        if (has_functor_to_run) {
          num_pending_tasks.dec();
          
          if (localvid < graph.local_vertices()) {
            // otherwise. run the vertex
            // create the scope
            context.init(globalvid, threadid);
            // run the update function
            functor_to_run(context);
            // check if there are tasks to run
            update_counts[threadid]++;
          }
          else {
            rmi.remote_call(graph.localvid_to_owner(localvid),
                  &distributed_delta_engine<Graph, UpdateFunctor>::schedule_impl,
                  globalvid, functor_to_run);
          }
        }
      }
      if (num_pending_tasks.value == 0) {
        consensus.begin_done_critical_section();
        if (num_pending_tasks.value > 0){
          consensus.end_done_critical_section(false);
        }
        else {
          if (consensus.end_done_critical_section(true)) break;
        }
      }
    }
  }
 
  
  float elapsed_time() const { return lowres_time_seconds() - start_time; }

  
  /** Execute the engine */
  void start() {
    start_time = lowres_time_seconds();
    
    // generate colors then
    // wait for everyone to enter start    
    init_syncs();
    termination_reason = execution_status::UNSET;
    barrier_time = 0.0;
    force_stop = false;
    numsyncs.value = 0;
    num_dist_barriers_called = 0;
    std::fill(update_counts.begin(), update_counts.end(), 0);
    
    // two full barrers to complete flush replies
    rmi.dc().full_barrier();
    rmi.dc().full_barrier();

    // reset indices
    curidx.value = 0;
    ti.start();
    
    // spawn threads
    thread_group thrgrp; 
    for (size_t i = 0;i < ncpus; ++i) {
      size_t aff = use_cpu_affinity ? i : -1;
      thrgrp.launch(boost::bind(
                            &distributed_delta_engine<Graph, UpdateFunctor>::start_thread,
                            this, i), aff);
    }
    thrgrp.join();              
    rmi.barrier();
    

    
    // proc 0 gathers all update counts
    std::vector<size_t> procupdatecounts(rmi.numprocs(), 0);
    procupdatecounts[rmi.procid()] = thisproc_update_counts();
    rmi.gather(procupdatecounts, 0);
    
    std::vector<double> barrier_times(rmi.numprocs(), 0);
    barrier_times[rmi.procid()] = barrier_time;
    rmi.gather(barrier_times, 0);
    // get RMI statistics
    std::map<std::string, size_t> ret = rmi.gather_statistics();

    if (rmi.procid() == 0) {
      total_update_count = 0;
      for(size_t i = 0; i < procupdatecounts.size(); ++i) {
        total_update_count +=  procupdatecounts[i];
      }
      total_barrier_time = 0;
      for(size_t i = 0; i < barrier_times.size(); ++i) {
        total_barrier_time += barrier_times[i];
      }

      total_bytes_sent = ret["total_bytes_sent"];
    }
    
    
    
  }
  
  /**
   * Performs a sync immediately. This function requires that the shared
   * variable already be registered with the engine.
   * and that the engine is not currently running
   * All processes must call simultaneously
   */
  void sync_now(const std::string& key) {
    // TODO
  }
  
  /**
    * \brief Update the engine options.  
    *
    * Setting the engine options will cause all existing state,
    * including scheduled update functors, to be cleared.
    */
  void set_options(const graphlab_options& newopts) {
    opts = newopts;
    rmi.barrier();
  }

  //! \brief Get the current engine options for this engine
  const graphlab_options& get_options() { return opts; } 
  

  
  
  static void print_options_help(std::ostream &out) {
  };


  void stop() {
    force_stop = true;
  }

  
  size_t total_update_count;
  
  size_t get_tasks_done() const {
    return  total_update_count;
  }
  
  double total_barrier_time;
  double get_barrier_time() const {
      return total_barrier_time;
  }
  
  long long int total_bytes_sent;
  long long int get_total_bytes_sent() {
     return total_bytes_sent;
  }
};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED2_DELTA_ENGINE_HPP

