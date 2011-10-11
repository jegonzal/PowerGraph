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


#ifndef DISTRIBUTED_LOCKING_ENGINE_HPP
#define DISTRIBUTED_LOCKING_ENGINE_HPP

#include <functional>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <ext/functional> // for select1st
#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/mutable_queue.hpp>

#include <graphlab/engine/iengine.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/metrics/metrics.hpp>
#include <graphlab/schedulers/support/redirect_scheduler_callback.hpp>
#include <graphlab/schedulers/support/binary_vertex_task_set.hpp>
#include <graphlab/graph/write_only_disk_atom.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/async_consensus.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>
#include <graphlab/distributed2/graph/dgraph_scope.hpp>
#include <graphlab/distributed2/graph/graph_lock.hpp>
#include <graphlab/distributed2/graph/chandy_misra_lock.hpp>
#include <graphlab/distributed2/graph/distributed_mutex_lock.hpp>
#include <graphlab/distributed2/snapshot_task.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {



/**
All processes must receive the same options at the same time.
i.e. if set_cpu_affinities is called, all processges mus call it at the same time.
This is true for all set_* functions.
*/
template <typename Graph, typename Scheduler >
class distributed_locking_engine:public iengine<Graph> {
 public:
  typedef iengine<Graph> iengine_base;
  typedef typename iengine_base::update_task_type update_task_type;
  typedef typename iengine_base::update_function_type update_function_type;
  typedef typename iengine_base::termination_function_type termination_function_type;
  typedef typename iengine_base::iscope_type iscope_type;
  
  typedef typename iengine_base::sync_function_type sync_function_type;
  typedef typename iengine_base::merge_function_type merge_function_type;

  // unused
  typedef imonitor<Graph> imonitor_type;

  typedef redirect_scheduler_callback<Graph, 
                                      distributed_locking_engine<Graph, Scheduler> > callback_type;
                                      
  typedef icallback<Graph> icallback_type;

  /**
    A special add task redirector for the snapshotter. Provides a 
    serialized add_task interface, as well as access to some snapshotting parameters
  */
  struct snapshot2_scheduler_callback: public icallback<Graph> {
   public:
    distributed_locking_engine<Graph, Scheduler>* eng;
    size_t threadid;
    snapshot2_scheduler_callback(distributed_locking_engine<Graph, Scheduler>* eng,
                                 size_t threadid): eng(eng), threadid(threadid) { }
    
    void add_task(update_task_type task, double priority) {
      eng->snapshot2_add_task(task, priority);
    }
    void add_tasks(const std::vector<vertex_id_t>& vertices, 
                   update_function_type func,
                   double priority) {
      ASSERT_MSG(false, "Unimplemented");
    }
    void force_abort() {
      ASSERT_MSG(false, "Unimplemented");
    }
        
    bool get_snapshot_token(vertex_id_t vid) {
      return eng->snapshot2_tokens.get(eng->graph.globalvid_to_localvid(vid)) == 
                    eng->snapshot2_sense;
    }
    
    bool vertex_modified_since_last_snapshot(vertex_id_t vid) {
      return eng->graph.get_local_store()
                  .vertex_snapshot_req(eng->graph.globalvid_to_localvid(vid));
    }
    
    bool edge_modified_since_last_snapshot(edge_id_t eid) {
      return eng->graph.get_local_store().edge_snapshot_req(eid);
    }
    
    void set_and_synchronize_token(vertex_id_t vid) {
      eng->set_snapshot2_token(vid);
      eng->broadcast_snapshot2_token(vid);
    }
    
    void save_vertex(vertex_id_t vid, const typename Graph::vertex_data_type &vdata) {
      vertex_id_t localvid = eng->graph.globalvid_to_localvid(vid);
      eng->snapshot2_targets[threadid]->add_vertex_with_data(vid,
                                                            eng->graph.localvid_to_source_atom(localvid),
                                                            serialize_to_string(vdata));
                                                            
      eng->graph.get_local_store().set_vertex_snapshot_req(localvid, false);

    }
    
    void save_edge(edge_id_t eid, vertex_id_t srcvid, vertex_id_t targetvid,
                   const typename Graph::edge_data_type &edata) {
      eng->snapshot2_targets[threadid]->add_edge_with_data(srcvid,
                                                          eng->graph.globalvid_to_source_atom(srcvid),
                                                          targetvid,
                                                          eng->graph.globalvid_to_source_atom(targetvid),
                                                          serialize_to_string(edata));
      eng->graph.get_local_store().set_edge_snapshot_req(eid, false);
    }

  };
  
 private:
  // the local rmi instance
  dc_dist_object<distributed_locking_engine<Graph, Scheduler> > rmi;
  
  // the graph we are processing
  Graph &graph;

  // a redirect scheduler call back which 
  callback_type callback;

  
  // The manager will automatically attach to all the glshared variables
  distributed_glshared_manager glshared_manager; 
  
  /** Number of cpus to use */
  size_t ncpus; 

  /** Use processor affinities */
  bool use_cpu_affinity;

  
  /** Track the number of updates */
  std::vector<size_t> update_counts;

  atomic<size_t> numsyncs;

  /** terminators */ 
  std::vector<termination_function_type> term_functions;

  /** The timeout time in millis */
  size_t timeout_millis;
  timer ti;
  
  /// Used to identify when the engine is stopped through stop
  bool force_stop;
  
  /** The total number of tasks that should be executed */
  size_t task_budget;
  
  bool strength_reduction;
  size_t weak_color;
  
  /** The cause of the last termination condition */
  exec_status termination_reason;
  
  /** Parameters for snapshot algorithm 1: synchronous snapshotting */
  size_t snapshot_interval_updates;
  size_t last_snapshot;
  size_t snapshot_number;
  
  /** Parameters for snapshot algorithm 2: asynchronous snapshotting */
  mutex snapshot2_lock;
  size_t snapshot2_interval_updates;
  size_t last_snapshot2;
  size_t snapshot2_number;
  atomic<size_t> snapshot2_remaining_vertices;
  update_function_type snapshot2_update;
  std::vector<write_only_disk_atom*> snapshot2_targets;
  // this is EXTREMELY annoying. I need to tack on an additional bit of information
  // to each vertex. But other than requiring intrusive access to user data,
  // there is no easy way to this. Therefore I need to keep my own bitset and maintain
  // my own data consistency... URGH
  dense_bitset snapshot2_tokens;  
  // the current sense flag. if token == sense, snapshot has been taken
  bool snapshot2_sense;
  


  scope_range::scope_range_enum default_scope_range;
  scope_range::scope_range_enum sync_scope_range;


  /**
   * The set of tasks to have been pulled out of a scheduler
   * and is awaiting execution.
   */
  struct deferred_tasks {
    mutex lock;
    std::deque<update_function_type> updates;
    bool lockrequested;
    deferred_tasks() {
      lockrequested = false;
    }
  };
  
  /**
   * The set of tasks to have been pulled out of a scheduler
   * and is awaiting execution.
   */
  std::vector<deferred_tasks> vertex_deferred_tasks;
  atomic<size_t> num_deferred_tasks;
  size_t max_deferred_tasks;

  multi_blocking_queue<vertex_id_t> ready_vertices;
  
  
  double barrier_time;
  size_t num_dist_barriers_called;

  
  async_consensus consensus;
  
  struct sync_task {
    sync_function_type sync_fun;
    merge_function_type merge_fun;
    distributed_glshared_base::apply_function_type apply_fun;
    size_t sync_interval;
    size_t next_time;
    any zero;
    vertex_id_t rangelow;
    vertex_id_t rangehigh;
    distributed_glshared_base *sharedvariable;
    any mergeval;
    std::vector<any> thread_intermediate;
    sync_task() :
      sync_fun(NULL), merge_fun(NULL), apply_fun(NULL),
      sync_interval(-1),
      next_time(0), rangelow(0), 
      rangehigh(vertex_id_t(-1)), sharedvariable(NULL) { }
  };
  
  /// A list of all registered sync tasks
  std::vector<sync_task> sync_tasks;
  
  /// The list of tasks which are currently being evaluated
  std::vector<sync_task*> active_sync_tasks;

  // scheduler keeps a schedule over localvids
  Scheduler scheduler;
  graph_lock<Graph>* graphlock;
  int chandy_misra;
  /** the number of threads within the main loop when a thread has the
   * intention to leave, it must decrement this before entering the critical
   * section
   */ 
  atomic<size_t> threads_alive;
  
  binary_vertex_task_set<Graph> binary_vertex_tasks;

  
  
  size_t numtasksdone;

  metrics engine_metrics;
  
  size_t total_update_count;
  
  dc_services reduction_services;
 public:
  distributed_locking_engine(distributed_control &dc,
                                    Graph& graph,
                                    size_t ncpus = 1):
                            rmi(dc, this),
                            graph(graph),
                            callback(this),
                            glshared_manager(dc),
                            ncpus( std::max(ncpus, size_t(1)) ),
                            use_cpu_affinity(false),
                            update_counts(std::max(ncpus, size_t(1)), 0),
                            timeout_millis(0),
                            force_stop(false),
                            task_budget(0),
                            strength_reduction(false),
                            weak_color(0),
                            termination_reason(EXEC_UNSET),
                            snapshot_interval_updates(0), 
                            last_snapshot(0),
                            snapshot_number(0),
                            snapshot2_interval_updates(0), 
                            last_snapshot2(0),
                            snapshot2_number(0),
                            snapshot2_update(gl_impl::snapshot_update<Graph, snapshot2_scheduler_callback>),
                            snapshot2_sense(false),
                            default_scope_range(scope_range::EDGE_CONSISTENCY),
                            sync_scope_range(scope_range::VERTEX_CONSISTENCY),
                            vertex_deferred_tasks(graph.owned_vertices().size()),
                            max_deferred_tasks(1000),
                            ready_vertices(ncpus),
                            barrier_time(0.0),
                            consensus(dc, ncpus),
                            scheduler(this, graph, std::max(ncpus, size_t(1))),
                            graphlock(NULL),
                            chandy_misra(0),
                            threads_alive(ncpus),
                            binary_vertex_tasks(graph.local_vertices()),
                            engine_metrics("engine"),
                            total_update_count(0),
                            reduction_services(dc),
                            reduction_barrier(ncpus) { 
    graph.allocate_scope_callbacks();
    dc.barrier();
  }
  
  ~distributed_locking_engine() {
    rmi.barrier();
  }
  
  
  //! Get the number of cpus
  size_t get_ncpus() const { return ncpus; }


  /**
    Turns CPU affinity on or off.
    Must be called by all machines simultaneously
  */
  void set_cpu_affinities(bool value) {
    use_cpu_affinity = value;
    rmi.barrier();
  }


  /**
   * Set the default scope range.  The scope ranges are defined in
   * iscope.hpp
   */
  void set_default_scope(scope_range::scope_range_enum default_scope_range_) {
    default_scope_range = default_scope_range_;
    rmi.barrier();
  }
  
  /**
   * Set the default scope range.  The scope ranges are defined in
   * iscope.hpp
   */
  void set_sync_scope_range(scope_range::scope_range_enum sync_scope_range_) {
    sync_scope_range = sync_scope_range_;
    rmi.barrier();
  }
  
  using iengine<Graph>::exec_status_as_string;   
  

  /**
   * Return the reason why the engine last terminated
   */
  exec_status last_exec_status() const {
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
  void add_terminator(termination_function_type term) {
    term_functions.push_back(term);
    rmi.barrier();
  }


  /**
   * Clear all terminators from the engine
   * Must be called by all machines simultaneously.
   */
  void clear_terminators() {
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
  

  /**
   * \brief Adds an update task with a particular priority.
   * add_task on any vertex can be called by any machine.
   * The call is asynchronous and may not be completed until
   * a full_barrier is issued.
   */
  void add_task(update_task_type task, double priority) {
    if (graph.is_owned(task.vertex())) {
      // translate to local IDs
      task =  update_task_type(graph.globalvid_to_localvid(task.vertex()), task.function());
      ASSERT_LT(task.vertex(), vertex_deferred_tasks.size());
      if (binary_vertex_tasks.add(task)) {
        scheduler.add_task(task, priority);
        if (threads_alive.value < ncpus) {
          consensus.cancel_one();
        }
      }
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(task.vertex()),
                      &distributed_locking_engine<Graph, Scheduler>::add_task,
                      task,
                      priority);
    }
  }

  /**
   * \brief Creates a collection of tasks on all the vertices in
   * 'vertices', and all with the same update function and priority
   * This function is forwarded to the scheduler.
   */
  void add_tasks(const std::vector<vertex_id_t>& vertices,
                         update_function_type func, double priority) {
    // not the most efficient way to do it...
    for (size_t i = 0;i < vertices.size(); ++i) {
      add_task(update_task_type(vertices[i], func), priority);
    }
    if (threads_alive.value < ncpus) {
      consensus.cancel();
    }
  }

  /// \internal
  void add_task_to_all_from_remote(size_t func,
                                  double priority) {
    add_task_to_all_impl(reinterpret_cast<update_function_type>(func), priority);
  }
  
  /// \internal
  void add_task_to_all_impl(update_function_type func,
                            double priority) {
   std::vector<vertex_id_t> perm = graph.owned_vertices();
   random::shuffle(perm);
   for (size_t i = 0;i < perm.size(); ++i) {
      size_t localvid = graph.globalvid_to_localvid(perm[i]);      
      ASSERT_LT(localvid, vertex_deferred_tasks.size());
      if (binary_vertex_tasks.add(update_task_type(localvid, func))) {
        scheduler.add_task(update_task_type(localvid, func), priority);
      }
    }
    if (threads_alive.value < ncpus) {
      consensus.cancel();
    }
  }
 
  /**
   * \brief Creates a collection of tasks on all the vertices in the graph,
   * with the same update function and priority
   * Must be called by all machines simultaneously
   */
  void add_task_to_all(update_function_type func,
                               double priority) {
    add_task_to_all_impl(func,priority);
  }

  /**
    Registers a sync operation.
    Must be called by all machine simultaneously
  */
  void set_sync(glshared_base& shared,
                sync_function_type sync,
                glshared_base::apply_function_type apply,
                const any& zero,
                size_t sync_interval = 0,
                merge_function_type merge = NULL,
                vertex_id_t rangelow = 0,
                vertex_id_t rangehigh = -1) {
    ASSERT_MSG(merge != NULL, "merge is required for the distributed engine");
    sync_task st;
    st.sync_fun = sync;
    st.merge_fun = merge;
    st.apply_fun = apply;
    st.sync_interval = sync_interval;
    st.next_time = 0;
    st.zero = zero;
    st.rangelow = rangelow;
    st.rangehigh = rangehigh;
    st.sharedvariable = dynamic_cast<distributed_glshared_base*>(&shared) ;
    sync_tasks.push_back(st);
    rmi.barrier();
  }

  /************  Actual Execution Engine ****************/

 public: 
  
  struct termination_evaluation {
    size_t executed_tasks;
    bool terminator;
    bool timeout;
    bool force_stop;
    termination_evaluation(): executed_tasks(0),
                              terminator(false),
                              timeout(false),
                              force_stop(false) { }
                              
    void save(oarchive &oarc) const {
      oarc << executed_tasks << terminator << timeout << force_stop;
    }
    
    void load(iarchive &iarc) {
      iarc >> executed_tasks >> terminator >> timeout >> force_stop;
    }
  }; // end of termination evaluation struct



  /**
   * Initialize the sync tasks. Called by start()
   */
  void init_syncs() {
    active_sync_tasks.clear();
    // setup the intermediate values. initialize them to zero
    for (size_t i = 0;i < sync_tasks.size(); ++i) {
      sync_tasks[i].thread_intermediate.clear();
      sync_tasks[i].thread_intermediate.resize(ncpus, sync_tasks[i].zero);
      // everyone runs at the start even if scheduling interval is 0
      active_sync_tasks.push_back(&(sync_tasks[i]));
    }
  }

  /**
   * Called whenever a vertex is executed.  Accumulates the available
   * syncs
   */
  void eval_syncs(vertex_id_t curvertex, iscope_type& scope, size_t threadid) {
    // go through all the active sync tasks
    foreach(sync_task* task, active_sync_tasks) {
      // if in range, sync!
      if (task->rangelow <= curvertex && curvertex <= task->rangehigh) {
        task->sync_fun(scope, task->thread_intermediate[threadid]);
      }
    }
  }

  /** Called at the end of the iteration. Called by all threads after
      a barrier*/
  void sync_end_iteration(size_t threadid) {
    // merge and apply all the syncs. distribute the work among the
    // threads
    for (size_t curtask = threadid; curtask < active_sync_tasks.size(); 
         curtask += ncpus) {
      sync_task* task = active_sync_tasks[curtask];
      task->mergeval = task->thread_intermediate[0];
      task->thread_intermediate[0] = task->zero;
   
      for(size_t i = 1;i < task->thread_intermediate.size(); ++i) {
        task->merge_fun(task->mergeval, task->thread_intermediate[i]);
        task->thread_intermediate[i] = task->zero;
      }
      // for efficiency, lets merge each sync task to the prefered
      // machine
    }
    
    reduction_barrier.wait();

    // one thread of each machine participates in |active_sync_tasks| gathers
    if (threadid == 0) {
      for (size_t i = 0;i < active_sync_tasks.size(); ++i) {
        sync_task* task = active_sync_tasks[i];
        procid_t target = task->sharedvariable->preferred_machine();
        std::vector<any> gathervals(rmi.numprocs());
        gathervals[rmi.procid()] = task->mergeval;
        rmi.gather(gathervals, target);

        // now if I am target I need to do the final merge and apply
        if (target == rmi.procid()) {
          task->mergeval = gathervals[0];
          for (size_t i = 1; i < gathervals.size(); ++i) {
            task->merge_fun(task->mergeval, gathervals[i]);
          }
          // apply!!!
          task->sharedvariable->apply(task->apply_fun, task->mergeval);
          numsyncs.inc();
        }
      }
    }
  }

  /** clears the active sync tasks and figure out what syncs to run next.
      Called by one thread from each machine after sync_end_iteration */
  void compute_sync_schedule(size_t num_executed_tasks) {
    // update the next time variable
    for (size_t i = 0;i < active_sync_tasks.size(); ++i) {
      sync_tasks[i].next_time = num_executed_tasks + sync_tasks[i].sync_interval;
      // if sync interval of 0, this was the first iteration.
      // then I just set next time to infinity and it will never be run again
      if (sync_tasks[i].sync_interval == 0) {
        sync_tasks[i].next_time = size_t(-1);
      }
    }
    active_sync_tasks.clear();
    // figure out what to run next
    for (size_t i = 0;i < sync_tasks.size(); ++i) {
      if (sync_tasks[i].next_time < num_executed_tasks) {
        active_sync_tasks.push_back(&(sync_tasks[i]));
      }
    }
  }

  /** Checks all machines for termination and sets the termination reason.
      Also returns the number of update tasks completed globally */
  size_t check_global_termination() {
    std::vector<termination_evaluation> termination_test;
    termination_test.resize(rmi.numprocs());
    

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
    reduction_services.gather(termination_test, 0);
    // used to globally evaluate termination
    termination_evaluation aggregate;
    if (rmi.procid() == 0) {
      for (size_t i = 0;i < termination_test.size(); ++i) {
        aggregate.executed_tasks += termination_test[i].executed_tasks;
        aggregate.terminator |= termination_test[i].terminator;
        aggregate.timeout |= termination_test[i].timeout;
        aggregate.force_stop |= termination_test[i].force_stop;
      }
      
      if (task_budget > 0 && aggregate.executed_tasks >= task_budget) {
        termination_reason = EXEC_TASK_BUDGET_EXCEEDED;
      }
      else if (timeout_millis > 0 && aggregate.timeout) {
        termination_reason = EXEC_TIMEOUT;
      }
      else if (aggregate.terminator) {
        termination_reason = EXEC_TERM_FUNCTION;
      }
      else if (aggregate.force_stop) {
        termination_reason = EXEC_FORCED_ABORT;
      }
    }
    size_t treason = termination_reason;
    // note this is OK because only machine 0 will have the right value for
    // executed_tasks. And everyone is receiving from machine 0
    std::pair<size_t, size_t> reason_and_task(treason, aggregate.executed_tasks);
    reduction_services.broadcast(reason_and_task, rmi.procid() == 0);
    termination_reason = exec_status(reason_and_task.first);
    return reason_and_task.second;
  }

  
  void wake_up_reducer() {
    reduction_mut.lock();
    reduction_run = true;
    reduction_cond.broadcast();
    reduction_mut.unlock();
  }
  
  mutex reduction_mut;
  conditional reduction_cond;
  bool reduction_stop;
  bool reduction_run;
  barrier reduction_barrier;
  bool reduction_locks_already_acquired;
  
  void reduction_thread(size_t threadid) {
    dgraph_scope<Graph> scope;    
    while(true) {
      reduction_mut.lock();
      while(reduction_stop == false && reduction_run == false) 
        reduction_cond.wait(reduction_mut);
      if (reduction_stop) {
        reduction_mut.unlock();
        break;
      }
      reduction_mut.unlock();
      reduction_barrier.wait();
      reduction_run = false;
      
      // The first thread then evaluates the termination conditions
      // and rescheduels any syncs
      if (threadid == 0) {
        reduction_locks_already_acquired = false;
        //std::cout << rmi.procid() << ": End of all colors" << std::endl;
        numtasksdone = check_global_termination();        
        if (snapshot_interval_updates > 0  && 
            numtasksdone >= last_snapshot + snapshot_interval_updates) {
          reduction_locks_already_acquired = (active_sync_tasks.size() > 0);
          perform_snapshot(active_sync_tasks.size() == 0);
          last_snapshot = numtasksdone;
        }
        
        if (last_snapshot2 > numtasksdone) {
          // this means a snapshot was initialized before
          std::vector<size_t> remainingv(rmi.numprocs());
          remainingv[rmi.procid()] = snapshot2_remaining_vertices.value;
          reduction_services.all_gather(remainingv);
          size_t total_remaining_v = 0;
          for (size_t i = 0;i < remainingv.size(); ++i) {
            total_remaining_v += remainingv[i];
          }
          if (rmi.procid() == 0) {
            logstream(LOG_DEBUG) << "Un-snapshotted vertices: " << total_remaining_v << std::endl;
          }
          if (total_remaining_v == 0) {
            last_snapshot2 = numtasksdone;
          }
        }
        
        if (snapshot2_interval_updates > 0 && 
            numtasksdone >= last_snapshot2 + snapshot2_interval_updates) {
                        
          initialize_snapshot2(!snapshot2_sense, snapshot2_number + 1);
          reduction_services.barrier();
          schedule_snapshot2();
          // set last snapshot2 so it will never be triggered
          last_snapshot2 = std::numeric_limits<size_t>::max() - snapshot2_interval_updates;
        }
        
      }

      reduction_barrier.wait();
      
      if (active_sync_tasks.size() > 0) {
        // Lock all the vertices on this machine that will be operated
        // on by this thread
        if (reduction_locks_already_acquired == false) {
          for (size_t i = threadid; i < graph.owned_vertices().size(); i += ncpus) {
            vertex_deferred_tasks[i].lock.lock();        
          }
          if (threadid == 0) reduction_services.barrier();
        }
        if (threadid == 0) {
          // Synchronize the graph prior to running the sync
          if (sync_scope_range == scope_range::EDGE_CONSISTENCY) {
            if (reduction_locks_already_acquired == false) {
              graph.synchronize_all_edges(true);
              graph.wait_for_all_async_syncs();
            }
          } else if (sync_scope_range == scope_range::FULL_CONSISTENCY) {
            graph.synchronize_all_scopes(true);
            graph.wait_for_all_async_syncs();
          }        
        }
        reduction_barrier.wait();
        // Evaluate all the fold step of all the syncs that should be
        // run at this point
        for (size_t i = threadid; i < graph.owned_vertices().size(); 
             i += ncpus) {
          scope.init(&graph, graph.owned_vertices()[i]);
          eval_syncs(graph.owned_vertices()[i], scope, threadid);
        }
        // Wait for all threads to finish folding
        reduction_barrier.wait();
        // Release the locks
        for (size_t i = threadid;i < graph.owned_vertices().size(); i += ncpus) {
          vertex_deferred_tasks[i].lock.unlock();
        }        
        // wait for all sync threads to finish unlocking
        reduction_barrier.wait();
        // complet the final merge and apply apply operations
        sync_end_iteration(threadid);
      }
      // Wait for all threads to complete the previous syncs (if any)
      reduction_barrier.wait();

      // The first thread then evaluates the termination conditions
      // and rescheduels any syncs
      if (threadid == 0) {
        //std::cout << rmi.procid() << ": End of all colors" << std::endl;
        compute_sync_schedule(numtasksdone);        
        // if I am thread 0 on processor 0, I need to wake up the
        // the main thread which is waiting on a timer
        if (rmi.procid() == 0) {
          std::cout << numtasksdone << " tasks done" << std::endl;
          reduction_complete_signal();
        }
      }
    }
  }
  
  std::string snapshot_filename(size_t snap_number, 
                                size_t threadid, size_t numthreads) {
      std::stringstream strm;
      strm << "snapshot_" << std::setw(3) << std::setfill('0') 
            << snap_number << "_p"
            << std::setw(3) << rmi.procid()
            << ".part_" << threadid+1 << "_of_" << numthreads
            << ".dump";
      return strm.str();
  }

  // if release_locks is false, this will not release the locks 
  // acquired on the vertex_deferred_tasks
  void perform_snapshot(bool release_locks = true) {
    // Get the graph local store which is a local graph containing all
    // the vertices stored on this machine.
    typedef typename Graph::graph_local_store_type 
      graph_local_store_type;
    typedef typename graph_local_store_type::vertex_id_type 
      local_vertex_id_type;
    typedef typename graph_local_store_type::edge_id_type 
      local_edge_id_type;
    graph_local_store_type& graph_local_store = graph.get_local_store();
    // If the snapshot system is no longer needed then return
    // immediately
  // Determine if the proper ammount of time has elapsed
    // Lock all the vertices on this machine
    int nvertex_deferred_tasks = int(vertex_deferred_tasks.size());
    
#pragma omp parallel for
    for (int i = 0; i < nvertex_deferred_tasks; ++i) 
      vertex_deferred_tasks[i].lock.lock();        
    
    reduction_services.barrier();
    // Synchronize the local data in the graph
    graph.synchronize_all_edges(true);
    graph.wait_for_all_async_syncs();
    
    int items_added = 0;
    #pragma omp parallel reduction(+:items_added)
    {
      size_t thread_id = omp_get_thread_num();
      size_t numthreads = omp_get_num_threads();
      // Go ahead and open files for the snapshot
      const std::string filename = snapshot_filename(snapshot_number,thread_id,numthreads);
      write_only_disk_atom atom(filename, rmi.procid(), true);
      // Wait for the syncs to finish

      // Save the local vertex data. Only take owned vertices
      for(local_vertex_id_type localvid = thread_id; 
          localvid < graph_local_store.num_vertices(); localvid+=numthreads) {
        if(graph.localvid_is_ghost(localvid) == false &&
          graph_local_store.vertex_snapshot_req(localvid)) {
          atom.add_vertex_with_data(graph.localvid_to_globalvid(localvid),
                                    graph.localvid_to_source_atom(localvid),
                                    serialize_to_string(graph_local_store.vertex_data(localvid)));
          graph_local_store.set_vertex_snapshot_req(localvid, false);
          ++items_added;
        }
      } // end of for loop over local vids
      
      // Save the local edge data. Only take owned edges (target is owned)
      for(local_edge_id_type localeid = thread_id; 
          localeid < graph_local_store.num_edges(); localeid +=numthreads) {
        if(graph.localvid_is_ghost(graph_local_store.target(localeid)) == false &&
          graph_local_store.edge_snapshot_req(localeid)) {
          vertex_id_t localsource = graph_local_store.source(localeid);
          vertex_id_t localtarget = graph_local_store.target(localeid);
          atom.add_edge_with_data(graph.localvid_to_globalvid(localsource),
                                  graph.localvid_to_source_atom(localsource),
                                  graph.localvid_to_globalvid(localtarget),
                                  graph.localvid_to_source_atom(localtarget),
                                  serialize_to_string(graph_local_store.edge_data(localeid)));
          graph_local_store.set_edge_snapshot_req(localeid, false);
          ++items_added;
        }
      } // end of for loop over local eids
      atom.synchronize();
    }
    std::cout << "Snapshot "<< snapshot_number << " Items added in snapshot: " << items_added << std::endl;
    ++snapshot_number;
    // free all the vertices on this machine
    if (release_locks) {
      for (size_t i = 0; i < vertex_deferred_tasks.size(); ++i) 
        vertex_deferred_tasks[i].lock.unlock();      
    }
  } // end of snapshot thread


/**
  Asynchronous Snapshot implementation.
  1: Node 0 must call initialize_snapshot2 on all nodes and wait for a reply
  2: schedule_snapshot2() is then called on all nodes
*/
  void initialize_snapshot2(bool sense, size_t snapshot_number) {  
    snapshot2_lock.lock();
    for (size_t i = 0;i < snapshot2_targets.size(); ++i) {
      snapshot2_targets[i] = new write_only_disk_atom(
                                      snapshot_filename(snapshot_number, i, ncpus), 
                                      rmi.procid(), 
                                      true);
    }
    snapshot2_number = snapshot_number;
    snapshot2_remaining_vertices.value = graph.owned_vertices().size();
    // flip the sense last
    snapshot2_lock.unlock();
    snapshot2_sense = sense;

    
  }
  void schedule_snapshot2() {
      // schedule all vertices with a fake task
    add_task_to_all_impl(snapshot2_update, 
                         100.0);
  }  
  
  void broadcast_snapshot2_token(vertex_id_t globalvid) {
    const std::vector<procid_t>& replicas = graph.globalvid_to_replicas(globalvid);
    unsigned char prevkey = rmi.dc().set_sequentialization_key((globalvid % 254) + 1); 
    for (size_t i = 0; i < replicas.size(); ++i) {
      if (replicas[i] != rmi.procid()) {
        rmi.remote_call(replicas[i],
                        &distributed_locking_engine<Graph, Scheduler>::set_snapshot2_token,
                        globalvid);
      }
    }
    rmi.dc().set_sequentialization_key(prevkey);
  }
  
  void set_snapshot2_token(vertex_id_t globalvid) {
    snapshot2_tokens.set(graph.globalvid_to_localvid(globalvid), snapshot2_sense);
  }
  
  // a special version of add_task that ensures that the task is injected before
  // the update function releases.
  void snapshot2_add_task(update_task_type task, double priority) {
    if (graph.is_owned(task.vertex())) {
      // translate to local IDs
      task =  update_task_type(graph.globalvid_to_localvid(task.vertex()), task.function());
      ASSERT_LT(task.vertex(), vertex_deferred_tasks.size());
      if (binary_vertex_tasks.add(task)) {
        scheduler.add_task(task, priority);
        if (threads_alive.value < ncpus) {
          consensus.cancel_one();
        }
      }
    }
    else {
      unsigned char prevkey = rmi.dc().set_sequentialization_key((task.vertex() % 254) + 1);
      rmi.remote_call(graph.globalvid_to_owner(task.vertex()),
                      &distributed_locking_engine<Graph, Scheduler>::snapshot2_add_task,
                      task,
                      priority);
      rmi.dc().set_sequentialization_key(prevkey);
    }
  }
  

  /** Vertex i is ready. put it into the ready vertices set */
  void vertex_is_ready(vertex_id_t v) {
    //logstream(LOG_DEBUG) << "Enqueue: " << v << std::endl;
    ready_vertices.enqueue(graph.globalvid_to_localvid(v));
  }

  bool try_to_quit(size_t threadid, 
                   sched_status::status_enum& stat, 
                   update_task_type &task) {
    threads_alive.dec();
    consensus.begin_done_critical_section();
    stat = scheduler.get_next_task(threadid, task);
    if (stat == sched_status::EMPTY) {
      bool ret = consensus.end_done_critical_section(true);
      threads_alive.inc();
      return ret;
    } else {
      consensus.end_done_critical_section(false);
      threads_alive.inc();
      return false;
    }
  } // end of try to quit
  
  
  /**
   * Executed by a thread.
   *  - Begin deferred task
   *  - check in the "available" task set for a job to do
   *  - take the vertex, run all the jobs on it
   *  - check local termination
   *  - loooop
   */
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_scope<Graph> scope;
    update_task_type task;
    snapshot2_scheduler_callback snapshot2_callback(this, threadid);
    
    boost::function<void(vertex_id_t)> handler = 
      boost::bind(&distributed_locking_engine<Graph, Scheduler>::vertex_is_ready, this, _1);
    while(1) {
      if (termination_reason != EXEC_UNSET) {
        consensus.force_done();
        break;
      }
      // pick up a deferred task 
      if (num_deferred_tasks.value < max_deferred_tasks) {
        sched_status::status_enum stat = scheduler.get_next_task(threadid, task);
        // if there is nothing in the queue, and there are no deferred tasks to run
        // lets try to quit
        if (stat == sched_status::EMPTY && num_deferred_tasks.value == 0) {
          bool ret = try_to_quit(threadid, stat, task);
          if (ret == true) break;
          if (ret == false && stat == sched_status::EMPTY) continue;
        }
        
        //if scheduler game me a task
        if (stat != sched_status::EMPTY) {
          //added a deffered task
          num_deferred_tasks.inc();
          // translate the task back to globalids
          vertex_id_t globalvid = graph.localvid_to_globalvid(task.vertex());
  
          // acquire the lock
          ASSERT_LT(task.vertex(), vertex_deferred_tasks.size());
          vertex_deferred_tasks[task.vertex()].lock.lock();
          // insert the task
          vertex_deferred_tasks[task.vertex()].updates.push_back(task.function());
          // if a lock was not requested. request for it
          if (vertex_deferred_tasks[task.vertex()].lockrequested == false) {
            vertex_deferred_tasks[task.vertex()].lockrequested = true;
            if (strength_reduction == false || graph.color(globalvid) != weak_color) {
              graphlock->scope_request(globalvid, handler, default_scope_range);
            }
            else {
              graphlock->scope_request(globalvid, handler, scope_range::VERTEX_CONSISTENCY);
            }
          }
          vertex_deferred_tasks[task.vertex()].lock.unlock();
        }
      }
      // pick up a job to do
      std::pair<vertex_id_t, bool> job = ready_vertices.try_dequeue(threadid);
      
      while (job.second) {
        // lets do it
        // curv is a localvid
        vertex_id_t curv = job.first;
        vertex_id_t globalvid = graph.localvid_to_globalvid(curv);
        //logstream(LOG_DEBUG) << "Dequeue: " << globalvid << std::endl;

        ASSERT_LT(curv, vertex_deferred_tasks.size());

        vertex_deferred_tasks[curv].lock.lock();
        // check for snapshot task first
        if (snapshot2_interval_updates > 0 && 
            binary_vertex_tasks.get(update_task_type(curv, snapshot2_update))) {
          // am I scheduled?
          if (snapshot2_tokens.get(curv) != snapshot2_sense) {
            scope.init(&graph, globalvid);
            ASSERT_TRUE(snapshot2_targets[threadid] != NULL);
            snapshot2_update(scope, snapshot2_callback);
            ASSERT_EQ(snapshot2_tokens.get(curv), snapshot2_sense);
            if (snapshot2_remaining_vertices.dec() == 0) {
              snapshot2_lock.lock();
              for (size_t i = 0;i < snapshot2_targets.size(); ++i) {
                delete snapshot2_targets[i];
                snapshot2_targets[i] = NULL;
              }
              snapshot2_lock.unlock();
              logger(LOG_DEBUG, "Local Snapshot complete!");
            }
          }
          binary_vertex_tasks.remove(update_task_type(curv, snapshot2_update));
        }
        while (!vertex_deferred_tasks[curv].updates.empty()) {
          update_function_type ut = vertex_deferred_tasks[curv].updates.front();
          vertex_deferred_tasks[curv].updates.pop_front();
          //  vertex_deferred_tasks[curv].lock.unlock();

          scope.init(&graph, globalvid);
          
          binary_vertex_tasks.remove(update_task_type(curv, ut));
          if (ut != snapshot2_update) {
            // run the update function
            ut(scope, callback);          
            update_counts[threadid]++;
          }
          //vertex_deferred_tasks[curv].lock.lock();

          num_deferred_tasks.dec();
        }
        vertex_deferred_tasks[curv].lockrequested = false;
        
        if (strength_reduction == false || graph.color(globalvid) != weak_color) {
          graphlock->scope_unlock(globalvid, default_scope_range);
        }
        else {
          graphlock->scope_unlock(globalvid, scope_range::VERTEX_CONSISTENCY);
        }
        vertex_deferred_tasks[curv].lock.unlock();
	
	job = ready_vertices.try_dequeue(threadid);
      }
      sched_yield();
    }
  }
  
  /*
  These variables protect the reduction completion 
  flag on processor 0
  */
  mutex reduction_started_mut;
  conditional reduction_started_cond;
  bool proc0_reduction_started;
  
  size_t get_tasks_done() const {
    return numtasksdone;
  }
  
  void reduction_complete_signal() {
    reduction_started_mut.lock();
    proc0_reduction_started = false;
    reduction_started_cond.signal();
    reduction_started_mut.unlock();
  }
  
  /** Execute the engine */
  void start() {
    if (chandy_misra) {
      graphlock = new chandy_misra_lock<Graph>(rmi.dc(), graph, true);
    }
    else {
      graphlock = new distributed_mutex_lock<Graph>(rmi.dc(), graph, true);
    }
    // generate colors then
    // wait for everyone to enter start    
    if (sync_scope_range == scope_range::FULL_CONSISTENCY &&
        sync_scope_range == scope_range::EDGE_CONSISTENCY) {
      logstream(LOG_WARNING) << "Currently sync operations which require edge or full consistency incurs massive synchronization penalties" << std::endl;
    }
    init_syncs();
    termination_reason = EXEC_UNSET;
    barrier_time = 0.0;
    num_deferred_tasks.value = 0;
    force_stop = false;
    numsyncs.value = 0;
    num_dist_barriers_called = 0;
    numtasksdone = 0;
    last_snapshot = 0;
    snapshot_number = 0;
    last_snapshot2 = 0;
    snapshot2_number = 0;
    snapshot2_remaining_vertices.value = 0;
    snapshot2_sense = false;
    // using snapshot 2
    if (snapshot2_interval_updates > 0) {
      snapshot2_tokens.resize(graph.get_local_store().num_vertices());
      snapshot2_targets.resize(ncpus, NULL);
      snapshot2_tokens.clear();
    }
    reduction_stop = false; 
    reduction_run = false;
    threads_alive.value = ncpus;
    proc0_reduction_started = false;
    std::fill(update_counts.begin(), update_counts.end(), 0);
    rmi.dc().full_barrier();
    // reset indices
    ti.start();
    // spawn worker threads
    thread_group thrgrp; 
    for (size_t i = 0;i < ncpus; ++i) {
      size_t aff = use_cpu_affinity ? i : -1;
      thrgrp.launch(boost::bind(&distributed_locking_engine::start_thread,
                                this, i), aff);
    }
    // spawn reduction threads
    thread_group thrgrp_reduction; 
    for (size_t i = 0;i < ncpus; ++i) {
      size_t aff = use_cpu_affinity ? i : -1;
      thrgrp_reduction.
        launch(boost::bind(&distributed_locking_engine::reduction_thread,
                           this, i), aff);
    }
    
    rmi.barrier();
    
    std::map<double, size_t> upspertime;
    timer ti;
    ti.start();
    if (rmi.procid() == 0) {
      while(consensus.done_noblock() == false) {
        reduction_started_mut.lock();
        proc0_reduction_started = true;
        reduction_started_mut.unlock();
        for (size_t i = 1;i < rmi.numprocs(); ++i) {
          rmi.remote_call(i,
                          &distributed_locking_engine<Graph, Scheduler>::wake_up_reducer);
        }
        wake_up_reducer();
        
        reduction_started_mut.lock();
        while (proc0_reduction_started) reduction_started_cond.wait(reduction_started_mut);
        upspertime[ti.current_time()] = numtasksdone;
        reduction_started_mut.unlock();
        sleep(1);
      }
    }
    thrgrp.join();          
    reduction_mut.lock();
    reduction_stop = true;
    reduction_cond.broadcast();
    reduction_mut.unlock();
    thrgrp_reduction.join();    
    rmi.barrier();
    
    if (termination_reason == EXEC_UNSET) termination_reason = EXEC_TASK_DEPLETION;


    
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
      engine_metrics.add("runtime",
                        ti.current_time(), TIME);
      for(size_t i = 0; i < procupdatecounts.size(); ++i) {
        engine_metrics.add_vector_entry("updatecount", i,  procupdatecounts[i]);
        total_update_count += procupdatecounts[i];
      }
      for(size_t i = 0; i < barrier_times.size(); ++i) {
        engine_metrics.add_vector_entry("barrier_time", i, barrier_times[i]);
      }

      std::map<double, size_t>::const_iterator iter = upspertime.begin();
      while(iter != upspertime.end()) {
        engine_metrics.add_to_vector("updatecount_vector_t", iter->first);
        engine_metrics.add_to_vector("updatecount_vector_v", iter->second);        
        ++iter;
      }
      engine_metrics.set("termination_reason", 
                        exec_status_as_string(termination_reason));
      engine_metrics.add("dist_barriers_issued",
                        num_dist_barriers_called, INTEGER);

      engine_metrics.set("num_vertices", graph.num_vertices(), INTEGER);
      engine_metrics.set("num_edges", graph.num_edges(), INTEGER);
      engine_metrics.add("num_syncs", numsyncs.value, INTEGER);
      engine_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
      engine_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
      total_bytes_sent = ret["total_bytes_sent"];
    }
    
    
    threads_alive.value = ncpus;
    delete graphlock;
    graphlock = NULL;
  }
  
  /**
   * Performs a sync immediately. This function requires that the shared
   * variable already be registered with the engine.
   * and that the engine is not currently running
   */
  void sync_now(glshared_base& shared) {
    // TODO
  }
  
    /** \brief Update the scheduler options.  */
  void set_engine_options(const scheduler_options& opts) {
    opts.get_int_option("max_deferred_tasks_per_node", max_deferred_tasks);
    opts.get_int_option("chandy_misra", chandy_misra);
    opts.get_int_option("snapshot_interval", snapshot_interval_updates);
    opts.get_int_option("snapshot2_interval", snapshot2_interval_updates); 
    rmi.barrier();
  }
  
  void set_scheduler_options(const scheduler_options& opts) {
    scheduler.set_options(opts);
  }
  
  void set_strength_reduction(bool strength_reduction_, size_t weak_color_ = 0) {
    strength_reduction = strength_reduction_;
    // TODO: More intelligent picking of the color to weaken
    weak_color = weak_color_;
    rmi.barrier();
  }
  
  
  void set_max_deferred(size_t max_deferred) {
    max_deferred_tasks = max_deferred;
    rmi.barrier();
  }
  
  static void print_options_help(std::ostream &out) {
    out << "max_deferred_tasks_per_node = [integer, default = 1000]\n";
    out << "chandy_misra = [int, default = 0, If non-zero, uses the chandy misra locking method. Only supports edge scopes]\n";
    out << "snapshot_interval = [integer, default = 0, If non-zero, snapshots approximately this many updates]\n";
    out << "snapshot2_interval = [integer, default = 0, Fully asynchronous snapshotting. If non-zero, snapshots approximately this many updates]\n";
  };


  void stop() {
    force_stop = true;
  }

  void register_monitor(imonitor_type* listener) {
    logger(LOG_FATAL, "distributed engine does not support register monitor");
  }   
  
  
  long long int total_bytes_sent;
  long long int get_total_bytes_sent() {
     return total_bytes_sent;
  }


  metrics get_metrics() {
    return engine_metrics;
  }


  void reset_metrics() {
    engine_metrics.clear();
  }

};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED_LOCKING_ENGINE_HPP

