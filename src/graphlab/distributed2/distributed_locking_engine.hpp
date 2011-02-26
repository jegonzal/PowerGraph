#ifndef DISTRIBUTED_CHROMATIC_ENGINE_HPP
#define DISTRIBUTED_CHROMATIC_ENGINE_HPP

#include <functional>
#include <algorithm>
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

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/async_consensus.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>
#include <graphlab/distributed2/graph/dgraph_scope.hpp>
#include <graphlab/distributed2/graph/graph_lock.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {



/**
All processes must receive the same options at the same time.
i.e. if set_cpu_affinities is called, all processes mus call it at the same time.
This is true for all set_* functions.
*/
template <typename Graph, typename Scheduler>
class distributed_locking_engine:public iengine<Graph> {
 public:
  typedef iengine<Graph> iengine_base;
  typedef typename iengine_base::update_task_type update_task_type;
  typedef typename iengine_base::update_function_type update_function_type;
  typedef typename iengine_base::termination_function_type termination_function_type;
  typedef typename iengine_base::iscope_type iscope_type;
  typedef typename iengine_base::ishared_data_type ishared_data_type;
  typedef typename iengine_base::ishared_data_manager_type ishared_data_manager_type;
  
  typedef typename iengine_base::sync_function_type sync_function_type;
  typedef typename iengine_base::merge_function_type merge_function_type;

  // unused
  typedef imonitor<Graph> imonitor_type;

  typedef redirect_scheduler_callback<Graph, 
                                      distributed_locking_engine<Graph> > callback_type;
  typedef icallback<Graph> icallback_type;
private:

 private:
  // the local rmi instance
  dc_dist_object<distributed_locking_engine<Graph> > rmi;
  
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

  /** Use schedule yielding when waiting on the scheduler*/
  bool use_sched_yield;
 
  
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
  
 
  
  /** The cause of the last termination condition */
  exec_status termination_reason;

  scope_range::scope_range_enum default_scope_range;


  /**
   * The set of tasks to have been pulled out of a scheduler
   * and is awaiting execution.
   */
  struct deferred_tasks {
    spinlock lock;
    std::deque<update_function_type> updates;
    bool lockrequested;
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
  
  // other optimizations
  bool const_nbr_vertices, const_edges;

  async_consensus consensus;
  
  struct sync_task {
    sync_function_type sync_fun;
    merge_function_type merge_fun;
    distributed_glshared_base::apply_function_type apply_fun;
    size_t sync_interval;
    size_t next_time;
    any zero;
    size_t rangelow;
    size_t rangehigh;
    distributed_glshared_base *sharedvariable;
    any mergeval;
    std::vector<any> thread_intermediate;
    sync_task() :
      sync_fun(NULL), merge_fun(NULL), apply_fun(NULL),
      sync_interval(-1),
      next_time(0), rangelow(0), 
      rangehigh(size_t(-1)), sharedvariable(NULL) { }
  };
  
  /// A list of all registered sync tasks
  std::vector<sync_task> sync_tasks;
  
  /// The list of tasks which are currently being evaluated
  std::vector<sync_task*> active_sync_tasks;

  Scheduler scheduler;
  graph_lock<Graph> graphlock;

  /** the number of threads within the main loop when a thread has the
   * intention to leave, it must decrement this before entering the critical
   * section
   */ 
  atomic<size_t> threads_alive;
  
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
                            use_sched_yield(true),
                            update_counts(std::max(ncpus, size_t(1)), 0),
                            timeout_millis(0),
                            force_stop(false),
                            task_budget(0),
                            termination_reason(EXEC_UNSET),
                            vertex_deferred_tasks(graph.owned_vertices().size()),
                            max_deferred_tasks(-1),
                            ready_vertices(ncpus),
                            barrier_time(0.0),
                            const_nbr_vertices(true),
                            const_edges(false),
                            consensus(dc, this),
                            scheduler(this, graph, std::max(ncpus, size_t(1))),
                            graphlock(dc, graph),
                            thread_color_barrier(ncpus) { 
    rmi.barrier();
  }
  
  ~distributed_locking_engine() {
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
  void set_default_scope(scope_range::scope_range_enum default_scope_range_) {
    default_scope_range = default_scope_range_;
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
  size_t last_update_count() const {
    size_t sum = 0;
    for(size_t i = 0; i < update_counts.size(); ++i)
      sum += update_counts[i];
    return sum;
  } // end of last_update_count


  /**
   * Add a terminator to the engine.
   */
  void add_terminator(termination_function_type term) {
    term_functions.push_back(term);
    rmi.barrier();
  }


  /**
   * Clear all terminators from the engine
   */
  void clear_terminators() {
    term_functions.clear();
    rmi.barrier();
  }



  /**
   * Timeout. Default - no timeout. 
   */
  void set_timeout(size_t timeout_seconds = 0) {
    timeout_millis = timeout_seconds * 1000;
    rmi.barrier();
  }
  
  /**
   * Task budget - max number of tasks to allow
   */
  void set_task_budget(size_t max_tasks) {
    task_budget = max_tasks;
    rmi.barrier();
  }
  

  /**
   * \brief Adds an update task with a particular priority.
   * This function is forwarded to the scheduler.
   */
  void add_task(update_task_type task, double priority) {
    if (graph.is_owned(task.vertex())) {
      scheduler.add_task(task, priority);
      if (threads_alive.value < ncpus) {
        consensus.cancel_one();
      }
    }
    else {
      rmi.remote_call(graph.globalvid_to_owner(task.vertex()),
                      &distributed_chromatic_engine<Graph>::add_task,
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
  }


  void add_task_to_all_from_remote(size_t func,
                                  double priority) {
    add_task_to_all_impl(reinterpret_cast<update_function_type>(func), priority);
  }
  
  void add_task_to_all_impl(update_function_type func,
                            double priority) {
    for (size_t i = 0;i < graph.owned_vertices().size(); ++i) {
      scheduler.add_task(update_task_type(graph.owned_vertices()[i], func), priority);
  }
 
  /**
   * \brief Creates a collection of tasks on all the vertices in the graph,
   * with the same update function and priority
   * This function is forwarded to the scheduler.
   */
  void add_task_to_all(update_function_type func,
                               double priority) {
    add_task_to_all_impl(func,priority);
    for (size_t i = 0;i < rmi.numprocs(); ++i) {
      if (i != rmi.procid()) {
        rmi.remote_call(i,
                        &distributed_chromatic_engine<Graph>::add_task_to_all_from_remote,
                        reinterpret_cast<size_t>(func),
                        priority);
      }
    }
  }

  void set_sync(glshared_base& shared,
                sync_function_type sync,
                glshared_base::apply_function_type apply,
                const any& zero,
                size_t sync_interval = 0,
                merge_function_type merge = NULL,
                size_t rangelow = 0,
                size_t rangehigh = -1) {
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
 private:

  atomic<size_t> curidx;
  barrier thread_color_barrier;
 public: 
  
  struct termination_evaluation{
    size_t executed_tasks;
    bool terminator;
    bool timeout;
    bool force_stop;
    termination_evaluation(): executed_tasks(0),
                              terminator(false),
                              timeout(false),
                              force_stop(false) { }
                              
    void save(oarchive &oarc) const {
      oarc << executed_tasks
           << terminator
           << timeout
           << force_stop;
    }
    
    void load(iarchive &iarc) {
      iarc >> executed_tasks
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
    // setup the intermediate values. initialize them to zero
    for (size_t i = 0;i < sync_tasks.size(); ++i) {
      sync_tasks[i].thread_intermediate.clear();
      sync_tasks[i].thread_intermediate.resize(ncpus, sync_tasks[i].zero);
      // everyone runs at the start even if scheduling interval is 0
      active_sync_tasks.push_back(&(sync_tasks[i]));
    }
  }

  /**
   * Called whenever a vertex is executed.
   * Accumulates the available syncs
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

  /** Called at the end of the iteration. Called by all threads after a barrier*/
  void sync_end_iteration(size_t threadid) {
    // merge and apply all the syncs. distribute the work among the threads
    for (size_t curtask = threadid; curtask < active_sync_tasks.size(); curtask += ncpus) {
      sync_task* task = active_sync_tasks[curtask];
      task->mergeval = task->thread_intermediate[0];
      for(size_t i = 1;i < task->thread_intermediate.size(); ++i) {
        task->merge_fun(task->mergeval, task->thread_intermediate[i]);
      }
      // for efficiency, lets merge each sync task to the prefered machine
    }
    
    thread_color_barrier.wait();

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
      if (term_functions[i](NULL)) {
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
        aggregate.executed_tasks += termination_test[i].executed_tasks;
        aggregate.terminator |= termination_test[i].terminator;
        aggregate.timeout |= termination_test[i].timeout;
        aggregate.force_stop |= termination_test[i].force_stop;
      }
      
      else if (task_budget > 0 && aggregate.executed_tasks >= task_budget) {
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
    rmi.broadcast(reason_and_task, rmi.procid() == 0);
    termination_reason = exec_status(reason_and_task.first);
    return reason_and_task.second;
  }

  /** Vertex i is ready. put it into the ready vertices set */
  void vertex_is_ready(vertex_id_t v) {
    
  }

  bool try_to_quit(size_t threadid) {
    //
    threads_alive.dec();
    consensus.begin_done_critical_section();
    //check the scheduler again
    sched_status::status_enum stat = scheduler->get_next_task(cpuid, task);
    if (stat == sched_status::EMPTY) {
      return consensus.end_done_critical_section(true);
    }
    else {
      consensus.end_done_critical_section(false);
      threads_alive.inc();
      return false;
    }
      
  }
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
    update_task task;
    
    boost::function<void(vertex_id_t)> handler = boost::bind(&distributed_locking_engine<Graph, Scheduler>::vertex_is_ready, this, _1);
    while(1) {
      // pick up a deferred task 
      if (num_deferred_tasks.value < max_deferred_tasks) {
        sched_status::status_enum stat = scheduler->get_next_task(cpuid, task);
        if (stat == sched_status::EMPTY && num_deferred_tasks.value == 0) {
          if (try_to_quit(threadid)) {
            break;
          }
        }
        //
        num_deferred_tasks.inc();
        graph_lock.scope_request(task.vertex(), handler, default_scope_range);
      }

      // pick up a job to do
      std::pair<vertex_id_t, bool> job = ready_vertices.try_dequeue(threadid);
      if (job.second) {
        // lets do it
        size_t curv = job.first;
        vertex_deferred_tasks[curv].lock.lock();
        while (!vertex_deferred_tasks[curv].updates.empty()) {
          update_function_type ut = vertex_deferred_tasks[curv].updates.front();
          vertex_deferred_tasks[curv].updates.pop_front();
          vertex_deferred_tasks[curv].lock.unlock();

          scope.init(&graph, curv);
          // run the update function
          update_function(scope, callback, NULL);
          // check if there are tasks to run
          scope.commit_async_untracked();
          update_counts[threadid]++;
        }
        foreach((deferred_tasks& dt vertex_deferred_tasks[job.first]
      }
    }
  }
  
  void set_const_edges(bool const_edges_ = true) {
    const_edges = const_edges_;
  }

  void set_const_nbr_vertices(bool const_nbr_vertices_ = true) {
    const_nbr_vertices = const_nbr_vertices_;
  }

  
  /** Execute the engine */
  void start() {
    assert(update_function != NULL);
    // generate colors then
    // wait for everyone to enter start    
    generate_color_blocks();
    init_syncs();
    termination_reason = EXEC_UNSET;
    barrier_time = 0.0;
    num_deferred_tasks.value = 0;
    force_stop = false;
    numsyncs.value = 0;
    num_dist_barriers_called = 0;
    threads_alive = ncpus;
    
    std::fill(update_counts.begin(), update_counts.end(), 0);
    rmi.dc().full_barrier();
    // reset indices
    curidx.value = 0;
    ti.start();
    // spawn threads
    thread_group thrgrp; 
    for (size_t i = 0;i < ncpus; ++i) {
      size_t aff = use_cpu_affinity ? i : -1;
      launch_in_new_thread(thrgrp, 
                         boost::bind(
                            &distributed_locking_engine<Graph>::start_thread,
                            this, i), aff);
    }
    
    thrgrp.join();              
    rmi.barrier();
    

    
    // proc 0 gathers all update counts
    std::vector<size_t> procupdatecounts(rmi.numprocs(), 0);
    procupdatecounts[rmi.procid()] = last_update_count();
    rmi.gather(procupdatecounts, 0);
    
    std::vector<double> barrier_times(rmi.numprocs(), 0);
    barrier_times[rmi.procid()] = barrier_time;
    rmi.gather(barrier_times, 0);
    // get RMI statistics
    std::map<std::string, size_t> ret = rmi.gather_statistics();

    if (rmi.procid() == 0) {
      metrics& engine_metrics = metrics::create_metrics_instance("engine", true);
      engine_metrics.set("runtime", 
                        ti.current_time(), TIME);
      for(size_t i = 0; i < procupdatecounts.size(); ++i) {
        engine_metrics.add("updatecount", 
                            procupdatecounts[i], INTEGER);
      }
      for(size_t i = 0; i < barrier_times.size(); ++i) {
        engine_metrics.add("barrier_time", 
                            barrier_times[i], TIME);
      }

      engine_metrics.set("termination_reason", 
                        exec_status_as_string(termination_reason));
      engine_metrics.set("dist_barriers_issued", 
                        num_dist_barriers_called, INTEGER);

      engine_metrics.set("num_vertices", graph.num_vertices(), INTEGER);
      engine_metrics.set("num_edges", graph.num_edges(), INTEGER);
      engine_metrics.set("num_syncs", numsyncs.value, INTEGER);
      engine_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
      engine_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
  
    }
    
    
    
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
  void set_scheduler_options(const scheduler_options& opts) {
    opts.get_int_option("max_deferred_tasks_per_node", max_deferred_tasks);
    rmi.barrier();
  }
  
  static void print_options_help(std::ostream &out) {
    out << "max_deferred_tasks_per_node = [integer, default = unsigned word max]\n";
  };

  void set_shared_data_manager(ishared_data_manager_type* manager) { 
    logger(LOG_FATAL, "distributed engine does not support set shared data manager");
  }

  void stop() {
    force_stop = true;
  }

  void register_monitor(imonitor_type* listener) {
    logger(LOG_FATAL, "distributed engine does not support register monitor");
  }   
};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP
