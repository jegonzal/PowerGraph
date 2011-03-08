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
#include <graphlab/distributed2/distributed_glshared_manager.hpp>
#include <graphlab/distributed2/graph/dgraph_scope.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {



/**
All processes must receive the same options at the same time.
i.e. if set_cpu_affinities is called, all processes mus call it at the same time.
This is true for all set_* functions.
*/
template <typename Graph>
class distributed_chromatic_engine : public iengine<Graph> {
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
                                      distributed_chromatic_engine<Graph> > callback_type;
  typedef icallback<Graph> icallback_type;

 private:
  // the local rmi instance
  dc_dist_object<distributed_chromatic_engine<Graph> > rmi;
  
  // the graph we are processing
  Graph &graph;
  
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
  
  size_t randomize_schedule;
  
  /** If dynamic scheduling is used, the number of scheduled tasks */
  atomic<size_t> num_pending_tasks;
  
  
  /** The cause of the last termination condition */
  exec_status termination_reason;

  scope_range::scope_range_enum default_scope_range;
 
  std::vector<std::vector<vertex_id_t> > color_block; // set of localvids in each color
  dense_bitset scheduled_vertices;  // take advantage that local vertices
                                    // are always the first N
  
  update_function_type update_function;
  size_t max_iterations;
  double barrier_time;
  size_t num_dist_barriers_called;
  
  // other optimizations
  bool const_nbr_vertices, const_edges;
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
  

  

 public:
  distributed_chromatic_engine(distributed_control &dc,
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
                            randomize_schedule(0),
                            termination_reason(EXEC_UNSET),
                            scheduled_vertices(graph.owned_vertices().size()),
                            update_function(NULL),
                            max_iterations(0),
                            barrier_time(0.0),
                            const_nbr_vertices(true),
                            const_edges(false),
                            thread_color_barrier(ncpus) { 
    rmi.barrier();
  }
  
  ~distributed_chromatic_engine() {
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
    if (update_function != NULL) assert(update_function == task.function());
    else update_function = task.function();
    
    if (graph.is_owned(task.vertex())) {
      num_pending_tasks.inc(!
              scheduled_vertices.set_bit(graph.globalvid_to_localvid(task.vertex())) 
                          );
      //std::cout << "add task to " << task.vertex() << std::endl;
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
    if (update_function != NULL) assert(update_function == func);
    else update_function = func;
    
    scheduled_vertices.fill();
    num_pending_tasks.value = graph.owned_vertices().size();

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

  
  void generate_color_blocks() {
    // construct for each color, the set of vertices as well as the 
    // number of replicas for that vertex.
    // the number of replicas - 1 is the amount of communication
    // we have to perform to synchronize modifications to that vertex


    std::vector<std::vector<std::pair<size_t, vertex_id_t> > > color_block_and_weight;
    const size_t num_colors(graph.recompute_num_colors());
    // the list of vertices for each color
    color_block_and_weight.resize(num_colors);
    
    foreach(vertex_id_t v, graph.owned_vertices()) {
      color_block_and_weight[graph.get_color(v)].push_back(
                                    std::make_pair(graph.globalvid_to_replicas(v).size(), 
                                              graph.globalvid_to_localvid(v)));
    }
    color_block.clear();
    color_block.resize(num_colors);
    if (randomize_schedule) {
      for (size_t i = 0; i < color_block_and_weight.size(); ++i) {
         std::random_shuffle(color_block_and_weight[i].begin(),
                            color_block_and_weight[i].end());
      }
    }
    else {
      // optimize ordering. Sort in descending order
      // put all those which need a lot of communication in the front
      // to give communication the maximum amount if time possible.
      for (size_t i = 0; i < color_block_and_weight.size(); ++i) {
        std::sort(color_block_and_weight[i].rbegin(),
                  color_block_and_weight[i].rend());
      }
    }

    // insert the sorted vertices into the final color_block
    for (size_t i = 0;i < color_block_and_weight.size(); ++i ) {  
      std::transform(color_block_and_weight[i].begin(),
                    color_block_and_weight[i].end(), 
                    std::back_inserter(color_block[i]),
                    __gnu_cxx::select2nd<std::pair<size_t, vertex_id_t> >());
    }
    
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
      task->thread_intermediate[0] = task->zero;
      for(size_t i = 1; i < task->thread_intermediate.size(); ++i) {
        task->merge_fun(task->mergeval, task->thread_intermediate[i]);
        task->thread_intermediate[i] = task->zero;
      }
      // zero out the intermediate
      task->thread_intermediate.clear();
      task->thread_intermediate.resize(ncpus, sync_tasks[curtask].zero);
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
  size_t check_global_termination(bool check_dynamic_schedule) {
    std::vector<termination_evaluation> termination_test;
    termination_test.resize(rmi.numprocs());
    
    if (check_dynamic_schedule) {
      termination_test[rmi.procid()].pending_tasks = num_pending_tasks.value;
    }

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
        aggregate.pending_tasks += termination_test[i].pending_tasks;
        aggregate.executed_tasks += termination_test[i].executed_tasks;
        aggregate.terminator |= termination_test[i].terminator;
        aggregate.timeout |= termination_test[i].timeout;
        aggregate.force_stop |= termination_test[i].force_stop;
      }
      
      if (check_dynamic_schedule && aggregate.pending_tasks == 0) {
        termination_reason = EXEC_TASK_DEPLETION;
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
 
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_scope<Graph> scope;
    timer ti;

    // loop over iterations
    size_t iter = 0;
    bool usestatic = max_iterations > 0;
    while(1) {
      // if max_iterations is defined, quit
      if (usestatic && iter >= max_iterations) {
        termination_reason = EXEC_TASK_DEPLETION;
        break;
      }
      bool hassynctasks = active_sync_tasks.size() > 0;
      // loop over colors    
      for (size_t c = 0;c < color_block.size(); ++c) {
        // internal loop over vertices in the color
        while(1) {
          // grab a vertex  
          size_t i = curidx.inc_ret_last();  
          // if index out of scope, we are done with this color. break
          if (i >= color_block[c].size()) break;
          // otherwise, get the local and globalvid
          vertex_id_t localvid = color_block[c][i];
          vertex_id_t globalvid = graph.localvid_to_globalvid(color_block[c][i]);
          if (usestatic || scheduled_vertices.clear_bit(localvid)) {
            if (!usestatic) num_pending_tasks.dec();
            // otherwise. run the vertex
            // create the scope
            scope.init(&graph, globalvid);
            // run the update function
            update_function(scope, callback, NULL);
            // check if there are tasks to run
            if (hassynctasks) eval_syncs(globalvid, scope, threadid);
            scope.commit_async_untracked();
            update_counts[threadid]++;
          }
          else {
            // ok this vertex is not scheduled. But if there are syncs
            // to run I will still need to get the scope
            scope.init(&graph, globalvid);
            if (hassynctasks) eval_syncs(globalvid, scope, threadid);
            scope.commit_async_untracked();
          }
        }
        // wait for all threads to synchronize on this color.
        thread_color_barrier.wait();
        curidx.value = 0;
        // full barrier on the color
        // this will complete synchronization of all add tasks as well
        if (threadid == 0) {
          ti.start();
          graph.wait_for_all_async_syncs();
          // TODO! If synchronize() calls were made then this barrier is necessary
          // but the time needed to figure out if a synchronize call is required 
          // could be as long as the barrier itself
          if (const_nbr_vertices == false || const_edges == false)  rmi.dc().barrier();
          rmi.dc().full_barrier();
          num_dist_barriers_called++;
          //std::cout << rmi.procid() << ": Full Barrier at end of color" << std::endl;
          barrier_time += ti.current_time();
        }
        thread_color_barrier.wait();

      }


      sync_end_iteration(threadid);
      thread_color_barrier.wait();
      if (threadid == 0) {
        ti.start();
        //std::cout << rmi.procid() << ": End of all colors" << std::endl;
        size_t numtasksdone = check_global_termination(!usestatic);

        //std::cout << numtasksdone << " tasks done" << std::endl;
        compute_sync_schedule(numtasksdone);
        barrier_time += ti.current_time();
      }
      // all threads must wait for 0
      thread_color_barrier.wait();

      ++iter;
      if (termination_reason != EXEC_UNSET) {
        //std::cout << rmi.procid() << ": Termination Reason: " << termination_reason << std::endl;
        break;
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
      launch_in_new_thread(thrgrp, 
                         boost::bind(
                            &distributed_chromatic_engine<Graph>::start_thread,
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
      total_update_count = 0;
      for(size_t i = 0; i < procupdatecounts.size(); ++i) {
        engine_metrics.add("updatecount", 
                            procupdatecounts[i], INTEGER);
        total_update_count +=  procupdatecounts[i];
      }
      total_barrier_time = 0;
      for(size_t i = 0; i < barrier_times.size(); ++i) {
        engine_metrics.add("barrier_time", 
                            barrier_times[i], TIME);
        total_barrier_time += barrier_times[i];
      }

      engine_metrics.set("termination_reason", 
                        exec_status_as_string(termination_reason));
      engine_metrics.set("dist_barriers_issued", 
                        num_dist_barriers_called, INTEGER);

      engine_metrics.set("num_vertices", graph.num_vertices(), INTEGER);
      engine_metrics.set("num_edges", graph.num_edges(), INTEGER);
      engine_metrics.set("num_syncs", numsyncs.value, INTEGER);
      engine_metrics.set("isdynamic", max_iterations == 0, INTEGER);
      engine_metrics.set("iterations", max_iterations, INTEGER);
      engine_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
      engine_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
      total_bytes_sent = ret["total_bytes_sent"];
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
    opts.get_int_option("max_iterations", max_iterations);
    opts.get_int_option("randomize_schedule", randomize_schedule);

    any uf;
    if (opts.get_any_option("update_function", uf)) {
      update_function = uf.as<update_function_type>();
    }
    rmi.barrier();
  }
  
  void set_randomize_schedule(bool randomize_schedule_) {
    randomize_schedule = randomize_schedule_;
    rmi.barrier();
  }

  
  
  static void print_options_help(std::ostream &out) {
    out << "max_iterations = [integer, default = 0]\n";
    out << "update_function = [update_function_type,"
      "default = set on add_task]\n";
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
  
  // Temp hack.
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

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP
