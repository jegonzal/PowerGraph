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
class distributed_chromatic_engine:public iengine<Graph> {
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
  
  struct sync_task {
    sync_function_type sync_fun;
    merge_function_type merge_fun;
    distributed_glshared_base::apply_function_type apply_fun;
    size_t sync_interval;
    size_t next_time;
    any zero;
    mutex lock;
    size_t rangelow;
    size_t rangehigh;
    distributed_glshared_base *sharedvariable;
    sync_task() :
      sync_fun(NULL), merge_fun(NULL), apply_fun(NULL),
      sync_interval(-1),
      next_time(0), rangelow(0), 
      rangehigh(size_t(-1)), sharedvariable(NULL) { }
  };
  
  /// A list of all registered sync tasks
  std::vector<sync_task> sync_tasks;
  

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
                            termination_reason(EXEC_UNSET),
                            scheduled_vertices(graph.owned_vertices().size()),
                            update_function(NULL),
                            max_iterations(0),
                            thread_color_barrier(ncpus){ 
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
      std::cout << "add task to " << task.vertex() << std::endl;
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

    std::vector<std::vector<std::pair<vertex_id_t, size_t> > > color_block_and_weight;
    // the list of vertices for each color
    color_block_and_weight.resize(graph.num_colors());
    
    foreach(vertex_id_t v, graph.owned_vertices()) {
      color_block_and_weight[graph.get_color(v)].push_back(
                                    std::make_pair(graph.globalvid_to_localvid(v), 
                                                  graph.globalvid_to_replicas(v).size()));
    }
    color_block.clear();
    color_block.resize(graph.num_colors());
    // optimize ordering. Sort in descending order
    // put all those which need a lot of communication in the front
    // to give communication the maximum amount if time possible.
    for (size_t i = 0; i < color_block_and_weight.size(); ++i) {
      std::sort(color_block_and_weight[i].rbegin(),
                color_block_and_weight[i].rend());
      // insert the sorted vertices into the final color_block

      std::transform(color_block_and_weight[i].begin(),
                     color_block_and_weight[i].end(), 
                     std::back_inserter(color_block[i]),
                     __gnu_cxx::select1st<std::pair<vertex_id_t, size_t> >());

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
  
  
  
  void check_global_termination(bool check_dynamic_schedule) {
    std::vector<termination_evaluation> termination_test;
    termination_test.resize(rmi.numprocs());
    
    if (check_dynamic_schedule) {
      termination_test[rmi.procid()].pending_tasks = num_pending_tasks.value;
    }

    if (task_budget > 0) {
      size_t numupdates = 0;
      for (size_t i = 0; i < update_counts.size(); ++i) numupdates += update_counts[i];
      termination_test[rmi.procid()].executed_tasks = numupdates;
    }
    
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
    if (rmi.procid() == 0) {
      termination_evaluation aggregate;
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
    rmi.broadcast(treason, rmi.procid() == 0);
    termination_reason = exec_status(treason);
  }
 
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_scope<Graph> scope;
                  

    // loop over iterations
    size_t iter = 0;
    bool usestatic = max_iterations > 0;
    while(1) {
      // if max_iterations is defined, quit
      if (usestatic && iter >= max_iterations) {
        termination_reason = EXEC_TASK_DEPLETION;
        break;
      }
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
            scope.commit_async_untracked();
            update_counts[threadid]++;
          }
        }      
        // wait for all threads to synchronize on this color.
        thread_color_barrier.wait();
        curidx.value = 0;
        // full barrier on the color
        // this will complete synchronization of all add tasks as well
        rmi.dc().full_barrier();
      }
      ++iter;
      if (threadid == 0) {
        check_global_termination(!usestatic);
      }
      // all threads must wait for 0
      thread_color_barrier.wait();
      if (termination_reason != EXEC_UNSET) break;
    }
  }
  
  /** Execute the engine */
  void start() {
    assert(update_function != NULL);
    // generate colors then
    // wait for everyone to enter start    
    generate_color_blocks();
    termination_reason = EXEC_UNSET;
    force_stop = false;
    numsyncs.value = 0;
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
                            &distributed_chromatic_engine<Graph>::start_thread,
                            this, i), aff);
    }
    
    thrgrp.join();              
    rmi.barrier();
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
    any uf;
    if (opts.get_any_option("update_function", uf)) {
      update_function = uf.as<update_function_type>();
    }
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

};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP
