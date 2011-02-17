#ifndef DISTRIBUTED_CHROMATIC_ENGINE_HPP
#define DISTRIBUTED_CHROMATIC_ENGINE_HPP

#include <functional>

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
#include <graphlab/scheduler/support/redirect_scheduler_callback.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>
#include <graphlab/distributed2/graph/dgraph_scope.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {

#include <graphlab/macros_def.hpp>


/**
All processes must receive the same options. 
options are NOT synchronized. (i.e. if timeout is set, 
all processes must and should set the same timeout value, etc)
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

  typedef redirect_scheduler_callback<Graph, 
                                      distributed_chromatic_engine<Graph> > callback_type;
 private:

 private:
  // the local rmi instance
  dc_dist_object<distributed_chromatic_engine<Graph, Scheduler> > rmi;
  
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



  /** The timeout time in millis */
  size_t timeout_millis;
  timer ti;
  
  /** The total number of tasks that should be executed */
  size_t task_budget;
  
  /** If dynamic scheduling is used, the number of scheduled tasks */
  atomic<size_t> num_pending_tasks;
  
  /** Boolean that determins whether the engine is active */
  bool active; 

  /** The cause of the last termination condition */
  exec_status termination_reason;

  scope_range::scope_range_enum default_scope_range;
 
  std::vector<std::vector<vertex_id_t> > color_blocks;
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
  void distributed_chromatic_engine(distributed_control &dc,
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
                            start_time_millis(lowres_time_millis()),
                            timeout_millis(0),
                            last_check_millis(0),
                            task_budget(0),
                            active(false),
                            termination_reason(EXEC_UNSET),
                            thread_color_barrier(ncpus),
                            update_function(NULL),
                            max_iterations(0) { 
    scheduled_vertices.resize(graph.owned_vertices().size());
    scheduled_vertices.clear();
  }
  
  void ~distributed_chromatic_engine() {
  }
  
  
  //! Get the number of cpus
  size_t get_ncpus() const { return ncpus; }

  //! set sched yield
  void set_sched_yield(bool value) {
    use_sched_yield = value;
  }

  void set_cpu_affinities(bool value) {
    use_cpu_affinity = value;
  }


  /**
   * Set the default scope range.  The scope ranges are defined in
   * iscope.hpp
   */
  void set_default_scope(scope_range::scope_range_enum default_scope_range_) {
    default_scope_range = default_scope_range_;
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

  /** This function provides an approximation to the last update count.
   * This is a faster version of last_update_count and may be off
   by at most (APX_INTERVAL+1)*Nthreads */
  inline size_t approximate_last_update_count() const {
    return apx_update_counts.value;
  }


  /**
   * Add a terminator to the engine.
   */
  void add_terminator(termination_function_type term) {
    term_functions.push_back(term);
  }


  /**
   * Clear all terminators from the engine
   */
  void clear_terminators() {
    term_functions.clear();
  }



  /**
   * Timeout. Default - no timeout. 
   */
  void set_timeout(size_t timeout_seconds = 0) {
    timeout_millis = timeout_seconds * 1000;
  }
  
  /**
   * Task budget - max number of tasks to allow
   */
  virtual void set_task_budget(size_t max_tasks) {
    task_budget = max_tasks;
  }
  


  /**
   * \brief Adds an update task with a particular priority.
   * This function is forwarded to the scheduler.
   */
  virtual void add_task(update_task_type task, double priority) {
    if (update_function == NULL) assert update_function == task.function();
    if (graph.is_owned(task.vertex())) {
      num_pending_tasks.inc(!
              scheduled_vertices.set_bit(globalvid_to_localvid(task.vertex())) 
                          );
    }
    else {
      rmi.remote_call(globalvid_to_owner(task.vertex()),
                      &distributed_chromatic_engine<Graph>::add_task,
                      reinterpret_cast<size_t>(task),
                      priority);
    }
  }

  /**
   * \brief Creates a collection of tasks on all the vertices in
   * 'vertices', and all with the same update function and priority
   * This function is forwarded to the scheduler.
   */
  virtual void add_tasks(const std::vector<vertex_id_t>& vertices,
                         update_function_type func, double priority) {
    for (size_t i = 0;i < vertices.size(); +i) {
      add_task(update_task_type(vertices[i], func), priority);
    }
  }

  /**
   * \brief Creates a collection of tasks on all the vertices in the graph,
   * with the same update function and priority
   * This function is forwarded to the scheduler.
   */
  virtual void add_task_to_all(update_function_type func,
                               double priority) {
      scheduled_vertices.fill();
      num_pending_tasks.value = scheduled_vertices.size();
      for (size_t i = 0;i < rmi.numprocs(); ++i) {
        if (i != rmi.procid()) {
          rmi.remote_call(i,
                          &distributed_chromatic_engine<Graph>::add_task_to_all,
                          func,
                          priority);
        }
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
    use_fake_shared_data = true;
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
    var2synctask[&shared] = sync_tasks.size() - 1;
  }

  
  void generate_color_blocks() {
    // construct for each color, the set of vertices as well as the 
    // number of replicas for that vertex.
    // the number of replicas - 1 is the amount of communication
    // we have to perform to synchronize modifications to that vertex

    std::vector<std::pair<vertex_id_t, size_t> > color_block_and_weight;
    // the list of vertices for each color
    color_blocks_and_weight.resize(graph.num_colors());
    size_t numghostv = 0;
    foreach(vertex_id_t v, graph.owned_vertices()) {
      color_block_and_weight[graph.get_color(v)].push_back(
                                    std::make_pair(v, 
                                    g.globalvid_to_replicas(v)));
    }
    color_block.clear();
    color_block.resize(graph.num_color());
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
                     __gnu_cxx::select1st<std::pair<vertex_id_t, size_t> >);
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
    termination_evaluation(): pending_tasks(0),
                              executed_tasks(0),
                              terminator(false),
                              timeout(false) { }
  };
  void check_global_termination(bool check_dynamic_schedule) {
    std::vector<termination_evaluation> termination_test;
    termination_test.resize(rmi.numprocs());
    
    if (check_dynamic_schedule) {
      termination_test[rmi.procid()].pending_tasks = pending_tasks.value;
    }

    if (task_budget > 0) {
      size_t numupdates = 0;
      for (size_t i = 0; i < update_counts.size(); ++i) numupdates += update_counts[i];
      termination_test[rmi.procid()].executed_tasks = numupdates;
    }
    
    if (timeout_millis > 0 & ti.current_time_millis() > timeout_millis) {
      termination_test[rmi.procid()].timeout = true;
    }
    
    for (size_t i = rmi.procid(); i < term_functions.size(); i += rmi.numprocs()) {
      if (term_functions[i](NULL)) {
        termination_test[rmi.procid()].terminator = true;
        break;
      }
    }
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
      }
      
      if (check_dynamic_schedule && aggregate.pending_tasks == 0) {
        termination_reason = EXEC_TASK_DEPLETION;
      }
      else if (task_budget > 0 && aggregate.executed_tasks > task_budget) {
        termination_reason = EXEC_TASK_BUDGET_EXCEEDED;
      }
      else if (timeout_millis > 0 && aggregate.timeout) {
        termination_reason = EXEC_TIMEOUT;
      }
      else if (aggregate.terminator) {
        termination_reason = EXEC_TERM_FUNCTION;
      }
    }
    size_t treason = termination_reason;
    rmi.broadcast(treason, 0);
    termination_reason = treason;
  }
 
  void start_thread(size_t threadid) {
    // create the scope
    dgraph_scope<Graph> scope;
          typename Scheduler::callback_type& scallback =
                  scheduler->get_callback(cpuid);
                  

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
          if (usestatic || scheduled_vertices.clear_bit(i)) {
            // if index out of scope, we are done with this color. break
            if (i >= color_block[c].size()) break;
            // otherwise. run the vertex
            // create the scope
            scope.init(&graph, color_block[c][i]);
            // run the update function
            task.function()(scope, callback, NULL);
            graph.push_owned_scope_to_replicas(color_block[c][i],
                                               true,  // modified only
                                               true,  // clear modified
                                               true,  // async
                                               true);  // untracked
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
        test_global_termination(!usestatic);
      }
      // all threads must wait for 0
      thread_color_barrier.wait();
      if (termination_reason != EXEC_UNSET) break;
    }
  }
  
  /** Execute the engine */
  void start() {
    // generate colors then
    // wait for everyone to enter start    
    generate_color_blocks();
    termination_reason = EXEC_UNSET;
    rmi.barrier();
    // reset indices
    curidx.value = 0;
    ti.start();
    // spawn threads
    thread_group thrgrp; 
    for (size_t i = 0;i < ncpus; ++i) {
      size_t aff = use_cpu_affinity ? i : -1;
      launch_in_new_thread(thrgrp, 
                         boost::bind(
                            distributed_chromatic_engine<Graph>::start_thread,
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
  virtual void set_scheduler_options(const scheduler_options& opts) {
    opts.get_int_option("max_iterations", max_iterations);
    any uf;
    if (opts.get_any_option("update_function", uf)) {
      update_function = uf.as<update_function_type>();
    }
  }
  
  static void print_options_help(std::ostream &out) {
    out << "max_iterations = [integer, default = 0]\n";
    out << "update_function = [update_function_type,"
      "default = set on add_task]\n";
  };

  

};

} // namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP
