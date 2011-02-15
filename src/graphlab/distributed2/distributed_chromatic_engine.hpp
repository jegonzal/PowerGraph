#ifndef DISTRIBUTED_CHROMATIC_ENGINE_HPP
#define DISTRIBUTED_CHROMATIC_ENGINE_HPP

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/util/mutable_queue.hpp>

#include <graphlab/engine/iengine.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/metrics/metrics.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>


namespace graphlab {

#include <graphlab/macros_def.hpp>

template <typename Graph, typename Scheduler>
class distributed_chromatic_engine:public iengine<Graph> {
 public:
  typedef iengine<Graph> iengine_base;
  typedef typename iengine_base::update_task_type update_task_type;
  typedef typename iengine_base::update_function_type update_function_type;
  typedef typename iengine_base::ischeduler_type ischeduler_type;
  typedef typename iengine_base::termination_function_type termination_function_type;
  typedef typename iengine_base::iscope_type iscope_type;
  typedef typename iengine_base::ishared_data_type ishared_data_type;
  typedef typename iengine_base::ishared_data_manager_type ishared_data_manager_type;
  
  typedef typename iengine_base::sync_function_type sync_function_type;
  typedef typename iengine_base::merge_function_type merge_function_type;

  using base::apply_scheduler_options;
  using base::release_scheduler_and_scope_manager;
  using base::get_scheduler;
  using base::get_scope_manager;

 private:

 private:
  // the local rmi instance
  dc_dist_object<distributed_chromatic_engine<Graph, Scheduler> > rmi;
  
  // the graph we are processing
  Graph &graph;
  
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
  /** track an approximation to the number of updates. This 
      is only updated every (APX_INTERVAL+1) updates per thread.
      i.e. this can be off by at most (APX_INTERVAL+1)*Nthreads */
  atomic<size_t> apx_update_counts;
  
  /// frequency where apx_update_counts is updated. Must be a power of 2 - 1`
  static const size_t APX_INTERVAL = 127;

  atomic<size_t> numsyncs;


 /** The time in millis that the engine was started */
  size_t start_time_millis;

  /** The timeout time in millis */
  size_t timeout_millis;

  /** The last time a check was run */
  size_t last_check_millis;
  
  /** The total number of tasks that should be executed */
  size_t task_budget;
  
  /** Boolean that determins whether the engine is active */
  bool active; 

  /** The cause of the last termination condition */
  exec_status termination_reason;

  scope_range::scope_range_enum default_scope_range;
 
    
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
                            termination_reason(EXEC_UNSET) { }
  
  void ~distributed_chromatic_engine() {
  }
  
  
  //! Get the number of cpus
  size_t get_ncpus() const { return ncpus; }

  //! set sched yield
  void enable_sched_yield(bool value) {
    use_sched_yield = value;
  }

  void enable_cpu_affinities(bool value) {
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

  /** Execute the engine */
  void start() {
    \\ TODO
  }
  
  /**
   * Performs a sync immediately. This function requires that the shared
   * variable already be registered with the engine.
   * and that the engine is not currently running
   */
  void sync_now(glshared_base& shared) {
    \\ TODO
  }
  
  

};

} // namespace graphlab

#endif // DISTRIBUTED_CHROMATIC_ENGINE_HPP
