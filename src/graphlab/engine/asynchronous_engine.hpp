#ifndef GRAPHLAB_ASYNCHRONOUS_ENGINE_HPP
#define GRAPHLAB_ASYNCHRONOUS_ENGINE_HPP

#include <cmath>
#include <cassert>
#include <algorithm>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/random.hpp>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/schedulers/support/scheduler_option_cache.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * This class defines a basic asynchronous engine
   */
  template<typename Graph, typename Scheduler, typename ScopeFactory>
  class asynchronous_engine : 
    public iengine<Graph> {

  public: // Type declerations
    enum execution_type {THREADED, SIMULATED};
    typedef iengine<Graph> base;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::ischeduler_type ischeduler_type;
    typedef typename base::imonitor_type imonitor_type;
    typedef typename base::termination_function_type termination_function_type;
    typedef typename base::iscope_type iscope_type;
    typedef typename base::ishared_data_type ishared_data_type;
    typedef typename base::ishared_data_manager_type ishared_data_manager_type;


    /** The internal worker thread class used for the threaded engine */
    class engine_thread : public runnable {      
      asynchronous_engine* engine;
      size_t workerid;
    public:
      engine_thread() : engine(NULL), workerid(0) {  }
      void init(asynchronous_engine* _engine, size_t _workerid) {
        engine = _engine;
        workerid = _workerid;
      } // End of init      
      void run() {
        assert(engine != NULL);
        logger(LOG_INFO, "Worker %d started.\n", workerid);        
        /* Start consuming tasks while the engine is active*/
        while(engine->active) {
          bool executed_task = engine->run_once(workerid);         
          // If this was nothing to execute then fail
          if (!executed_task) break;
        }   
        // // Do any remaining syncs if any
        if(engine->shared_data != NULL)
          engine->shared_data->signal_all();
        logger(LOG_INFO, "Worker %d finished.\n", workerid);
      }      
    }; // end of task worker

    

  private: // data members

    /** The graph that this engine is executing */
    Graph& graph;
    
    /** Number of cpus to use */
    size_t ncpus; 

    /** Specify the execution type */
    execution_type exec_type;

    /** Use processor affinities */
    bool use_cpu_affinity;

    /** Use schedule yielding when waiting on the scheduler*/
    bool use_sched_yield;
    
    /** Responsible for managing the update of scopes */
    ScopeFactory *scope_manager;

    /** Responsible for maintaining the schedule over tasks */
    Scheduler *scheduler;

    /** Track the number of updates */
    std::vector<size_t> update_counts;

    /** The monitor which tracks and records engine events */
    imonitor_type* monitor;

    /** The share data manager */
    ishared_data_manager_type* shared_data;

    /** The time in millis that the engine was started */
    size_t start_time_millis;

    /** The timeout time in millis */
    size_t timeout_millis;

    /** The last time a check was run */
    size_t last_check_millis;
    
    /** The total number of tasks that should be executed */
    size_t task_budget;

    /** The termination functions */
    std::vector<termination_function_type> term_functions;

    /** Boolean that determins whether the engine is active */
    bool active;

    /** The cause of the last termination condition */
    exec_status termination_reason;

    scope_range::scope_range_enum default_scope_range;

    scheduler_option_cache sched_opt_cache;
    
    void construct_scheduler() {
      if (scheduler == NULL) {
        scheduler = new Scheduler(this, graph, std::max(ncpus, size_t(1)));
        sched_opt_cache.apply_options(*scheduler);
      }
    }
    void construct_scope_manager() {
      if (scope_manager == NULL) {
        scope_manager = new ScopeFactory(graph, std::max(ncpus, size_t(1)));
      }
    }
    
  public:

    /**
     * Create an asynchronous engine.
     *
     * The execution type specifies whether the engine should simulate
     * threads or use actual threads. 
     */
    asynchronous_engine(Graph& graph,
                        size_t ncpus = 1,
                        execution_type exec_type = THREADED) :
      graph(graph),
      ncpus( std::max(ncpus, size_t(1)) ),
      exec_type(exec_type),
      use_cpu_affinity(false),
      use_sched_yield(true),
      scope_manager(NULL),
      scheduler(NULL),
      update_counts(std::max(ncpus, size_t(1)), 0),
      monitor(NULL),
      shared_data(NULL),
      start_time_millis(lowres_time_millis()),
      timeout_millis(0),
      last_check_millis(0),
      task_budget(0),
      active(false) { }

    ~asynchronous_engine() {
      if (scope_manager != NULL) delete scope_manager;
      if (scheduler != NULL) delete scheduler;
    }

    //! Get the number of cpus
    size_t get_ncpus() const { return ncpus; }

    void set_sched_option(std::stringstream &strm) {
      sched_opt_cache.parse_options(strm);
    }

    void set_sched_option(scheduler_options::options_enum opt, void* value) {
      sched_opt_cache.set_option(opt, value);
    }

    void set_sched_option(const scheduler_option_cache &cache) {
      sched_opt_cache = cache;
    }

    //! set sched yeild
    void enable_sched_yield(bool value) {
      use_sched_yield = value;
    }

    void enable_cpu_affinities(bool value) {
      use_cpu_affinity = value;
    }

    //! Get the scheduler associated with this engine
    ischeduler_type& get_scheduler() {
      construct_scheduler();
      return *scheduler;
    }
    
    //! Set the shared data manager for this engine
    void set_shared_data_manager(ishared_data_manager_type* _shared_data) {
      shared_data = _shared_data;
      // if the data manager is not null then the scope factory is
      // passed back
      if(shared_data != NULL) {
        shared_data->set_scope_factory(scope_manager);
      }
    } // end of set shared data manager


    /**
     * Set the default scope range.  The scope ranges are defined in
     * iscope.hpp
     */
    void set_default_scope(scope_range::scope_range_enum default_scope_range_) {
      default_scope_range = default_scope_range_;
    }


    /** Execute the engine */
    void start() {
      construct_scheduler();
      construct_scope_manager();
      sched_opt_cache.apply_options(*scheduler);
      scope_manager->set_default_scope(default_scope_range);
      /**
       * Prepare data structures for execution:
       * 1) finalize the graph.
       * 2) Reset engine fields
       */
      // Prepare the graph
      graph.finalize();      
      // Clear the update counts
      std::fill(update_counts.begin(), update_counts.end(), 0);
      // Reset timers
      start_time_millis = lowres_time_millis();
      last_check_millis = 0;
      // Reset active flag
      active = true;  
      // Reset the last exec status 
      termination_reason = EXEC_TASK_DEPLETION;
      // Ensure that the data manager has the correct scope_factory
      if(shared_data != NULL) shared_data->set_scope_factory(scope_manager);
      // Start any scheduler threads (if necessary)
      scheduler->start();
      
      /**
       * Depending on the execution type call the correct internal
       * start
       */
      if(exec_type == THREADED) run_threaded();
      else run_simulated();
      
      /**
       * Run any necessary cleanup
       */
      delete scheduler;
      // TODO: delete scope manager here
      //delete scope_manager;
      scheduler = NULL;
      //scope_manager = NULL;
    } // End of start


    /**
     * Stop the engine
     */
    void stop() {
      termination_reason = EXEC_FORCED_ABORT;
      active = false;
    }
    

    /**
     * Return the reason why the engine last terminated
     */
    exec_status last_exec_status() const {
      return termination_reason;
    }

    
    /**
     * This function computes the last update count by adding all the
     * update counts of the individual threads.  This will lead to a
     * potential data race when called while the engine is executing.
     * However, this data race will only produce an overly
     * conservative estimate.
     */
    size_t last_update_count() const {
      size_t sum = 0;
      for(size_t i = 0; i < update_counts.size(); ++i)
        sum += update_counts[i];
      return sum;
    } // end of last_update_count



    /**
     * Register a monitor with this engine.  Currently this engine
     * only supports a single monitor. 
     */
    void register_monitor(imonitor<Graph>* _monitor = NULL) {
      monitor = _monitor;
      scheduler->register_monitor(monitor);
      if(monitor != NULL) monitor->init(this);
    } // end of register monitor


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
    


  protected: // internal functions

    /**
     * Execute the engine using actual threads
     */
    void run_threaded() {
      /* Initialize a pool of threads */
      std::vector<engine_thread> workers(ncpus);
      thread_group threads;
      for(size_t i = 0; i < ncpus; ++i) {
        // Initialize the worker
        workers[i].init(this, i);
        // Start the worker thread using the thread group with cpu
        // affinity attached (CPU affinity currently only supported in
        // linux) since Mac affinity is set through the NX frameworks
        if(use_cpu_affinity)  {
          threads.launch(&(workers[i]), i);
        } else {
          threads.launch(&(workers[i]));
        }
      }
      threads.join();
    } // end of run threaded


    /**
     * Simulate the use of actual threads.
     */
    void run_simulated() {
      use_sched_yield = false;
      // repeatedly invoke run once as a random thread
      while(active) {
        // Pick a random cpu to run as
        size_t cpuid = 0;
        if(ncpus > 1) {
          cpuid = random::rand_int(ncpus - 1);
        }
        // Execute the update as that cpu
        active = run_once(cpuid);
      }
      // Do any remaining syncs if any
      if(shared_data != NULL) shared_data->signal_all();
    } // end of run simulated

    
    
    /**
     * Check all the terminators to identify if any termination
     * condition has been met. If a termination condition is met this
     * function returns true otherwise it will return false.  
     */
    bool satisfies_termination_condition() {
      /**
       * Timeout termination condition.
       *
       * If a timeout was set and the current time is greater than the
       * timeout point return true
       */
      if(timeout_millis > 0 && 
         start_time_millis + timeout_millis < lowres_time_millis()) {
        termination_reason = EXEC_TIMEOUT;
        // Deactivate the engine
        return true;
      }

      /**
       * Task budget termination condition
       *
       * If the task budget is greater than zero and the last update
       * count exceeds the task budget then terminate
       */
      if(task_budget > 0 && last_update_count() > task_budget) {
        termination_reason = EXEC_TASK_BUDGET_EXCEEDED;
        // Deactivate the engine
        return true;
      }

      /**
       * Check all the terminators.
       *
       * If a shared data manager is not provided then this will fail.
       */
      for (size_t i = 0; i < term_functions.size(); ++i) {
        if (term_functions[i](shared_data)) {
          termination_reason = EXEC_TERM_FUNCTION;
          return true;
        }
      }
      // No termination condition was met
      return false;
    } // End of satisfies termination condition




    bool run_once(size_t cpuid) {
      // Loop until we get a task for recieve a termination signal
      while(active) {       
        /**
         * Run any pending syncs and then test all termination
         * conditions.
         */
        if (last_check_millis < lowres_time_millis()) {
          last_check_millis = lowres_time_millis();
          // If a data manager is available try and run any pending
          // synchronizations.
          if(shared_data != NULL) shared_data->signal_all();
          // Check all termination conditions
          if(satisfies_termination_condition()) {
            active = false;
            return false;
          }
        }
        
        /**
         * Get and execute the next task from the scheduler.
         */
        update_task_type task;
        sched_status::status_enum stat = scheduler->get_next_task(cpuid, task);
        
        if (stat == sched_status::EMPTY) {
          // check the schedule terminator
          scheduler->get_terminator().begin_critical_section(cpuid);
          stat = scheduler->get_next_task(cpuid, task);
          if (stat == sched_status::NEWTASK) {
            scheduler->get_terminator().cancel_critical_section(cpuid);
          }
          else {
            if (scheduler->get_terminator().end_critical_section(cpuid)) {
              termination_reason = EXEC_TASK_DEPLETION;
              active = true;
              return false;
            }
            else {
              if(use_sched_yield) sched_yield();
              if (exec_type != SIMULATED) continue;
              else return true;
            }
          }
        }
        
        if (stat == sched_status::NEWTASK) {
          // If the status is new task than we must execute the task
          const vertex_id_t vertex = task.vertex();
          assert(vertex < graph.num_vertices());
          assert(task.function() != NULL);
          // Lock the vertex to ensure that no other processor tries
          // to take it build a scope
          iscope_type* scope = scope_manager->get_scope(cpuid, vertex);
          assert(scope != NULL);                    
          // get the callback for this cpu
          typename Scheduler::callback_type& scallback =
                  scheduler->get_callback(cpuid);
          // execute the task
          task.function()(*scope, scallback, shared_data);
          // Commit any changes to the scope
          scope->commit();
          // Release the scope
          scope_manager->release_scope(scope);
          // Mark the task as completed in the scheduler
          scheduler->completed_task(cpuid, task);
          // record the successful execution of the task
          update_counts[cpuid]++;
          return true;
          break;
        } // end of switch
      } // end of while(true)
      // If this point is reached 
      return false;
    } // End of run once
    
   
  }; // end of asynchronous engine


}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
