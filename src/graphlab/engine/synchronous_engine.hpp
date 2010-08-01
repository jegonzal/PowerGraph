
#ifndef GRAPHLAB_SYNCHRONOUS_ENGINE_HPP
#define GRAPHLAB_SYNCHRONOUS_ENGINE_HPP

#include <cassert>
#include <cstdio>


#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/scope/synchronous_scope_factory.hpp>
#include <graphlab/schedulers/support/binary_scheduler_callback.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>

#include <graphlab/shared_data/ishared_data.hpp>
#include <graphlab/shared_data/ishared_data_manager.hpp>
#include <graphlab/monitoring/imonitor.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {
 
  /**
   * This class defines a synchronous engine. Writes are visible
   * only on next iteration.
   **/  
 
  template<typename Graph>
  class synchronous_engine : 
    public iengine<Graph> {

  public:   
    typedef iengine<Graph> base;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::ischeduler_type ischeduler_type;
    typedef typename base::imonitor_type imonitor_type;
    typedef typename base::termination_function_type termination_function_type;
    typedef typename base::iscope_type iscope_type;
    typedef typename base::ishared_data_type ishared_data_type;
    typedef typename base::ishared_data_manager_type ishared_data_manager_type;

    typedef synchronous_scope_factory<Graph> synch_scope_factory_type;
    typedef synchronous_scope<Graph> synch_scope_type;
    typedef binary_scheduler_callback<Graph> binary_callback_type;

  private:
    

    /** The number of cpus to use */
    size_t ncpus;
    std::vector<size_t> task_counts;
    
    Graph& src;
    Graph dest;
    update_function_type updatefunc;
    
    synch_scope_factory_type scope_manager;
    binary_callback_type callback;
    
    timer _timer;
    size_t taskbudget;
    atomic<size_t> taskcount;
    double timeout;
    bool aborted;
    bool cloneedges;
    barrier iterationbarrier;
    size_t niterations;
    
    imonitor_type* listener;
    ishared_data_manager_type* data_manager;
    
    size_t fixediterations; 
    std::vector<termination_function_type> term_functions;
    float lasttermcheck;

    exec_status termination_reason;

  public:
    /** Initialize the multi threaded engine */
    synchronous_engine(Graph& g, size_t num_cpus = thread::cpu_count()) :
      ncpus(num_cpus),
      task_counts(num_cpus,0),
      src(g),      
      scope_manager(src, dest, num_cpus), 
      taskcount(0), 
      iterationbarrier(num_cpus),
      listener(NULL), 
      data_manager(NULL) {
      timeout = 0;
      taskbudget = 0;
      aborted = false;
      updatefunc = NULL;
      assert(num_cpus >= 1);
      
      // Todo: fix this?
      // dest.copy_from(src);
      dest = src;


      niterations = 0;
      cloneedges = true;
    }

    size_t get_ncpus() const { return ncpus; } 

    /** register the listener */
    void register_monitor(imonitor_type* _listener) {
      if(_listener == NULL) return;
      this->listener = _listener;
      listener->init(this);
    }

    void set_default_scope(scope_range::scope_range_enum default_scope_range) {
      scope_manager.set_default_scope(default_scope_range);
    }


    /**
     * Timeout. Default - no timeout. 
     */
    void set_timeout(size_t timeout_secs) {
      timeout = timeout_secs;
    }
    
    /**
     * Task budget - max number of tasks to allow
     */
    virtual void set_task_budget(size_t max_tasks) {
      taskbudget = max_tasks;
    }

    void set_clone_edges(bool clone_edges) {
      if (clone_edges == false) {
        // TODO: support disabling edge cloning
        logger(LOG_WARNING, "Edge cloning cannot be disabled.(TODO)");
      }
      cloneedges=clone_edges;
    }
    
    bool check_all_terminators() {
      for (size_t i = 0;i < term_functions.size();++i) {
        if (term_functions[i](data_manager)) return true;
      }
      return false;
    }
 
    
    /** The internal worker thread class */
    class synchronous_worker : public runnable {      
      synchronous_engine* engine;
      size_t workerid;
    public:
      synchronous_worker() : engine(NULL), workerid(0) {  }
      void init(synchronous_engine* _engine, size_t _workerid) {
        engine = _engine;
        workerid = _workerid;
      } // End of init      
      void run() {
        assert(engine != NULL);
        logger(LOG_INFO, "Worker %d initialized\n", workerid);        
        logger(LOG_INFO, "Worker %d started\n", workerid);        
        /* Start consuming tasks */
        while(engine->iteration(workerid)) { }
        // Record worker death with the listener
        if (engine->listener != NULL)
          engine->listener->
            engine_worker_dies(workerid, 
                               engine->task_counts[workerid]);
        logger(LOG_INFO, "Worker %d died\n", workerid);
      }      
    }; // end of task worker

 
    bool iteration(int cpuid) {
      if (cpuid == 0) {
        ++niterations;
      }
      iterationbarrier.wait();
      for(vertex_id_t vertex = cpuid; vertex < src.num_vertices(); 
          vertex += ncpus) {
        task_counts[cpuid] = task_counts[cpuid]++;
        // build a scope
        iscope_type* scope = scope_manager.get_scope(cpuid, vertex);
        
        assert(scope != NULL);

        update_task_type task(vertex, updatefunc);
        
        if (listener != NULL) {
          listener->scheduler_task_scheduled(task, 0.0);
          listener->engine_task_execute_start(task, scope, cpuid);
        }
        // execute the task
        // task.execute(*scope, callback, data_manager);
        assert(task.function() != NULL);
        task.function()(*scope, callback, data_manager);

        if (listener != NULL)
          listener->engine_task_execute_finished(task, scope, cpuid);      
     
        scope->commit();

        scope_manager.release_scope(scope);

        // commit the callback ad update the state of the scheduler
      }
      iterationbarrier.wait();
      // abortion check
      if (cpuid == 0) {
        logger(LOG_INFO, "Iteration %d complete.", niterations);
      }
      if (cpuid == 0 && 
          ((fixediterations == 0 && callback.add_task_called == false) ||
           (fixediterations > 0 && fixediterations == niterations)  ||
           (taskbudget > 0 && niterations * src.num_vertices() > taskbudget) || 
           (timeout > 0 && _timer.current_time() > timeout)  || 
           check_all_terminators()) ) {
        aborted = true;
        logger(LOG_INFO, "Aborting Synchronous Engine");
      }

      iterationbarrier.wait();
      if (cpuid == 0) {
        scope_manager.swap_graphs();
        callback.reset();
      }
      
      return !aborted;
    }
   
    void set_update_function(update_function_type u) {
      updatefunc = u;
    }
   
    /** get a reference to the scheduler */
    ischeduler_type& get_scheduler() { 
      logger(LOG_FATAL, "get_scheduler() not supported for synchronous engine");
      exit(0);
      //      return *reinterpret_cast<ischeduler_type*>(NULL);
    }

    void set_shared_data_manager(ishared_data_manager_type* manager) {
      data_manager = manager;
      if(data_manager != NULL) {
        data_manager->set_scope_factory(&scope_manager);
      }
    }


    void add_terminator(termination_function_type term) {
      term_functions.push_back(term);
    }
    void clear_terminators() {
      term_functions.clear();
    }

    void start() {
      start_with_iteration_limit(0);
    }

    void stop() {
      assert(false); // Unsupported
    }

    /** Execute the preloaded tasks on the graph */
    void start_with_iteration_limit(size_t numiterations = 0) {
      termination_reason = EXEC_TASK_DEPLETION;
      assert(updatefunc != NULL);
      //! Finalize the graph (this could take a while so you should do
      //! it before calling start for timing purposes)
      src.finalize();

      // Ensure that the data manager has the correct scope_factory
      if(data_manager != NULL) {
        data_manager->set_scope_factory(&scope_manager);
      }


      dest = src;


      niterations = 0;
      fixediterations = 0;
      if (numiterations>0) fixediterations= numiterations;
      /* Timing */
      _timer.start();
      lasttermcheck = lowres_time_seconds();
      callback.reset();
      /* Enable scheduler to clean up in restarts */
      
      /* Initialize a pool of threads */
      std::vector<synchronous_worker> workers(ncpus);
      thread_group threads;
      for(size_t i = 0; i < ncpus; ++i) {
        // Initialize the worker
        workers[i].init(this, i);
        // Start the worker thread using the thread group with cpu
        // affinity attached (CPU affinity currently only supported in
        // linux) since Mac affinity is set through the NX frameworks
        threads.launch(&(workers[i]));
        if (listener != NULL)
          listener->engine_worker_starts(i);
      }
      /* Wait for all threads to return */
      logger(LOG_INFO, "Wait until finished...");
      threads.join();
      //double running_time = _timer.current_time();
      // ok. Now we need to know which is the right "final" graph
      // The src graph in the scope manager is the right graph
      if (scope_manager.get_src_graph() != &src) {
        src = *scope_manager.get_src_graph();
      }
      // now I need to clone the vertex data over (if necessary)
      if (scope_manager.get_vertex_data_graph() != &src) {             
        Graph* vdatagraph = scope_manager.get_vertex_data_graph();
        for(size_t i = 0; i < src.num_vertices(); ++i) {
          src.vertex_data(i) = vdatagraph->vertex_data(i);
        }
      }
      /**
       * Log task counts. It is useful to see worke-specific task
       * counts to see if work was distributed evenly
       */
      size_t total_counts = 0;
      for(size_t wid = 0; wid < task_counts.size(); ++wid) {
        total_counts += task_counts[wid];
        logger(LOG_INFO,
               "Worker %d finished: task count = %d",
               wid, task_counts[wid]);
      } // end of loop over task_counts
      logger(LOG_INFO, "=== Total task count: %d", total_counts);
      
      
      
      if (taskbudget > 0 && taskcount.inc() > taskbudget) {
        termination_reason = EXEC_TASK_BUDGET_EXCEEDED;
      } else {
        termination_reason = EXEC_TIMEOUT;
      }
   
      
    } // end of start

    exec_status last_exec_status() const { 
      return termination_reason; 
    }

    /**
     * Get the total number of updates executed by this engine.
     */
    size_t last_update_count() const { 
      size_t total_updates = 0;
      foreach(size_t tc, task_counts) {
        total_updates += tc;
      }
      return total_updates;
    } // end of get profiling info

  }; // end of synchronous_engine

  
      
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
