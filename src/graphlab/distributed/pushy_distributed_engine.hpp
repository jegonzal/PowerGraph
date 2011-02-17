/**
 * Distributed engine that pushes updated data across the computing
 * network.
 **/
#ifndef GRAPHLAB_PUSHY_DISTRIBUTED_ENGINE_HPP
#define GRAPHLAB_PUSHY_DISTRIBUTED_ENGINE_HPP

#include <cassert>
#include <cstdio>



#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/distributed/metrics/distributed_metrics.hpp> 

#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  
  
  
  /**
   * This class defines a simple multi threaded engine. 
   */
  template<typename DistGraph, typename Scheduler, typename ScopeFactory>
  class pushy_distributed_engine : 
    public iengine<DistGraph> {

  public:   
    typedef iengine<DistGraph> base;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::ischeduler_type ischeduler_type;
    typedef typename base::imonitor_type imonitor_type;
    typedef typename base::termination_function_type termination_function_type;
    typedef typename base::iscope_type iscope_type;
    typedef typename base::ishared_data_manager_type ishared_data_manager_type;

  private:
    
    /** The internal worker thread class */
    class task_worker : public runnable {      
      pushy_distributed_engine* engine;
      size_t workerid;
    public:
      task_worker() : engine(NULL), workerid(0) {  }
      void init(pushy_distributed_engine* _engine, size_t _workerid) {
        engine = _engine;
        workerid = _workerid;
      } // End of init      
      void run() {
        assert(engine != NULL);
        logger(LOG_INFO, "Worker %d started\n", workerid);        
        /* Start consuming tasks */
        while(true) {
          bool executed_task = engine->run_next_task(workerid);
         
          // If this was nothing to execute then fail
          if (!executed_task)  break;
        }
        // // Do any remaining syncs if any
        if(engine->data_manager != NULL)
          engine->data_manager->signal_all();
        
        // Record worker death with the listener
        if (engine->listener != NULL)
          engine->listener->
            engine_worker_dies(workerid, 
                               engine->task_counts[workerid]);
        logger(LOG_INFO, "Worker %d died\n", workerid);
      }      
    }; // end of task worker
    
    /** Responsible for managing the update of scopes */
    ScopeFactory scope_manager;
    distributed_control &dc;

    /** Responsible for maintaining the schedule over tasks */
    Scheduler scheduler;

    /** The number of cpus to use */
    size_t ncpus;

    /* Worker bookkeeping */
    std::vector<size_t> task_counts;
    std::vector<size_t> worker_works;

    
    /* Listener */
    imonitor_type* listener;

    ishared_data_manager_type* data_manager;

    DistGraph& _distgraph;
    
    timer * _timer;
    
    size_t taskbudget;
    atomic<size_t> taskcount;

    double timeout;
    bool aborted;
    bool terminatoraborted;

    exec_status termination_cause;
	
    
    std::vector<termination_function_type> term_functions;
    float lasttermcheck;
  public:
    /** Initialize the multi threaded engine */
    pushy_distributed_engine(distributed_control &_dc, DistGraph& g, 
                             size_t num_cpus = thread::cpu_count()) :
      scope_manager(g, num_cpus),
      dc(_dc),
      scheduler(this, _dc, g, num_cpus),
      ncpus(num_cpus),
      task_counts(num_cpus,0),
      worker_works(num_cpus,0),
      listener(NULL),
      data_manager(NULL),
      _distgraph(g), taskcount(0) {
      scheduler.register_monitor(NULL);
      _timer = NULL;
    
      timeout = 0;
      taskbudget = 0;
      aborted = false;
      terminatoraborted = false;
      termination_cause = EXEC_TASK_DEPLETION;

      assert(num_cpus > 0);
    }

    size_t get_ncpus() const { return ncpus; } 

    void set_default_scope(scope_range::scope_range_enum default_scope_range) {
      scope_manager.set_default_scope(default_scope_range);
    }

    exec_status last_exec_status() const {
      return termination_cause;
    }

    /** register the listener */
    void register_monitor(imonitor_type* _listener) {
      if(_listener == NULL) return;
      this->listener = _listener;
      scheduler.register_monitor(listener);
      listener->init(this);
    }
    
    void set_sched_yield(bool value) {
      logger(LOG_INFO, "distributed engine does not support set_sched_yield()");    
    }

    void set_cpu_affinities(bool value) {
      logger(LOG_INFO, "distributed engine does not support set_cpu_affinities()");    
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

    bool check_all_terminators() {
      if(data_manager != NULL) {        
        for (size_t i = 0;i < term_functions.size();++i) {
          if (term_functions[i](data_manager)) return true;
        }
        return false;
      } else if (term_functions.size() > 0) {

        logger(LOG_WARNING,
               "Assessing termination without a data_manager!");
      }
      return false;

    }
    
    /**
     * Pushes all changes to neighboring partitions.
     * Unfortunately currently we do not track changes,
     * so all changes must be sent...
     */
    void push_changes(iscope_type *scope) { 
      // Replicate the vertex 
      _distgraph.update_vertex(scope->vertex());
     	
      // Replicate remote edges. Distribute graph
      // takes care that edge is update to proper owner.
      if (!_distgraph.has_constant_edges()) {
        foreach(edge_id_t e, scope->out_edge_ids()) {
          _distgraph.update_edge(e);
        }
      }
    }
        
    /**
     * Called by workers. out_of_work flag is passed as a pointer so
     * task's employment status can be updated immediatelly.  This
     * function must be thread safe.
     */
    bool run_next_task(int cpuid) {
      update_task_type task;

      while(true) {
        
        /* Check for task budget and timeout (only cpuid 0 checks) */
        if (timeout > 0 && cpuid == 0) {
          if (_timer->current_time() > timeout) {
            aborted = true;
          }
        }
        if (lowres_time_seconds() - lasttermcheck >= 0.1) {
          if (check_all_terminators()) terminatoraborted = true;
          if(data_manager != NULL) data_manager->signal_all();
          lasttermcheck = lowres_time_seconds();
          
        }
        
        if (aborted || terminatoraborted) return false;
        
        // Get the next task along with the status of the scheduler
        sched_status::status_enum stat = scheduler.get_next_task(cpuid, task);
          
        if (stat == sched_status::WAITING) { // If the status is waiting sleep and
          usleep(10);          // then try again
        } else if (stat == sched_status::COMPLETE) {
          /* Scheduler says we terminate */
          return false;
        } else if (stat == sched_status::NEWTASK) {
          const vertex_id_t vertex = task.vertex();
          // There is a new task to schedule 
          // Ensure that a valid task function was passed
          assert(task.function() != NULL);
          // Lock the vertex to ensure that no other processor tries
          // to take it
          // build a scope
          iscope_type* scope = scope_manager.get_scope(cpuid, vertex);
          
          assert(scope != NULL);
          
          // Update task counts and "work". Work is indegree+outdegree          
          task_counts[cpuid] = task_counts[cpuid]++;
          worker_works[cpuid] += scope->in_edge_ids().size() +
            scope->out_edge_ids().size();
          
          // get the callback for this cpu
          typename Scheduler::callback_type& scallback =
            scheduler.get_callback(cpuid);
          
          scallback.enable_buffering();

          if (listener != NULL)
            listener->engine_task_execute_start(task, scope, cpuid);

          // execute the task
          // task.execute(*scope, scallback, data_manager);
          assert(task.function() != NULL);
          scheduler.started_task(cpuid, task);

          task.function()(*scope, scallback, data_manager);

          if (listener != NULL)
            listener->engine_task_execute_finished(task, scope, cpuid);      
    	  push_changes(scope);
    	  

    	  scallback.commit();

          scope->commit();

        
          // Let shared data manager do its work
          if (data_manager != NULL) 
            data_manager->progress(cpuid, scope);
        
          scope_manager.release_scope(scope);

          scheduler.completed_task(cpuid, task);
         
          if (taskbudget > 0 && (taskcount.inc() > taskbudget)) {
            aborted = true;
          }
          
          return true;
        } // end of else if (stat == newtask)
      } // end of while(true)


      
    } // end of run_next_task

   
    /** get a reference to the scheduler */
    ischeduler_type& get_scheduler() { return scheduler; }

  
    void set_shared_data_manager(ishared_data_manager_type* manager) {
      data_manager = manager;      
      manager->set_scope_factory(&scope_manager);
    }

    //     ishared_data_type* get_shared_data() {
    //       return data_manager;
    //     }

    void add_terminator(termination_function_type term) {
      term_functions.push_back(term);
    }

    void clear_terminators() {
      term_functions.clear();
    }

    /** Execute the preloaded tasks on the graph */
    void start() {
      //! Finalize the graph (this could take a while so you should do
      //! it before calling start for timing purposes)
      _distgraph.finalize();
      
      
      // First do a bandwith test. TODO: remove
      distributed_metrics::instance(&dc)->execute_bandwith_test();

      // Ensure that the data manager has the correct scope_factory
      /*TODO when sync is ready:
        if(data_manager != NULL) {
        data_manager->set_scope_factory(this);
        }*/


      
      aborted = false;
      terminatoraborted = false;
      /* Timing */
      if (_timer == NULL) {
        /* If graphlab called in a loop, timer is started only once */
        _timer = new timer();
        _timer->start();
      }
      lasttermcheck = lowres_time_seconds();
      /* Enable scheduler to clean up in restarts */
      scheduler.start();
      /* Initialize a pool of threads */
      std::vector<task_worker> workers(ncpus);
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
     
      double running_time = _timer->current_time();
      logstream(LOG_INFO) << "Running time: " << running_time << std::endl;
      /**
       * Log task counts. It is useful to see worke-specific task
       * counts to see if work was distributed evenly
       */
      size_t total_counts = 0;
      size_t total_work = 0;
      for(size_t wid = 0; wid < task_counts.size(); ++wid) {
        total_counts += task_counts[wid];
        total_work += worker_works[wid];
        logger(LOG_INFO,
               "Worker %d finished: task count = %d, work = %ld",
               wid, task_counts[wid], worker_works[wid]);
      } // end of loop over task_counts
      logger(LOG_INFO, "=== Total task count: %d,   work=%ld", total_counts, 
             total_work);

      
      /// ===== START STATS OUTPUT ===== ///

      /* Output running time into a file (used by benchmark app) */
      char statsfilename[255];
      sprintf(statsfilename, ".runstats_%d_%d.R", dc.procid(), dc.numprocs());
      FILE * F = fopen(statsfilename, "w");
      fprintf(F, "engine=\"pushy\"\nexecution_time=%lf\nncpus=%d\ntaskcount=%ld\nwork=%ld\nresidual=0\nmemory_writes_mb=0\nmemory_reads_mb=0\n", running_time, (int) ncpus,total_counts, total_work);
      fclose(F);
      /// ===== END STATS OUTPUT ===== ///
      distributed_metrics::instance(&dc)->set_value("execution_time", running_time); 
      distributed_metrics::instance(&dc)->set_value("taskcount", total_counts);
      distributed_metrics::instance(&dc)->set_value("ncpus", ncpus);
      dc.report_stats();
    
      if (!aborted) {
      
        termination_cause =  EXEC_TASK_DEPLETION;
        return;
      } else {
        if (taskbudget > 0 && taskcount.inc() > taskbudget) {
          termination_cause = EXEC_TASK_BUDGET_EXCEEDED;
          return;
        } else {
          termination_cause = EXEC_FORCED_ABORT;
          return;
        }
      }
      
    } // end of start


    size_t last_update_count() const {
      size_t total_updates = 0;
      foreach(size_t tc, task_counts) {
        total_updates += tc;
      }
      return total_updates;    
    }

    void stop() {
      assert(false);
      // This should kill the engine execution
    }


    /**
     * Get the profiling info associated with a particular string
     * (key)
     */
    int64_t get_profiling_info(const std::string &prof) {
      if (prof == "update_count") {
        size_t total_updates = 0;
        foreach(size_t tc, task_counts) {
          total_updates += tc;
        }
        return total_updates;
      }
      logger(LOG_WARNING,
             "multi thread engine does not have profiling"
             "information for \'%s\'",
             prof.c_str());
      return 0;
    } // end of get profiling info

    std::vector<std::string> supported_profiling_info() {
      std::vector<std::string> ret;
      ret.push_back("update_count");
      return ret;
    }

    /////////////////Stuff not implemented ////////////////////////
    void add_task(update_task_type task, double priority) {}

    /**
     * Creates a collection of tasks on all the vertices in
     * 'vertices', and all with the same update function and priority
     * This function is forwarded to the scheduler.
     */
    void add_tasks(const std::vector<vertex_id_t>& vertices,
                   update_function_type func, double priority) {}

    /**
     * Creates a collection of tasks on all the vertices in the graph,
     * with the same update function and priority
     * This function is forwarded to the scheduler.
     */
    void add_task_to_all(update_function_type func,
                         double priority) {}
    scheduler_options unused;
    void set_scheduler_options(const scheduler_options& opts) { 
      unused = opts;
    }






  }; // end of pushy_distributed_engine

  
      
}; // end of namespace graphlab


#include <graphlab/macros_undef.hpp>

#endif
