#ifndef DISTRIBUTED_ENGINE_HPP
#define DISTRIBUTED_ENGINE_HPP
#include <cassert>
#include <cstdio>

#include <boost/unordered_map.hpp>
#include <boost/bind.hpp>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/synchronized_unordered_map.hpp>
#include <graphlab/util/synchronized_unordered_map2.hpp>
#include <graphlab/distributed/graph_lock_manager.hpp>
#include <graphlab/distributed/distributed_scope.hpp>
#include <graphlab/distributed/distributed_scheduler_wrapper.hpp>
#include <graphlab/distributed/distributed_shared_data.hpp>
#include <graphlab/distributed/metrics/distributed_metrics.hpp> 

#include <graphlab/macros_def.hpp>

namespace graphlab {
  
 
  /**
   * This class defines a simple multi threaded engine. 
   */
  template<typename Graph, typename DistributedScheduler>
  class distributed_engine : 
    public iengine<Graph>, public iscope_factory<Graph> {

  public:   
    typedef iengine<Graph> base;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::ischeduler_type ischeduler_type;
    typedef typename base::imonitor_type imonitor_type;
    typedef typename base::termination_function_type termination_function_type;
    typedef typename base::ishared_data_manager_type ishared_data_manager_type;
    typedef distributed_scope<Graph, graph_lock_manager<Graph> > scope_type;
    typedef std::pair<mutex, std::queue<update_task_type> > locked_pending_queue_type;
    typedef synchronized_unordered_map2<locked_pending_queue_type> pending_tasks_map_type;
  private:
    
    /** The internal worker thread class */
    class task_worker {      
      distributed_engine<Graph, DistributedScheduler>* engine;
      size_t workerid;

    public:
      task_worker() : engine(NULL), workerid(0) {  }
      void init(distributed_engine<Graph, DistributedScheduler>* _engine, size_t _workerid) {
        engine = _engine;
        workerid = _workerid;
      } // End of init      
      void run() {
        assert(engine != NULL);
        logger(LOG_INFO, "Worker %d started\n", workerid);        
        /* Start consuming tasks */
        while(true) {
          bool executed_task = engine->run_next_task(workerid, engine->lock_block_id, engine->pending_tasks);
         
          // If this was nothing to execute then fail
          if (!executed_task)  break;
        }
        
        logger(LOG_INFO, "Worker %d died\n", workerid);
      }      
    }; // end of task worker
    

    /** Responsible for maintaining the schedule over tasks */
    distributed_control &dc;
    DistributedScheduler scheduler;
    /** The number of cpus to use */
    size_t ncpus;
    size_t maximum_backlog;
    /* Worker bookkeeping */
    std::vector<size_t> task_counts;
    std::vector<size_t> worker_works;
    size_t max_back_log;
    atomic<size_t> backlog;
    pending_tasks_map_type pending_tasks;
    size_t lock_block_id;
    
    /** Pointer to the data manager */
    Graph& graph;
    distributed_lock_manager<Graph> dlm;
    graph_lock_manager<Graph> lock_manager;
    distributed_shared_data<Graph>* data_manager;


    bool printed;

    timer _timer;   
    double timeout;
    bool aborted;
    volatile size_t sync_barrier_on;
    
    std::vector<termination_function_type> term_functions;
    float lasttimercheck;
    scope_range::scope_range_enum scoperange;
    
    std::vector<scope_type*> scopes;

    exec_status termination_cause;
    
  public:
    /** Initialize the multi threaded engine */
    distributed_engine(distributed_control &dc,
                       Graph& g, 
                       size_t local_num_cpus = thread::cpu_count()) :
      iscope_factory<Graph>(g, local_num_cpus),
      dc(dc),
      scheduler(this, dc, g, local_num_cpus),
      ncpus(local_num_cpus),
      maximum_backlog(100),
      task_counts(local_num_cpus,0),
      worker_works(local_num_cpus,0),
      max_back_log(0),
      backlog(0),
      pending_tasks(131071),
      graph(g), dlm(dc, g), lock_manager(dc, dlm, g), data_manager(NULL),
      scoperange(scope_range::USE_DEFAULT){

      aborted = false;
      scheduler.register_monitor(NULL);
      timeout = 0;
      lasttimercheck = 0;
      assert(local_num_cpus > 0);
      printed = false;
      sync_barrier_on = 0;
      // allocate the scopes
      scopes.resize(local_num_cpus);
      for (size_t i = 0;i < scopes.size(); ++i) {
        scopes[i] = new scope_type(lock_manager);
      }
      
      /*for (size_t i = 0 ;i < g.my_vertices().size(); ++i) {
        pending_tasks.insert(i, locked_pending_queue_type());
        }*/
      // create the lockblock
    }
    
    ~distributed_engine() {
      for (size_t i = 0;i < scopes.size(); ++i) {
        delete scopes[i];
      }
      scopes.clear();
    }
    
    /** This returns the LOCAL number of threads started */
    size_t get_ncpus() const { return ncpus; } 


    void set_default_scope(scope_range::scope_range_enum default_scope_range) {
      scoperange = default_scope_range;
    }

    void set_max_backlog(size_t _maximum_backlog) {
      maximum_backlog =_maximum_backlog;
    }

    /** register the listener */
    void register_monitor(imonitor_type* _listener) {
      ASSERT_TRUE(false);
    }
    

    /**
     * Timeout. Default - no timeout. 
     */
    void set_timeout(size_t timeout_secs) {
      timeout = timeout_secs;
    }
    


    void set_caching(bool _caching) {
      lock_manager.set_caching(_caching);
      if (_caching) lock_manager.init();
    }

    void set_constant_edges(bool _const_edges) {
      lock_manager.set_constant_edges(_const_edges);
    }

    void set_use_adjacent_vertices(bool _use_adjacent_vertices) {
      lock_manager.set_use_adjacent_vertices(_use_adjacent_vertices);
    }

    // caching must be true for this to have effect
    void set_vertex_scope_pushed_updates(bool _vertex_scope_pushed_updates) {
      lock_manager.set_vertex_scope_pushed_updates(_vertex_scope_pushed_updates);
    }


    void sync_barrier() {
      std::cout << dc.procid() << " entering sync barrier " << std::endl;
      sync_barrier_on = 1;
      while (sync_barrier_on) {
        sched_yield();
      }
      std::cout << dc.procid() << " MPI barrier " << std::endl;
      dc.mpi_barrier();
      std::cout << dc.procid() << " Leaving sync barrier" << std::endl;
    }
    /**
     * Task budget - max number of tasks to allow
     */
    virtual void set_task_budget(size_t max_tasks) {
      ASSERT_TRUE(false);
    }

    bool check_all_terminators() {
      if(data_manager != NULL) {        
        for (size_t i = 0;i < term_functions.size();++i) {
          if (term_functions[i](*data_manager)) return true;
        }
        return false;
      } else if (term_functions.size() > 0) {

        logger(LOG_WARNING,
               "Assessing termination without a data_manager!");
      }
      return false;
    }
        
    /**
     * Called by workers. out_of_work flag is passed as a pointer so
     * task's employment status can be updated immediatelly.  This
     * function must be thread safe.
     */
    bool run_next_task(int cpuid, size_t lock_block_id, pending_tasks_map_type &pending_tasks) {
      if (aborted) return false;
      /* Check for timeout (only cpuid 0 checks) */
      if (timeout > 0 && cpuid == 0) {
        if (_timer.current_time() > timeout) {
          aborted = true;
        }
      }
      bool didstuff = false;
      
      if (dc.procid() == 0 && lowres_time_seconds() - lasttimercheck >= 0.1) {
        if(data_manager != NULL) data_manager->signal_all();
        lasttimercheck = lowres_time_seconds();  
        //logger(LOG_INFO, "MaxBacklog: %d", max_back_log[cpuid]);
      }
      /////////////////////////////////////////////////////////////
      // Phase 3 check distributed shared data for tasks         //
      /////////////////////////////////////////////////////////////
      

      if (data_manager != NULL) {
        if (data_manager->progress(cpuid)) {
          didstuff = true;
        }
      }

      if (sync_barrier_on) {
        if (data_manager->has_pending_tasks() == false) { 
          sync_barrier_on = 0;
          std::cout << "Sync barrier flag cleared" << std::endl;
          return true;
        }
      }

      max_back_log = std::max(max_back_log, (size_t)(backlog.value));
      
      /////////////////////////////////////////////////////////////
      // Phase 1 first check if any task is ready to be executed //
      /////////////////////////////////////////////////////////////
      if (backlog.value > 0) {
        dist_scope_request scopereqdone;
        lock_manager.block_status(lock_block_id, scopereqdone);
        while (scopereqdone.vertex != vertex_id_t(-1)) {
          // we have work to do!
          
          scope_type* scope = get_distributed_scope(cpuid, scopereqdone.vertex);
          
          DCHECK_NE(scope, NULL);
          // get the task we are supposed to run
          pending_tasks.write_critical_section(scopereqdone.vertex);
          std::pair<bool, locked_pending_queue_type*> pt = pending_tasks.find(scopereqdone.vertex);
          if (pt.first == false) {
            ASSERT_MSG(false, "Successful scope lock without a pending task");
          }
          locked_pending_queue_type & pendingqueue = *(pt.second);
          update_task_type task = pendingqueue.second.front();
          pendingqueue.second.pop();
          if (pendingqueue.second.size() == 0) {
            pending_tasks.erase(scopereqdone.vertex);
          }
          pending_tasks.release_critical_section(scopereqdone.vertex);
          // Update task counts and "work". Work is indegree+outdegree          
          task_counts[cpuid] = task_counts[cpuid]++;
          worker_works[cpuid] += scope->in_edge_ids().size() +
            scope->out_edge_ids().size();
          

          // execute the task
          assert(task.function() != NULL);
          // NOTE: IMPORTANT! the task callback MUST immediately issue
          // and cannot delay. By socket streaming rules, task creations 
          // MUST reach the other side, before the unlock happens or we
          // will lose causal consistency
          scheduler.started_task(cpuid, task);
          typename DistributedScheduler::callback_type & cback = scheduler.get_callback(cpuid);
          task.function()(*scope, cback, data_manager);
          
          // we never commit this scope. This scope is managed through a block
          
          lock_manager.block_release_partial(lock_block_id, scopereqdone);
          backlog.dec();
          cback.commit();
          scheduler.completed_task(cpuid, task);
          didstuff = true;
          lock_manager.block_status(lock_block_id, scopereqdone);
        }
      }
      /////////////////////////////////////////////////////////////
      // Phase 2 check scheduler for new tasks to do             //
      /////////////////////////////////////////////////////////////

      if (backlog.value < maximum_backlog) {
        
        // we now check the scheduler
        update_task_type newtask;              
        // Get the next task along with the status of the scheduler
        sched_status::status_enum stat = scheduler.get_next_task(cpuid, newtask);
          
        if (stat == sched_status::COMPLETE && backlog.value == 0) {
          logger(LOG_INFO, "%d Complete", dc.procid());
          logger(LOG_INFO, "Backlog: %d", backlog.value);
          logger(LOG_INFO, "MaxBacklog: %d", max_back_log);
          lock_manager.print_stats();
          /* Scheduler says we terminate */
          return false;
        }
        else if (stat == sched_status::NEWTASK) {
          didstuff = true;
          const vertex_id_t vertex = newtask.vertex();
          // There is a new task to schedule 
          // Ensure that a valid task function was passed
          assert(newtask.function() != NULL);
          // put into the queue
          pending_tasks.write_critical_section(vertex);
          std::pair<bool, locked_pending_queue_type*> pt = pending_tasks.find(vertex);
          if (pt.first == false) {
            //            logger(LOG_WARNING, "Adding remote task. Performance could be adversely affected");
            // create the task queue
            pt = pending_tasks.insert_with_failure_detect(vertex, locked_pending_queue_type());
          }
          locked_pending_queue_type & pendingqueue = *(pt.second);
          pendingqueue.second.push(newtask);
          pending_tasks.release_critical_section(vertex);
          // issue a deferred lock for this task
          lock_manager.block_add_deferred_lock(lock_block_id,
                                               dist_scope_request(vertex, scoperange));
          backlog.inc();
        
        }
        else if (stat == sched_status::WAITING) {
          //logger(LOG_INFO, "wait");
          didstuff = false;
        }
      }
      //std::cout << ".";
      //std::cout.flush();
      /*   else {
           logger(LOG_INFO, "backlogged");
           }*/
   
   

      if (didstuff == false) {
        //logger(LOG_INFO, "yield");
        sched_yield();
      }
      return true;
    } // end of run_next_task

   
    void set_sched_yield(bool value) {
      logger(LOG_INFO, "distributed engine does not support set_sched_yield()");    
    }

    void set_cpu_affinities(bool value) {
      logger(LOG_INFO, "distributed engine does not support set_cpu_affinities()");    
    }
   
    /** get a reference to the scheduler */
    ischeduler_type& get_scheduler() { return scheduler; }

    void set_shared_data_manager(ishared_data_manager_type* manager) {
      data_manager = dynamic_cast<distributed_shared_data<Graph>*>(manager);

      if(data_manager != NULL) {
        data_manager->set_scope_factory(this);
        data_manager->set_lock_manager(&graph, &lock_manager);
      }
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
    
      // First do a bandwith test. TODO: remove
      distributed_metrics::instance(&dc)->execute_bandwith_test();
    
      //! Finalize the graph (this could take a while so you should do
      //! it before calling start for timing purposes)
      termination_cause = EXEC_TASK_DEPLETION;
      
      aborted = false;
      /* Timing */
      _timer.start();
      lasttimercheck = lowres_time_seconds();
      /* Enable scheduler to clean up in restarts */
      dc.barrier();
      scheduler.start();
      dc.barrier();
      /* Initialize a pool of threads */
      std::vector<task_worker> workers(ncpus);
      thread_group threads;
      lock_block_id = lock_manager.block_deferred_lock(std::vector<dist_scope_request>(), 2*maximum_backlog);
      for(size_t i = 0; i < ncpus; ++i) {
        // Initialize the worker
        workers[i].init(this, i);
        // Start the worker thread using the thread group with cpu
        // affinity attached (CPU affinity currently only supported in
        // linux) since Mac affinity is set through the NX frameworks
        threads.launch(boost::bind(&task_worker::run, &(workers[i])));
      }
      /* Wait for all threads to return */
      logger(LOG_INFO, "Wait until finished...");
      threads.join();
     
      double running_time = _timer.current_time();
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
      fprintf(F, "engine=\"pull\"\nexecution_time=%lf\nncpus=%d\ntaskcount=%ld\nwork=%ld\nresidual=0\nmemory_writes_mb=0\nmemory_reads_mb=0\n", running_time, (int) ncpus,total_counts, total_work);
      fclose(F);
      /// ===== END STATS OUTPUT ===== ///
      distributed_metrics::instance(&dc)->set_value("execution_time", running_time); 
      distributed_metrics::instance(&dc)->set_value("taskcount", total_counts);
      distributed_metrics::instance(&dc)->set_value("ncpus", ncpus);
      dc.report_stats();

      if (!aborted) {
        termination_cause =  EXEC_TIMEOUT;
      }
      
    } // end of start



    void stop() {
      // Unsupported:
      assert(false);
    }
    

    exec_status last_exec_status() const { 
      return termination_cause;
    }


    /**
     * Get the total number of updates
     */
    size_t last_update_count() const { 
      size_t total_updates = 0;
      foreach(size_t tc, task_counts) {
        total_updates += tc;
      }
      return total_updates;
    } // end of get profiling info
    

    /** It is not safe to call get_distributed_scope on a vertex without a lock on the 
        vertex */
    scope_type* get_distributed_scope(size_t cpuid,
                                      vertex_id_t vertex) {
      ASSERT_LE(cpuid, scopes.size());
      scopes[cpuid]->init(&graph, vertex);
      return scopes[cpuid];
    }
    
    
    // these functions allow the distributed engine to also
    // satisfy the scope_factory interface
    size_t num_vertices() const{
      return graph.num_vertices();
    }
    
    void release_scope(iscope<Graph>* s){
      // UNUSED
    }
    
    iscope<Graph>* get_scope(size_t cpuid,
                             vertex_id_t vertex,
                             scope_range::scope_range_enum s = scope_range::USE_DEFAULT) {
      return get_distributed_scope(cpuid, vertex);
    }



    /////////////////Stuff not implemented ////////////////////////
    /**
     * Adds an update task with a particular priority.
     * This function is forwarded to the scheduler.
     */
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
                         double priority) {};

    scheduler_options unused;
    void set_scheduler_options(const scheduler_options& opts) { 
      unused = opts;
    }

  
  }; // end of distributed_engine

  
      
}; // end of namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif
