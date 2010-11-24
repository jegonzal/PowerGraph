#ifndef GRAPHLAB_SET_SCHEDULER_HPP
#define GRAPHLAB_SET_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>
#include <vector>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/synchronized_circular_queue.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/schedulers/set_scheduler/ivertex_set.hpp>
#include <graphlab/schedulers/set_scheduler/vertex_set.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/schedulers/support/unused_scheduler_callback.hpp>
#include <graphlab/schedulers/support/vertex_task_set.hpp>
#include <graphlab/schedulers/set_scheduler/set_generic_typedefs.hpp>
#include <graphlab/util/shared_termination.hpp>
#include <graphlab/schedulers/set_scheduler/execution_plan.hpp>

#include <graphlab/logger/logger.hpp>

namespace graphlab {

  /**
     This is the advanced set scheduler framework.
     Usage:

     User must specify a scheduling function of the type
     schedule_function.  The scheduling function must first define a
     set of vertex sets which are connected to each other via trigger
     functions. There is a root_set() which is the set of all
     vertices.
  */
  template<typename Graph>
  class set_scheduler : 
    public ischeduler<Graph> {

  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;
    typedef ivertex_set<Graph> ivertex_set_type;
    typedef execution_plan<Graph> execution_plan_type;

    typedef void (*schedule_function_type)(set_scheduler &sched);

  private:
    using base::monitor;

  public:
    set_scheduler(iengine_type* engine,
                  Graph& g_, 
                  size_t ncpus):
      g(g_),
      activeset(-1),
      warningprinted(false),
      executing(false),
      numcpus(ncpus),
      terminator(ncpus+1),
      callback(engine) {
      cureplan = NULL;
      nextschedpos.resize(ncpus);
      rootset.init(&g, NULL, this, ncpus);
      startschedule = new barrier(2);
    }
  
    /// destructor. deletes all vertex sets
    ~set_scheduler() { 
      for (size_t i = 0; i < vertexsets.size(); ++i) {
        delete vertexsets[i];
      }
      delete startschedule;
    }

    /// Returns a task from the active set (set by execute_rep)
    bool get_task_from_active_set(size_t cpuid, update_task_type &ret_task) {
      if (activeset == -1) return false;
      
      int localactiveset = activeset;
      
      if (localactiveset == -1) {
        return false;
      }
      if (vertexsets[localactiveset]->next(nextschedpos[cpuid], cpuid) == false) {
        return false;
      }
      size_t tmp = nextschedpos[cpuid];
      ret_task = update_task_type(tmp, schedupdate);
      return true;
    }

    /** Called by the engine to get a new task */
    sched_status::status_enum get_next_task(size_t cpuid,
                                           update_task_type &ret_task) {
      // check for pending tasks from execute()
      if (cureplan != NULL) {
        if (cureplan->get_next_task(cpuid, ret_task)) { 
          executedtaskctr.inc();
          return sched_status::NEWTASK;
        }
        else {
          schedposlock.lock();
          schedposcond.signal();
          schedposlock.unlock();
          return sched_status::EMPTY;
        }
      }
      else if(pendingtasks.safepop(&ret_task)) {
        pendingtaskctr.dec();
        executedtaskctr.inc();
        if (monitor != NULL) 
          monitor->scheduler_task_scheduled(ret_task, 0.0);
        return sched_status::NEWTASK;
      }
      // no pending tasks. check for tasks from the active set
      else if(get_task_from_active_set(cpuid, ret_task)) {
        executedtaskctr.inc();
        if (monitor != NULL) 
          monitor->scheduler_task_scheduled(ret_task, 0.0);
        return sched_status::NEWTASK;
      }
      // otherwise, are we complete?
      else if (complete) return sched_status::EMPTY;
      // If not, we are waiting
      else {
        schedposlock.lock();
        schedposcond.signal();
        schedposlock.unlock();
        return sched_status::EMPTY;
      }
    }


    /// triggers changes in the root set
    void update_state(size_t cpuid,
                      const std::vector<vertex_id_t> &updated_vertices,
                      const std::vector<edge_id_t> &updated_edges) {
      if (rootsetevents & MODIFY_VERTEX_EVENT) {
        for (size_t i = 0;i < updated_vertices.size(); ++i) {
          rootset.modify_vertex(NULL, updated_vertices[i]);
        }
      }
      if (rootsetevents & MODIFY_EDGE_EVENT) {
        for (size_t i = 0;i < updated_edges.size(); ++i) {
          rootset.modify_edge(NULL, updated_edges[i]);
        }
      }
    }

    size_t num_cpus() const {
      return numcpus;
    }

    void scoped_modifications(size_t cpuid, vertex_id_t rootvertex,
                              const std::vector<edge_id_t>& updated_edges) {
      if (rootsetevents & MODIFY_VERTEX_EVENT) {
        rootset.modify_vertex(NULL, rootvertex);
      }
      if (rootsetevents & MODIFY_EDGE_EVENT) {
        for (size_t i = 0;i < updated_edges.size(); ++i) {
          rootset.modify_edge(NULL, updated_edges[i]);
        }
      }
    }

    void completed_task(size_t cpuid, const update_task_type& task) {
      if (cureplan != NULL) {
        cureplan->completed_task(cpuid);
      }
      finishedtaskctr.inc();
    }


    /*========================================================================*/
    /*              Functions Called by the scheduling function               */
    /*========================================================================*/
    /// schedule access function. Returns the root set
    ivertex_set_type& root_set() { return rootset; }

    /// schedule access function. Attaches a set to another set through a handler
    template <typename SetType>
    SetType& attach(const SetType& destset, ivertex_set_type &srcset) {
      if (executing) {
        logger(LOG_FATAL, "You should not attach sets once execution begins");
        static SetType __UNUSED_SET__(destset);
        return __UNUSED_SET__;
      }

      // see if the vset already exists
      typename std::vector<ivertex_set_type*>::iterator iter = find(vertexsets.begin(),
                                                           vertexsets.end(),
                                                           &srcset);
    
      assert(iter != vertexsets.end() || (&srcset == &rootset));
      // create a new vertex set
      SetType* newset = new SetType(destset);
      vertexsets.push_back(newset);
    
      newset->init(&g, &srcset, this, numcpus);
      logger(LOG_INFO, "New Set type %s with id %d", newset->name().c_str(),
             vertexsets.size() - 1);
      logger(LOG_INFO, "Connecting %s to %s", srcset.name().c_str(),
             newset->name().c_str());
      return *newset;
    }

    void start() {
      begin_schedule();
    }
    /// Called after all the sets are constructed.
    void init() {
      ss_set_type allvertices;
      for (size_t i = 0; i < g.num_vertices(); ++i) {
        ss_insert(allvertices, i);
      }
      static ss_set_type unused;
      rootset.rebuild(NULL, allvertices);
      rootset.resolve_event_handlers();
      rootset.set_as_root_set();
      rootsetevents = rootset.all_dependent_events();
      logger(LOG_INFO, "Set Scheduler Initialization Complete");
      logger(LOG_INFO, "Set Scheduler Active Event Set: %d", rootsetevents);
      startschedule->wait();
    }
  
  
    /// schedule access function. Inserts the current state of the set
    /// into the task queue for execution using update function u
    void execute(ivertex_set_type& srcset, update_function_type u) {
      executing = true;
      ss_set_type sched = srcset.get();
      ss_set_type_iterator en = begin(sched);
      ss_set_type_iterator en_end = end(sched);
      while(en != en_end) {
        pendingtaskctr.inc();
        pendingtasks.push(update_task_type(*en, u));
        ++en;
      }
      terminator.new_job();
    }

    void execute(vertex_id_t v, update_function_type u) {
      executing = true;
      pendingtaskctr.inc();
      pendingtasks.push(update_task_type(v, u));
      terminator.new_job();
    }

    void execute_plan(execution_plan_type& eplan) {
      executing = true;
      schedposlock.lock();
      eplan.init_plan();
      cureplan = &eplan;
      schedposlock.unlock();
      terminator.new_job();
      wait();
    
      schedposlock.lock();
      cureplan = NULL;
      schedposlock.unlock();

    }
  
    /// Executes update_function u repeatedly on elements in the vertex_set
    /// until the vertex set is empty.
    void execute_rep(ivertex_set_type& srcset, 
                     update_function_type u) {
      schedposlock.lock();
      for (size_t i = 0; i < numcpus; ++i) srcset.first(nextschedpos[i], i);
      schedupdate = u;
      // find the activeset
      for (size_t i = 0; i < vertexsets.size(); ++i) {
        if (vertexsets[i] == (&srcset)) {
          activeset = i;
        }
      }
      assert(activeset != -1);
      schedposlock.unlock();
      terminator.new_job();
      wait();
    
      schedposlock.lock();
      activeset = -1;
      schedupdate = NULL;
      schedposlock.unlock();
    } // end of execute_rep


    /// Waits for the schedule to become empty
    void wait() {
      schedposlock.lock();
      while(true) {  
        if ((activeset == -1 || vertexsets[activeset]->size() == 0) && 
            pendingtaskctr.value == 0 && 
            executedtaskctr.value == finishedtaskctr.value &&
            (cureplan == NULL || cureplan->done())) {
          break;
        }
        schedposcond.wait(schedposlock);
      }
      schedposlock.unlock();
      
    }// end of wait

    /*========================================================================*/
    /*                                  END                                   */
    /*========================================================================*/


    /// Launches the scheduler using the scheduling function s. 
    /// This function should be called before the engine starts
    void begin_schedule()  {
      executing = false;
      complete = false;
      size_t nschedulethreads = 1;
      //  startbarrier = new graphlab::barrier(nschedulethreads);
      worker.resize(nschedulethreads);
      for (size_t i = 0; i < nschedulethreads; ++i) {
        worker[i].func = scheduling_function;
        worker[i].parent = this;
        scheduling_thread.launch(&worker[i]);
      }
      startschedule->wait();
    }
    
    void set_schedule_function(schedule_function_type s) {
      scheduling_function = s;
    }
    void set_option(scheduler_options::options_enum opt, void* value) { 
      if (opt == scheduler_options::SCHEDULING_FUNCTION) {
        set_schedule_function((schedule_function_type)(value));
      }
    }


    /// The background schedule worker. This worker basically calls
    /// the schedulng function and waits for the schedule to complete execution
    class scheduler_worker : public runnable {
    public:
      schedule_function_type func;
      set_scheduler* parent;
      size_t nattachcalls;  // number of calls to attach so far
    
      scheduler_worker() {nattachcalls = 0;}

      void run() {
        func(*parent);
        parent->wait();
        while(!parent->complete) {
          parent->terminator.begin_sleep_critical_section(parent->numcpus);
          if (parent->terminator.end_sleep_critical_section(parent->numcpus)) 
            break;
        }
      }
    };


 
    /// Whether this scheduler wants to know changes of scope
    bool observes_scope_changes() {
      return true;
    }
   
    bool need_vertex_locks() {
      return true;
    }

    /// UNUSED
    callback_type& get_callback(size_t cpuid) { return callback; }

    /// UNUSED
    void add_task(update_task_type task, double priority) { 
      print_warning(); 
    }

    /// UNUSED
    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority) { print_warning(); }

    /// UNUSED
    void add_task_to_all(update_function_type func, double priority) { 
      print_warning(); 
    }
    

    bool completed() const {
      return complete;
    }

 
    Graph& get_graph() {
      return g;
    }
  private:

    Graph& g;     ///  A reference to the engine's graph
  
    
    ///  A collection of all vertex sets constructed
    std::vector<ivertex_set_type*> vertexsets;


    /// The root vertex set comprising of all vertices
    vertex_set<Graph> rootset;
  
  
    //======= Datastructures to manage get_next_task =================
    /// this manages stuff added by execute(ivertex_set)
    synchronized_circular_queue<update_task_type> pendingtasks;
  
  
    //====Datastructures to manages stuff added by execute_rep(vertex_set) ====
    // set scheduling datastructures
    int activeset; ///The active set scheduled in schedule_rep
    std::vector<vertex_id_t> nextschedpos;    ///The position of the "next" iterator
    update_function_type schedupdate;  ///The update function to use
    mutex schedposlock; /// A spinlock on the next iterator
    conditional schedposcond; 
    execution_plan_type* cureplan; /// execution plan
  
    /** This class prints a warning if the "unused functions" are called. 
        the warnings will only be printed once*/
    bool warningprinted;
  
    /// Whether the scheduler is running
    bool executing;
  
    /// Whether the scheduler is complete
    bool complete;
  
    size_t numcpus; /// number of processors
    size_t rootsetevents; 

    thread_group scheduling_thread; /// The worker running the schedule
    std::vector<scheduler_worker> worker;  /// The worker running the schedule

    barrier *startschedule;
    /** barrier to start all threads simultaneously. This barrier is used twice.
        Once to wait for the master thread to finish creating all the sets, 
        And second at the end of init() */
  
  
    atomic<size_t> pendingtaskctr;  /// A counter of the number of tasks pending
    atomic<size_t> executedtaskctr; /// A counter of the number of tasks executed
    atomic<size_t> finishedtaskctr; /// A counter of the number of tasks completed

    schedule_function_type scheduling_function;

    
    
    shared_termination terminator;

    /// The unused callback returned by get_callback
    unused_scheduler_callback<Graph> callback;

    
    void print_warning() {
      if (!warningprinted) {
        logger(LOG_ERROR, "set scheduler does not support add_task");
      }
      warningprinted = true;
    }
  }; 


} // end of namespace graphlab

#endif
