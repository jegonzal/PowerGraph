#ifndef DISTRIBUTED_SCHEDULER_WRAPPER_HPP
#define DISTRIBUTED_SCHEDULER_WRAPPER_HPP
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/distributed/distributed_terminator.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/schedulers/support/binary_vertex_task_set.hpp>
namespace graphlab {


  // forward declaration
  template<typename Graph, typename Scheduler> class distributed_engine;

  template<typename Graph, typename Scheduler>
  class distributed_scheduler_wrapper: public ischeduler<Graph>{
  public:
    typedef Graph graph_type;
    typedef update_task<Graph> update_task_type;
    typedef typename update_task_type::update_function_type 
    update_function_type;
    typedef iengine<Graph> iengine_type;
    typedef iscope_factory<Graph> scope_factory_type;
    typedef icallback<Graph> callback_type;
    typedef imonitor<Graph> monitor_type;

  private:
    distributed_control &dc;
    distributed_terminator dterm;
    Scheduler sched;
    Graph &graph;
    /// The callbacks pre-created for each cpuid
    std::vector<direct_callback<Graph> > callbacks;
    binary_vertex_task_set<Graph> popped_but_incomplete_tasks;
    atomic<size_t> numtaskstransmitted;
    atomic<size_t> numtasksreceived;
    

  public:     
    distributed_scheduler_wrapper(iengine_type* engine,
                                  distributed_control &dc,
                                  Graph& graph,
                                  size_t ncpus):
      dc(dc),
      dterm(dc),
      sched(engine, graph, ncpus),
      graph(graph),
      callbacks(ncpus, direct_callback<Graph>(this, engine)),
      popped_but_incomplete_tasks(graph.num_vertices()) {
      sched_wrapper_target = this;
      for (size_t i  =0;i < callbacks.size(); ++i) callbacks[i].enable_buffering();
    }
    
    /// destructor
    virtual ~distributed_scheduler_wrapper() {};
        
    Scheduler& scheduler() {
      return sched;
    }	
   
    void start() { sched.start(); };

    static distributed_scheduler_wrapper<Graph, Scheduler> *sched_wrapper_target;
    
    static void add_task_handler(distributed_control& dc,
                                 procid_t source,
                                 void* ptr,   // ptr is a vertex_id_t list
                                 size_t len,
                                 update_function_type func,
                                 double priority) {
      size_t numel = len / sizeof(vertex_id_t);
      vertex_id_t *vidlist = (vertex_id_t*)(ptr);
      for (size_t i = 0;i < numel; ++i) {
        update_task_type newtask(vidlist[i], func);
        // if it is not popped an incomplete
        if (sched_wrapper_target->popped_but_incomplete_tasks.get(newtask) == false) {
          sched_wrapper_target->sched.add_task(newtask, priority);
        }
      }
      sched_wrapper_target->numtasksreceived.inc();
    }

    static void add_task_to_all_handler(distributed_control& dc,
                                        procid_t source,
                                        void* ptr,   // unused
                                        size_t len,
                                        update_function_type func,
                                        double priority) {
      sched_wrapper_target->sched.add_tasks(sched_wrapper_target->graph.my_vertices(),
                                            func,
                                            priority);
      sched_wrapper_target->numtasksreceived.inc();
    }


    void add_task(update_task_type task, double priority) {
      vertex_id_t vertexid = task.vertex();
      procid_t owner = graph.owner(vertexid);

      if (owner == dc.procid()) {
        if (sched_wrapper_target->popped_but_incomplete_tasks.get(task) == false) {
          sched.add_task(task, priority);
        }
      }
      else {
        dc.remote_callx(owner,
                        distributed_scheduler_wrapper<Graph, Scheduler>::add_task_handler,
                        &vertexid,
                        sizeof(vertex_id_t),
                        task.function(),
                        priority);
        numtaskstransmitted.inc();
      }
    }
    

    void add_tasks(const std::vector<vertex_id_t>& vertices, 
                   update_function_type func, double priority) {
      std::vector<std::vector<vertex_id_t> > procs;
      procs.resize(dc.numprocs());
      for (size_t i = 0;i < vertices.size(); ++i) {
        vertex_id_t vertexid = vertices[i];
        procid_t owner = graph.owner(vertexid);
        // if its for myself, just add it
        if (owner == dc.procid()) {
          update_task_type newtask(vertexid, func);
          if (sched_wrapper_target->popped_but_incomplete_tasks.get(newtask) == false) {
            sched.add_task(newtask, priority);
          }
        }
        else {
          procs[owner].push_back(vertexid);
        }
      }
      
      
      // issue all the remote calls
      for (size_t i = 0; i < procs.size(); ++i) {
        if (procs[i].empty() == false) {
          dc.remote_callx(i,
                          distributed_scheduler_wrapper<Graph, Scheduler>::add_task_handler,
                          &(procs[i][0]),
                          sizeof(vertex_id_t) * procs[i].size(),
                          func,
                          priority);
          numtaskstransmitted.inc();
        }
      }
    }
    

    void add_task_to_all(update_function_type func, 
                         double priority) {
      for (size_t i = 0; i < dc.numprocs(); ++i) {
        if (i == dc.procid()) {
          sched.add_tasks(graph.my_vertices(), func, priority);
        }
        else {
          dc.remote_callx(i,
                          distributed_scheduler_wrapper<Graph, Scheduler>::add_task_to_all_handler,
                          NULL,
                          0,
                          func,
                          priority);
          numtaskstransmitted.inc();
        }
      }
    }
    
    callback_type& get_callback(size_t cpuid) {
      return callbacks[cpuid];
    }

    sched_status::status_enum get_next_task(size_t cpuid, 
                                            update_task_type &ret_task) {
      sched_status::status_enum ret = sched.get_next_task(cpuid, ret_task);
      if (ret == sched_status::COMPLETE) {
        // if complete, we still have to check the distributed terminator
        if (dterm.done(numtaskstransmitted.value, numtasksreceived.value)) {
          return sched_status::COMPLETE;
        }
        else {
          return sched_status::WAITING;
        }
      }
      else if (ret == sched_status::WAITING){      
        // no tasks left to do
        return ret;
      }
      else {
        // return a new task
        popped_but_incomplete_tasks.add(ret_task);
        return ret;
      }
    }

    void started_task(size_t cpuid,
                      const update_task_type &task) {
      popped_but_incomplete_tasks.remove(task);
    }
    /**
     * This is called after a task has been executed
     */
    void completed_task(size_t cpuid, 
                        const update_task_type &task) {
      sched.completed_task(cpuid, task);
    }


    /** Installs a listener (done by the engine) */
    void register_monitor(monitor_type* monitor_) { 
      sched.register_monitor(monitor_);
    }

    void abort() { 
      sched.abort();
    }

    void set_option(scheduler_options_enum::options_enum opt, void* value) { 
      sched.set_option(opt,value);
    }
  };

  template <typename Graph, typename Scheduler>
  distributed_scheduler_wrapper<Graph, Scheduler>
  *distributed_scheduler_wrapper<Graph, Scheduler>::sched_wrapper_target;
}
#endif
