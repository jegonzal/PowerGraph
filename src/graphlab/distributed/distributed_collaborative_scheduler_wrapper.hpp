#ifndef DISTRIBUTED_COLLAB_SCHEDULER_WRAPPER_HPP
#define DISTRIBUTED_COLLAB_SCHEDULER_WRAPPER_HPP
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/distributed/distributed_terminator.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/schedulers/support/binary_vertex_task_set.hpp>


/**
  * Distributed wrapper for scheduler. Unlike the basic
  * distributed_scheduler_wrapper, this schedules a task
  * only when all tasks that were spun by the same task
  * previously have been executed.
  * TODO
  *  - support multiple functions
  *  - support priority
  */

namespace graphlab {


  // forward declaration
  template<typename Graph, typename Scheduler> class distributed_engine;

  template<typename Graph, typename Scheduler>
  class distributed_collaborative_scheduler_wrapper: public ischeduler<Graph>{
   public:
    typedef Graph graph_type;
    typedef update_task<Graph> update_task_type;
    typedef typename update_task_type::update_function_type 
                                        update_function_type;
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
    
    /* Collaborative bookkeeping. Supports only one task now. */
    std::vector< bool > book_is_scheduled;
    std::vector< atomic<int> > book_waiting_to_execute;
    std::vector< std::vector< vertex_id_t > > book_waitlist;
    std::vector< vertex_id_t > curvertex;
	
	atomic<size_t> numpendingtasks, num_handledpendingtask;
	
	std::vector< spinlock >  booklock;
	
	   
		
   public:     
    distributed_collaborative_scheduler_wrapper(distributed_control &dc,
                                  Graph& graph,
                                  size_t ncpus):
                                  ischeduler<Graph>(graph, ncpus),
                                  dc(dc),
                                  dterm(dc),
                                  sched(graph, ncpus),
                                  graph(graph),
                                  callbacks(ncpus, direct_callback<Graph>(this)),
                                  popped_but_incomplete_tasks(graph.num_vertices()),
                                  book_is_scheduled(graph.num_vertices(), false),
                                  book_waiting_to_execute(graph.num_vertices(), 0),
                                  book_waitlist(graph.num_vertices(), std::vector<vertex_id_t>())
                                  {
      sched_wrapper_target = this;
      curvertex.resize(ncpus);
      booklock.resize(graph.num_vertices());
    }
    
    /// destructor
    virtual ~distributed_collaborative_scheduler_wrapper() {};
        
        
   
    void start() { sched.start(); };

    static distributed_collaborative_scheduler_wrapper<Graph, Scheduler> *sched_wrapper_target;
    
    static void add_task_handler(distributed_control& dc,
                                procid_t source,
                                void* ptr,   // ptr is a vertex_id_t list
                                size_t len,
                                update_function_type func,
                                double priority) {
      size_t numel = len / sizeof(vertex_id_t);
      vertex_id_t *vidlist = (vertex_id_t*)(ptr);
	

      for (size_t i = 0;i < numel; ++i) {
        vertex_id_t newvid = vidlist[i];
        update_task_type newtask(newvid, func);
        // if it is not popped an incomplete
        sched_wrapper_target->booklock[newvid].lock();

        if (sched_wrapper_target->popped_but_incomplete_tasks.get(newtask) == false) {
			// is it ready for immediate execution?
			if (sched_wrapper_target->book_waiting_to_execute[newvid].value <= 0) {
				sched_wrapper_target->sched.add_task(newtask, 1.0);
			} else {
				if (sched_wrapper_target->book_is_scheduled[newvid] == false) {
					sched_wrapper_target->numpendingtasks.inc();
					sched_wrapper_target->book_is_scheduled[newvid] = true;
				}
			}
        }
        sched_wrapper_target->booklock[newvid].unlock();

      }
      sched_wrapper_target->numtasksreceived.inc();
    }
    
    
    void  task_finished(vertex_id_t vid, update_function_type func) {
    	 booklock[vid].lock();
    	 std::vector<vertex_id_t> tmplist = book_waitlist[vid];
    	 book_waitlist[vid].clear();
		 booklock[vid].unlock();

         std::vector<vertex_id_t>::iterator it = tmplist.begin();
          // Decrease counters of all tasks waiting for us.
          while(it != tmplist.end()) {
               
             booklock[*it].lock();

          	 if (book_waiting_to_execute[ *it].dec() == 0) {
          	 	if (book_is_scheduled[*it] == true) {
          	 		 book_is_scheduled[*it] = false;

  					 num_handledpendingtask.inc();
          	 		 sched.add_task(update_task_type(*it, func), 1.0);
          	 	} 
          	 }
          	 booklock[*it].unlock();

          	 //printf("%d %d %d\n", *it, book_waiting_to_execute[ *it], book_waitlist[vid].size());
          	 it++;
          }

    }
    
    static void task_finished_handler(distributed_control& dc,
                                        procid_t source,
                                        void* ptr,   
                                        size_t len,
                                        vertex_id_t vid,
                                        update_function_type func) {

         sched_wrapper_target->task_finished(vid, func);
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
       // Prune remote task. Note: this is not consistent - but perhaps still ok!

      vertex_id_t adding_vertex = curvertex[thread::thread_id()];
      vertex_id_t vertexid = task.vertex();
      procid_t owner = graph.owner(vertexid);


	  booklock[adding_vertex].lock();
      book_waiting_to_execute[adding_vertex].inc();
      booklock[adding_vertex].unlock();

       // Add to wait list
	  booklock[vertexid].lock();
      book_waitlist[vertexid].push_back(adding_vertex);
       
      if (owner == dc.procid()) {
        if (sched_wrapper_target->popped_but_incomplete_tasks.get(task) == false) {
			// Can we run it immediatelly?
			if (book_waiting_to_execute[vertexid].value <= 0) {
				 sched.add_task(task, 1.0);
			} else {
				if (book_is_scheduled[vertexid] == false) {
					numpendingtasks.inc();
					book_is_scheduled[vertexid] = true;
				}
			}    
\
        }
      }
      else {

            	dc.remote_callx(owner,
                        distributed_collaborative_scheduler_wrapper<Graph, Scheduler>::add_task_handler,
                        &vertexid,
                        sizeof(vertex_id_t),
                        task.function(),
                        priority);
            	numtaskstransmitted.inc();
            
      }
       booklock[vertexid].unlock();

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
                        distributed_collaborative_scheduler_wrapper<Graph, Scheduler>::add_task_handler,
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
                          distributed_collaborative_scheduler_wrapper<Graph, Scheduler>::add_task_to_all_handler,
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
      if (numpendingtasks.value <= num_handledpendingtask.value && dterm.done(numtaskstransmitted.value, numtasksreceived.value)) {
        return sched_status::COMPLETE;
      }
      else {
    //   printf("%d %d\n", numpendingtasks.value, num_handledpendingtask.value);
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
      curvertex[thread::thread_id()] = task.vertex();
    }
    /**
     * This is called after a task has been executed
     */
    void completed_task(size_t cpuid, 
                                const update_task_type &task) {
      sched.completed_task(cpuid, task);
      
      // Broadcast that I completed
      for (size_t i = 0; i < dc.numprocs(); ++i) {
        if (i != dc.procid()) {
          dc.remote_callx(i,
                          distributed_collaborative_scheduler_wrapper<Graph, Scheduler>::task_finished_handler,
                          NULL,
                          0,
                          task.vertex(),
                          task.function());
        }
      }
      task_finished(task.vertex(), task.function());
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
  distributed_collaborative_scheduler_wrapper<Graph, Scheduler>
        *distributed_collaborative_scheduler_wrapper<Graph, Scheduler>::sched_wrapper_target;
}
#endif