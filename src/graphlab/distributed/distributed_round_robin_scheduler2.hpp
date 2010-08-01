/**
 * This class defines a basic distributed_round_robin_scheduler2
 * written by Danny Bickson
 **/
#ifndef DISTR_RR2_SCHEDULER_HPP
#define DISTR_RR2_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>
#include <vector>
#include <climits>


#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/schedulers/support/vertex_task_set.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/util/synchronized_circular_queue.hpp>
#include <graphlab/util/task_count_termination.hpp>


#include <graphlab/macros_def.hpp>

namespace graphlab {

 
  template<typename Graph>
class distributed_round_robin_scheduler2: public ischeduler<Graph> {
 public:
   typedef Graph graph_type;
   typedef ischeduler<Graph> base;

    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;

    std::vector< vertex_id_t > checkpoints;
  std::list< vertex_id_t > checkpoint_iter;
  bool checkpoint_on;
  mutex checkpoint_lock;
  conditional checkpoint_cond;
  
    distributed_control * dc;

  
 private:
  size_t numvertices; /// Remember the number of vertices in the graph
  atomic<size_t> cur_task;  /// Which vertex I am executing now.
  std::vector< direct_callback<Graph> > callbacks;  /// callback for each processor
  std::vector<update_task_type> task_set;           /// collection of tasks
  atomic<size_t> iterations;           /// Number of complete iterations so far
  size_t maxiterations;               /// maximum number of iterations

  task_count_termination terminator;
  std::vector< vertex_id_t >  myvertices;

  double barrier_wait_secs;
  //================================
  //TODO GREAT BIG UGLY HACK
  //================================
  void (*barrier_fn)(void);
 public:
  distributed_round_robin_scheduler2(Graph &g, size_t ncpus):
          base(g, ncpus),
          numvertices(g.num_vertices()),
          callbacks(ncpus, direct_callback<Graph>(this)),
          task_set(g.num_vertices()) {
          
          cur_task.value = size_t(0);
          barrier_wait_secs = 0.0;
          
          maxiterations = 0;
          checkpoint_on = false;
          barrier_fn  = NULL;
          // Sort vertices
          myvertices = std::vector<vertex_id_t>(g.my_vertices());
          std::sort(myvertices.begin(), myvertices.end());
   }

 
  ~distributed_round_robin_scheduler2() { }

  callback_type& get_callback(size_t cpuid) {
    return callbacks[cpuid];
  }


  void completed_task(size_t cpuid, const update_task_type &task){ };

  void abort(){ terminator.abort(); }

  void restart() { terminator.restart(); }
  

  void set_max_iterations(size_t maxi) {
    maxiterations = maxi;
  }

  void add_task(update_task_type task, double priority) {
    if (task_set[task.vertex()].function() != NULL) {
      logger(LOG_WARNING, 
             "Adding task on vertex %d where a task already exists", 
              task.vertex());
    }
    task_set[task.vertex()] = update_task_type(task.vertex(), task.function());
  }

  void add_task_to_all(update_function_type func, double priority) {
     printf("Add task to all...\n");
    for (unsigned int i = 0; i < myvertices.size(); ++i){
      add_task(update_task_type(myvertices[i], func), priority);
    }
      printf("Add task to all... FINISH \n ");

  }


  void add_tasks(const std::vector<vertex_id_t> &vertices,
                                update_function_type func,
                                double priority) {
    foreach(vertex_id_t vertex, vertices) {
      task_set[vertex] = update_task_type(vertex, func);
    } 
  }
  
  void set_start_vertex(size_t v){
  
    assert(false); // Not supported yet
 }

 
 void set_barrier_function(void(*bar)(void)){
    barrier_fn = bar;
 }
  int get_iterations(){
     return (int)iterations.value;
  }
 
  void set_option(scheduler_options::options_enum opt, void* value) { 
    if (opt == scheduler_options::MAX_ITERATIONS) {
      set_max_iterations((size_t)(value));
    }
    else if (opt == scheduler_options::START_VERTEX) {
      set_start_vertex((size_t)(value));
    }
    else if (opt == scheduler_options::BARRIER) {
       vertex_id_t vid = (* (vertex_id_t *) (value));
      if (vid > myvertices[myvertices.size() - 1]) vid = myvertices[myvertices.size() - 1];
       if (checkpoints.size() > 0) {
        // Check points must be inserted in ascending order
        ASSERT_TRUE(vid > checkpoints[checkpoints.size()-1]);
      }
      checkpoints.push_back(vid);
      checkpoint_iter.push_back(vid);
    } else if (opt == scheduler_options::DISTRIBUTED_CONTROL) { 
        dc = (distributed_control *) value;
    } else {
      logger(LOG_WARNING, 
             "Round Robin Scheduler was passed an invalid option %d", 
             opt);
    }
  };

  /** Get the next element in the queue */
  sched_status::status_enum get_next_task(size_t cpuid, update_task_type &ret_task) {
    if (terminator.is_aborted()) return sched_status::COMPLETE;
    
         
   if (checkpoint_on) {
       return sched_status::WAITING;
   }
      
  
    size_t i = cur_task.inc() - 1; //DB: we want the increment to happen later
    i = i % myvertices.size();
    
    if (i == 0) {
          // Restore check points for next round.
          foreach(vertex_id_t w, checkpoints) checkpoint_iter.push_back(w);
    }
    
    if (maxiterations != 0) {
      if (iterations.value >= maxiterations) {
        return sched_status::COMPLETE;
      }
      if (i == myvertices.size()-1) {
          iterations.inc();
      }
    }

    vertex_id_t vid = myvertices[i];
  
   // Checkpoint checking
    if (checkpoint_iter.empty() == false) {
        checkpoint_lock.lock();
        if (checkpoint_iter.empty() == false) {
        vertex_id_t next_checkpoint = checkpoint_iter.front();
        if (vid >= next_checkpoint) {
            checkpoint_iter.pop_front();
            checkpoint_on = true;
            checkpoint_lock.unlock();
            
            timer t;
            t.start();
            
            printf("===== %d Activated checkpoint: %ld (%ld) %ld\n", (int) dc->procid(), (long int) vid,  (long int) i,  (long int) next_checkpoint);
            if (barrier_fn == NULL) dc->barrier();
            else barrier_fn();
            printf("===== %d Passed checkpoint: %ld (%lf secs)\n", (int) dc->procid(), (long int) next_checkpoint, t.current_time());
  
            checkpoint_lock.lock();
          
          // Release locals. 
            checkpoint_on = false;
                         
            barrier_wait_secs += t.current_time();
            distributed_metrics::instance(dc)->set_value("rr_barrier_wait", barrier_wait_secs);
            checkpoint_lock.unlock();


          } else checkpoint_lock.unlock();
         } else checkpoint_lock.unlock();
     }
    
    
    ret_task = task_set[vid];     
    return sched_status::NEWTASK;
  }

};

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
