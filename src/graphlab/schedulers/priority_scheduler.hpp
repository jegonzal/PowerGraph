/**
 * This class defines a basic priority scheduler.
 **/
#ifndef GRAPHLAB_PRIORITY_SCHEDULER_HPP
#define GRAPHLAB_PRIORITY_SCHEDULER_HPP

#include <cmath>
#include <cassert>
#include <queue>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/schedulers/support/vertex_task_set.hpp>
#include <graphlab/util/task_count_termination.hpp>

//#include <bitmagic/bm.h>

#include <graphlab/macros_def.hpp>
namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Graph>
  class priority_scheduler : 
    public ischeduler<Graph> {

  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;
    typedef task_count_termination terminator_type;

  private:
    using base::monitor;

  private:
    /** Remember the number of vertices in the graph */
    size_t num_vertices;  
    
    /** The queue over vertices */
    mutable_queue<size_t, double> task_queue;

    /** The lock on the priority queue */
    spinlock queuelock; 

    /** The vertex task set which maintains the tasks at each
        vertex */
    vertex_task_set<Graph> task_set;
  
    /** The callbacks pre-created for each cpuid */
    std::vector<direct_callback<Graph> > callbacks;
    
    /** Used to assess termination */
    task_count_termination terminator;
    
    
  public:
    
    priority_scheduler(iengine_type* engine,
                       Graph &g, 
                       size_t ncpus) :
      num_vertices(g.local_vertices()),
      task_set(g.local_vertices()),
      callbacks(ncpus, direct_callback<Graph>(this, engine) ) { }
    

    ~priority_scheduler() { }
    
    callback_type& get_callback(size_t cpuid) {
      return callbacks[cpuid];
    }
    

    void start() {};

    /** Get the next element in the queue */
    sched_status::status_enum get_next_task(size_t cpuid, update_task_type& ret_task) {
      queuelock.lock();
      if (task_queue.empty()) {
        queuelock.unlock();
        return sched_status::EMPTY;
      } else {
        // Get the highest priority vertex
        vertex_id_t vertex = task_queue.pop().first;

        // From the task set get the highest priority task associated
        // with that vertex
        double priority = 0; 
        bool success = task_set.pop(vertex, ret_task, priority);
        assert(success); 

        // Update the priority queue with the new value for the vertex
        update_task_type new_top_task;
        double new_priority(0);
        if(task_set.top(vertex, new_top_task, new_priority)) {
          task_queue.insert_max(vertex, new_priority);
        }
        queuelock.unlock();
        if (monitor != NULL)
          monitor->scheduler_task_scheduled(ret_task, priority);
        return sched_status::NEWTASK;
      }
    } // end of get next task
    
    
    void add_task(update_task_type task, double priority) {
      bool first_add(false);
      queuelock.lock();
      // Try and add the task to the task_set
      if(task_set.add(task, priority)) {
        // This was a unique add so then increment the terminator
        terminator.new_job();
        first_add = true;
      } 
      // Update the priority queue
      vertex_id_t vertex = task.vertex();
      task_queue.insert_max(vertex, priority);
      queuelock.unlock();
      // Update any listeners
      if(monitor != NULL) {
        if(first_add) {
          monitor->scheduler_task_added(task, priority);
        } else {
          monitor->scheduler_task_pruned(task);
        } // end of inner if
      } // end of outer if
    } // end of add_task
    

    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority) {
      foreach(vertex_id_t vertex, vertices) {
        add_task(update_task_type(vertex, func), priority);
      }
    } // end of add tasks 


    
    
    void add_task_to_all(update_function_type func, double priority) {
      for (vertex_id_t vertex = 0; vertex < num_vertices; ++vertex){
        add_task(update_task_type(vertex, func), priority);
      }
    } // end of add tasks to all

    void completed_task(size_t cpuid, const update_task_type &task) {
      terminator.completed_job();
    }

    terminator_type& get_terminator() {
      return terminator;
    };

    void set_options(const scheduler_options &opts) { }

    static void print_options_help(std::ostream &out) { };

  }; // end of priority_queue class


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
