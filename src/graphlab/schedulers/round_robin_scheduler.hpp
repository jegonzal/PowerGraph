/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * This class defines a basic round robin scheduler.
 * written by Danny Bickson
 **/
#ifndef RR_SCHEDULER_HPP
#define RR_SCHEDULER_HPP

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
#include <graphlab/schedulers/support/unused_scheduler_callback.hpp>
#include <graphlab/util/synchronized_circular_queue.hpp>
#include <graphlab/util/controlled_termination.hpp>


#include <graphlab/macros_def.hpp>

namespace graphlab {



  /**
     \ingroup group_schedulers
     This class defines a simple First-In-First_Out scheduler
  */
  template<typename Graph>
  class round_robin_scheduler: public ischeduler<Graph> {
  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;
    typedef controlled_termination terminator_type;
  
  private:
    size_t numvertices; /// Remember the number of vertices in the graph
    atomic<size_t> cur_task;  /// Which vertex I am executing now.
    unused_scheduler_callback<Graph> callback;
    std::vector<update_task_type> task_set;           /// collection of tasks
    atomic<size_t> iterations;           /// Number of complete iterations so far
    size_t maxiterations;               /// maximum number of iterations
    size_t startvertex;
    size_t endtask;     /// last vertex to run
    size_t endvertex;
    controlled_termination terminator;
    size_t step;
    
  public:
    round_robin_scheduler(iengine_type* engine,
                          Graph &g, size_t ncpus):
      numvertices(g.num_vertices()),
      callback(engine),
      task_set(numvertices),
      maxiterations(0),
      startvertex(0),
      endtask(numvertices),
      endvertex(numvertices),
      step(1){
      cur_task.value = size_t(0);
    }

 
    ~round_robin_scheduler() { }

    callback_type& get_callback(size_t cpuid) {
      return callback;
    }

    void start() {
      cur_task.value = startvertex;
      endtask = startvertex + maxiterations * numvertices - 1;
      endvertex = (startvertex + numvertices - step) % numvertices;
      std::cout << "max iterations = " << maxiterations << std::endl;
      std::cout << "step = " << step << std::endl;
    };
    
    void completed_task(size_t cpuid, const update_task_type &task){ };



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
      for (vertex_id_t vertex = 0; vertex < numvertices; ++vertex){
        add_task(update_task_type(vertex, func), priority);
      }
    }


    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority) {
      foreach(vertex_id_t vertex, vertices) {
        assert(vertex >=0 && vertex < numvertices);
        task_set[vertex] = update_task_type(vertex, func);
      } 
    }
  
    void set_start_vertex(size_t v){
      logstream(LOG_INFO) << "Round robin: Starting from " << v << std::endl;    
      startvertex = v;

    }

 
    /** Get the number of iterations that the scheduler has run */
    size_t get_iterations() {
      return iterations.value;
    }

    /** Get the next element in the queue */
    sched_status::status_enum get_next_task(size_t cpuid,
                                            update_task_type &ret_task) {
      while(1) {
        size_t oldtaskvid = cur_task.inc_ret_last(); //DB: we want the increment to happen later
        size_t taskvid = oldtaskvid  % numvertices;
        size_t alternate_taskvid = (taskvid * step) % numvertices;
        if (maxiterations != 0) {
          
          if (oldtaskvid > endtask) {
            terminator.complete();
            return sched_status::EMPTY;
          }
          if (taskvid == endvertex) {
            iterations.inc();
          }
        }
        
        ret_task = task_set[alternate_taskvid];
        if (ret_task.vertex() == vertex_id_t(-1)) continue;
        assert(ret_task.vertex() == alternate_taskvid);
        return sched_status::NEWTASK;
      }
    }

    terminator_type& get_terminator() {
      return terminator;
    };

    void set_options(const scheduler_options &opts) {
      opts.get_int_option("max_iterations", maxiterations);
      opts.get_int_option("start_vertex", startvertex);
      opts.get_int_option("step", step);
    }

    static void print_options_help(std::ostream &out) {
      out << "max_iterations = [integer, default = 0]\n";
      out << "start_vertex = [integer, default = 0]\n";
      out << "step = [integer which is either 1 or relatively prime to #vertices, default = 1]\n";
    }

    
    
  };

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

