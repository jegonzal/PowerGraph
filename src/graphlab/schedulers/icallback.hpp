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

#ifndef GRAPHLAB_ICALLBACK_HPP
#define GRAPHLAB_ICALLBACK_HPP

#include <vector>

#include <graphlab/scope/iscope_factory.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/monitoring/imonitor.hpp>



namespace graphlab {
  template <typename Graph> class iengine;
  /**
   * This is the class that is passed to the update functions. This
   * class allows the update functions to create and add new tasks to
   * be scheduled.  The commit() function should be called by the task
   * executor at the end of the update function call.  This allows the
   * scheduler to either "buffer" the tasks and apply them only on a
   * commit, or to do some intelligent pre-pruning of tasks, or to
   * even insert the task directly into the scheduler.
   */
  template<typename Graph>
  class icallback {
  public:
    typedef Graph graph_type;
    typedef update_task<Graph> update_task_type;
    typedef iengine<Graph> iengine_type;
    typedef typename update_task_type::update_function_type 
    update_function_type;


    icallback() {}
    virtual ~icallback() {}

    /**
     * Adds a task to execute the update function on the vertex with
     * the given priority.
     */
    void add_task(vertex_id_t vertex,
                   update_function_type update_fun,
                   double priority = 1.0) {
      update_task_type task(vertex, update_fun);
      add_task(task, priority);
    }
    
    /**
     * Adds an update task with a particular priority.  The exact
     * behavior of this function depends on the type of scheduler
     */
    virtual void add_task(update_task_type task, double priority) = 0;

    /**
     * Creates a collection of tasks on all the vertices in
     * 'vertices', and all with the same update function and priority
     */
    virtual void add_tasks(const std::vector<vertex_id_t>& vertices, 
                           update_function_type func, double priority) = 0;    
                         

    /**
     * Calling this function will force the engine to abort
     * immediately
     */
    virtual void force_abort() = 0;

    

    // TODO : FIX
     virtual void enable_buffering() {}
     virtual void commit() {}
     virtual int num_of_buffered_tasks() { return 0; }
  };

}; //end graphlab namespace

#endif
