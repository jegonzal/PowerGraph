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

#ifndef BINARY_SCHEDULER_CALLBACK_HPP
#define BINARY_SCHEDULER_CALLBACK_HPP

#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/schedulers/icallback.hpp>
#include <graphlab/logger/logger.hpp>

namespace graphlab {
  /**
     This is a special callback which only records whether or not
     add_tasks was called
  */
  template<typename Graph>
  class binary_scheduler_callback : 
    public icallback<Graph> {
 
    typedef Graph graph_type;
    typedef icallback<Graph> base;
    typedef ischeduler<Graph> scheduler_type;

    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
 
  public:
    binary_scheduler_callback() {add_task_called= false; }
    virtual ~binary_scheduler_callback() {}
  
    void add_task(update_task_type task, double priority) {
      add_task_called = true;
    }
    /** Creates a collection of tasks on all the vertices in
        'vertices', and all with the same update function and
        priority  */
    void add_tasks(const std::vector<vertex_id_t>& vertices, 
                   update_function_type func, double priority) {
      add_task_called = true;
    }
  
  
    void reset() {
      add_task_called = false;
    }

    void force_abort() {
      assert(false); // Unsupported
    }
    
    /// Commits the tasks added
    void commit() { };
    bool add_task_called;
  };

}
#endif
