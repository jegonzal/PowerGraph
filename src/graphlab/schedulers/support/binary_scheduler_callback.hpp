/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
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

