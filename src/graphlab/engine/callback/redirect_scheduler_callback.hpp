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


#ifndef GRAPHLAB_REDIRECT_SCHEDULER_CALLBACK_HPP
#define GRAPHLAB_REDIRECT_SCHEDULER_CALLBACK_HPP
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/schedulers/icallback.hpp>
#include <graphlab/logger/logger.hpp>

namespace graphlab {
  
  template<typename Graph, typename RedirectToType>
  class redirect_scheduler_callback : 
    public icallback<Graph> {
  public:
    typedef icallback<Graph> base;
    typedef typename base::vertex_id_type    vertex_id_type;
    typedef typename base::iengine_type      iengine_type;
    typedef typename base::update_task_type  update_task_type;
    typedef typename base::update_function_type update_function_type;
    RedirectToType* redirectto;
    
    
  public:
    redirect_scheduler_callback(RedirectToType* redirectobj = NULL) :
      redirectobj(redirectobj) { }
    
    virtual ~redirect_scheduler_callback() {}

    
    void add_task(update_task_type task, double priority) {
      redirectobj->add_task(task, priority);
    }

    /** Creates a collection of tasks on all the vertices in
        'vertices', and all with the same update function and
        priority  */
    void add_tasks(const std::vector<vertex_id_type>& vertices, 
                   update_function_type func, double priority) {
      redirectobj->add_tasks(vertices, func, priority);
    }

    /**
     * Force the engine to abort.
     */
    void force_abort() {
      logger(LOG_WARNING, "abort not supported");
    }


    void disable_buffering() {}

      
    /// Commits the tasks added
    void commit() { 
      logger(LOG_WARNING, "commit not supported");
    }
    
  private:
    RedirectToType* redirectobj;
  };

}
#endif

