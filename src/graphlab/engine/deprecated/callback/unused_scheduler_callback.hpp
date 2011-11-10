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


#ifndef GRAPHLAB_UNUSED_SCHEDULER_CALLBACK_HPP
#define GRAPHLAB_UNUSED_SCHEDULER_CALLBACK_HPP
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/schedulers/icallback.hpp>
#include <graphlab/logger/logger.hpp>

namespace graphlab {
  
  template<typename Graph>
  class unused_scheduler_callback : 
    public icallback<Graph> {
  public:
    typedef Graph graph_type;
    typedef icallback<Graph> base;
    typedef ischeduler<Graph> scheduler_type;

    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;

    
    
  public:
    unused_scheduler_callback(iengine_type* engine = NULL) :
      engine(engine), logoutput(false) { }
    
    virtual ~unused_scheduler_callback() {}

    
    void add_task(update_task_type task, double priority) {
      if (logoutput == false) {
        logger(LOG_WARNING, "This scheduler does not permit add task operations."
               "This warning will only be displayed once");
        logoutput = true;
      }
    }

    /** Creates a collection of tasks on all the vertices in
        'vertices', and all with the same update function and
        priority  */
    void add_tasks(const std::vector<vertex_id_type>& vertices, 
                   update_function_type func, double priority) {
      if (logoutput == false) {
        logger(LOG_WARNING, "This scheduler does not permit add task operations."
               "This warning will only be displayed once");
        logoutput = true;
      }
    }

    /**
     * Force the engine to abort.
     */
    void force_abort() {
      assert(engine != NULL);
      engine->stop();
    }


    void disable_buffering() {}

      
    /// Commits the tasks added
    void commit() { };
  private:
    iengine_type* engine;
    bool logoutput;
  };

}
#endif

