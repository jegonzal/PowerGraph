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




#ifndef DIRECT_CALLBACK_HPP
#define DIRECT_CALLBACK_HPP

#include <vector>
#include <algorithm>

#include <graphlab/graph/graph.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/schedulers/icallback.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {

  /**
     This is a callback class that is passed to the update functions, and allow the
     update functions to add tasks back into the scheduler
     \see ischeduler_callback */
  template<typename Graph>
  class direct_callback : 
    public icallback<Graph> {
  public:
    typedef Graph graph_type;
    typedef icallback<Graph> base;
    typedef ischeduler<Graph> scheduler_type;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
   
  protected: 
  

    typedef std::pair<update_task_type, double> task_priority_pair;
    
    scheduler_type* scheduler;    /// a pointer to the owning scheduler
    iengine_type* engine; 
    bool buffering_enabled;
    std::vector< task_priority_pair > tasks; /// The collection of tasks inserted
  public:
  
    direct_callback(scheduler_type* scheduler = NULL,
                    iengine_type* engine = NULL) :
      scheduler(scheduler), engine(engine),
      buffering_enabled(false) { }
  
    ~direct_callback() { }


    void set_engine(iengine_type* eng) {
      engine = eng;
    }
  
    void enable_buffering() {
      buffering_enabled = true;
    }

  
    void add_task(update_task_type task, double priority) {
      assert(task.function() != NULL);
      if (!buffering_enabled) {
      	scheduler->add_task(task, priority);  
      } else {
      	tasks.push_back(task_priority_pair(task, priority));
      }
    }
  
    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority) {
      foreach(vertex_id_t vertex, vertices) {
        add_task(update_task_type(vertex, func), priority);
      }    
    }
    
    
    int num_of_buffered_tasks() {
      return (int)tasks.size();
    }
    
    
    
    void commit() {
      if(buffering_enabled) {
        foreach(task_priority_pair tp, tasks) {
          scheduler->add_task(tp.first, tp.second);
        }
        tasks.clear();
      }
    }

    void force_abort() {
      assert(engine != NULL);
      engine->stop();
    }
    
    
  };
}

#include <graphlab/macros_undef.hpp>


#endif

