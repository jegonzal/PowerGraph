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

