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

#ifndef GRAPHLAB_IMULTIPLEXLISTENER_HPP
#define GRAPHLAB_IMULTIPLEXLISTENER_HPP

#include <vector>


#include <graphlab/graph/graph.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  
    /**
      * \brief multiplexing monitor
      * 
      * This monitor class allows use of multiple monitor classes
      * simultaneously.
      */
  
  template<typename Graph>
  class monitor_multiplexer : 
    public imonitor<Graph> {
  public:
    typedef imonitor<Graph> imonitor_type;
    typedef typename imonitor_type::update_task_type update_task_type;
    typedef typename imonitor_type::iengine_type iengine_type;
    typedef typename imonitor_type::iscope_type iscope_type;

   
  private:
    std::vector< imonitor_type* > children;
    
  public:
        
    ~monitor_multiplexer() {
      foreach (imonitor_type* child, children) {
        if(child != NULL)
          delete child;
      }
    }

    /**
     * Add an additional monitor to the multiplexer.  The imonitor
     * instance must have a copy constructor.
     */
    template<typename Monitor>
    void add(const Monitor& child) {
      imonitor_type* child_ptr = new Monitor(child);
      assert(child_ptr != NULL);
      children.push_back(child_ptr);
    }
    
    /* Initialization, called by the engine */
    void init(iengine_type* engine) {
      foreach (imonitor_type * child, children) child->init(engine);
    }
    
    /* Engine calls */
    void engine_task_execute_start(update_task_type task, 
                                   iscope_type* scope, 
                                   size_t cpuid) { 
      foreach (imonitor_type * child, children) {
        child->engine_task_execute_start(task, scope, cpuid);
      }                                      
    }
    
    void engine_task_execute_finished(update_task_type task, 
                                      iscope_type* scope, 
                                      size_t cpuid) { 
      foreach (imonitor_type * child, children) {
        child->engine_task_execute_finished(task, scope, cpuid);
      }                                      
    }
    
    
    void engine_worker_starts(size_t cpuid) { 
      foreach (imonitor_type * child, children) {
        child->engine_worker_starts(cpuid);
      }                                      
    }
    
    
    void engine_worker_dies(size_t cpuid, int taskcount){ 
      foreach (imonitor_type * child, children) {
        child->engine_worker_dies(cpuid, taskcount);
      }                                      
    }
    
    
    
    
    /* Scheduler calls */
    void scheduler_task_added(update_task_type task, double priority){ 
      foreach (imonitor_type * child, children) {
        child->scheduler_task_added(task, priority);
      }                                      
    }
    
    
    void scheduler_task_promoted(update_task_type task, 
                                 double diffpriority, 
                                 double totalpriority){ 
      foreach (imonitor_type * child, children) {
        child->scheduler_task_promoted(task, diffpriority, totalpriority);
      }                                      
    }
    
    
    void scheduler_task_scheduled(update_task_type task, 
                                  double current_max_priority){ 
      foreach (imonitor_type * child, children) {
        child->scheduler_task_scheduled(task, current_max_priority);
      }                                      
    }
    
    
    void scheduler_task_pruned(update_task_type task) { 
      foreach (imonitor_type * child, children) {
        child->scheduler_task_pruned(task);
      }                                      
    }
    
    /* Application calls */
    void app_set_vertex_value(vertex_id_t vid, double value) { 
      foreach (imonitor_type * child, children) {
        child->app_set_vertex_value(vid, value);
      }                                      
    }
    
    
    /* Called by application to help visualizers to scale values properly */
    void app_set_vertex_value_scale(double min, double max) { 
      foreach (imonitor_type * child, children) {
        child->app_set_vertex_value_scale(min, max);
      }                                      
    }
    


    
  };
}


#include <graphlab/macros_undef.hpp>

#endif

