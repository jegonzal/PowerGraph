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


#ifndef GRAPHLAB_CONSOLE_MONITOR_HPP
#define GRAPHLAB_CONSOLE_MONITOR_HPP

#include <iostream>
#include <graphlab/monitoring/imonitor.hpp>

namespace graphlab {

  template<typename Graph>
  class console_monitor : 
    public imonitor<Graph> {
    
  public:
    
    typedef typename imonitor<Graph>::update_task_type update_task_type;
    typedef typename imonitor<Graph>::iengine_type iengine_type;
    
  public:
    unsigned long output_priority_freq_tasks;
    unsigned long task_counter;
    double current_max_priority;
    
    console_monitor(unsigned long output_priority_freq_tasks = 1000) :
      output_priority_freq_tasks(output_priority_freq_tasks),
      task_counter(0),
      current_max_priority(0) { }
    
    virtual void init(iengine_type* engine) {
      std::cout << "====== CONSOLE LISTENER STARTED =====" << std::endl;
    }


    void scheduler_task_added(update_task_type task, double priority) {
      current_max_priority = std::max(current_max_priority, priority);
      // Not thread safe!!!
      if (++task_counter % output_priority_freq_tasks == 0) {
        std::cout << "Approx. task count "
                  << task_counter
                  << " with max priority "
                  << current_max_priority
                  << " adding task " << size_t(task.function()) 
                  << std::endl;
        current_max_priority = 0;
      }
    }

    


    
    void set_status_str(char * str) {
      printf("Status %s\n", str); 
    }
  };
  
  
}

#endif

