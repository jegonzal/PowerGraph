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

