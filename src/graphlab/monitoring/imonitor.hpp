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





//
// Generic listener interface for monitoring, debugging and profiling.
//

#ifndef GRAPHLAB_IMONITOR_HPP
#define GRAPHLAB_IMONITOR_HPP


#include <graphlab/graph/graph.hpp>


namespace graphlab {
  
  
  template<typename Engine>
  class imonitor {    
  public:
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
    typedef typename engine_type::update_functor_type update_functor_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;


    /* Initialization, called by the engine */
    virtual void init(engine_type* engine) { }
    
    /* Engine calls */
    virtual void engine_task_execute_start(vertex_id_type vid, 
                                           const update_functor_type& ufun,
                                           size_t cpuid) { }
    
    virtual void engine_task_execute_finished(vertex_id_type vid,
                                              const update_functor_type& ufun,
                                              size_t cpuid) { }
    
    virtual void engine_worker_starts(size_t cpuid) { }
    
    virtual void engine_worker_dies(size_t cpuid, int taskcount) { }
    
    
    
    /* Scheduler calls */
    virtual void scheduler_task_added(vertex_id_type vid,
                                      const update_functor_type& fun) { }
    
    virtual void scheduler_task_promoted(vertex_id_type vid,
                                         const update_functor_type& fun,
                                         double diffpriority, 
                                         double totalpriority) { }
    
    virtual void scheduler_task_scheduled(vertex_id_type vid,
                                          const update_functor_type& fun) { }
    
    
    /* Application calls */
    // virtual void app_set_vertex_value(vertex_id_type vid, double value) { }
    
    /* Called by application to help visualizers to scale values properly */
    // virtual void app_set_vertex_value_scale(double min, double max) { }
    
  };
}

#endif

