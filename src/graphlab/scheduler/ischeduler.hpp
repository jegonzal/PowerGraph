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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#ifndef GRAPHLAB_ISCHEDULER_HPP
#define GRAPHLAB_ISCHEDULER_HPP

#include <vector>
#include <sstream>
#include <ostream>
#include <graphlab/monitoring/imonitor.hpp>
#include <graphlab/options/options_map.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>

#include <graphlab/metrics/metrics.hpp>

namespace graphlab {
 
  /**
   * This is an enumeration for the possible return values for
   * get_next_tasks
   */
  struct sched_status {
    /// \brief the possible scheduler status.
    enum status_enum {
      NEW_TASK,      /**< The get_next_tasks function returned a new task 
                        to be executed */
      EMPTY,         /**< The schedule is empty. */      
    };
  };


  /**
   * \ingroup group_schedulers
   * This describes the interface/concept for the scheduler. The
   * engine will be passed the scheduler type as a template argument,
   * so the scheduler must inherit and satisfy this interface
   * EXACTLY. Note that all functions (with the exception of the
   * constructor and destructor) must be thread-safe.
   */
  template<typename Engine>
  class ischeduler {
  public:

    typedef Engine engine_type;
    typedef typename engine_type::graph_type           graph_type;
    typedef typename engine_type::update_functor_type  update_functor_type;

    typedef typename graph_type::vertex_id_type    vertex_id_type;
    typedef typename graph_type::edge_id_type      edge_id_type;
    typedef typename graph_type::vertex_color_type vertex_color_type;
    
    /// destructor
    virtual ~ischeduler() {};
        
    /** Called by engine before starting the schedule.
     *  This function will only be called once throughout the lifetime
     * of the scheduler.
     */
    virtual void start() = 0;


    /**
     * Adds an update task with a particular priority. 
     * This function may be called at anytime.
     */
    virtual void schedule(const vertex_id_type vid, 
                          const update_functor_type& fun) = 0;

    /**
     * Adds an update task with a particular priority. 
     * This function should be called only by the same threads consuming tasks
     */
    virtual void schedule_from_execution_thread(const vertex_id_type vid, 
                                          const update_functor_type& fun) {
      schedule(vid, fun);
    }
    
    
    /** 
     * Creates a collection of tasks on all the vertices in the graph,
     * with the same update function and priority
     * This function may be called at anytime.
     */
    virtual void schedule_all(const update_functor_type& fun) = 0;


    /**
     * This function is called by the engine to ask for new work to
     * do.  The update task to be executed is returned in ret_task.
     *
     *  \retval NEWTASK There is an update task in ret_task to be
     *   executed
     * 
     *  \retval EMPTY There are no tasks available in the scheduler.
     */
    virtual sched_status::status_enum 
    get_next(const size_t cpuid,
             vertex_id_type& ret_vid,
             update_functor_type& ret_fun) = 0;



    /**
     * This is called after a task has been executed
     */
    virtual void completed(const size_t cpuid,
                           const vertex_id_type vid,
                           const update_functor_type& fun) { }


    /**
     * Get the terminator associated with this scheduler
     */
    virtual iterminator& terminator() = 0;
    


    /**
     * Print a help string describing the options that this scheduler
     * accepts.
     */
    static void print_options_help(std::ostream& out) { };

  };

}
#endif

