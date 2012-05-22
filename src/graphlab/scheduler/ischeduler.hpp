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
#include <graphlab/options/options_map.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>


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
   * \ingroup group_schedulers This describes the interface/concept
   * for the scheduler. The engine will be passed the scheduler type
   * as a template argument, so the scheduler must inherit and satisfy
   * this interface EXACTLY. Note that all functions (with the
   * exception of the constructor and destructor) must be thread-safe.
   */
  template<typename Graph, typename MessageType>
  class ischeduler {
  public:

    typedef Graph graph_type;
    typedef MessageType message_type;
    typedef typename graph_type::vertex_id_type vertex_id_type;

    
    /// destructor
    virtual ~ischeduler() {};
        
    /** Called by engine before starting the schedule.
     *  This function will only be called once throughout the lifetime
     * of the scheduler.
     */
    virtual void start() = 0;


    /**
     * Adds a message destined to the vertex vid to the schedule.
     */
    virtual void schedule(const vertex_id_type vid, 
                          const message_type& message) = 0;

    /**
     * Adds a message destined to the vertex vid to the schedule
     */
    virtual void schedule_from_execution_thread(const size_t cpuid, 
                                                const vertex_id_type vid, 
                                                const message_type& message) {
      schedule(vid, message);
    }
    


    /** 
     * Schedule the message to be received by all vertices in the
     * graph.
     */
    virtual void schedule_all(const message_type& message,
                              const std::string& order = "sequential") = 0;


    /**
     * This function is called by the engine to ask for the next
     * message to process.  The message and receiving vertex are
     * returned in ret_msg and ret_vid respectively.
     *
     *  \retval NEWTASK There is a new message to process
     * 
     *  \retval EMPTY There are no messages to process
     */
    virtual sched_status::status_enum 
    get_next(const size_t cpuid, vertex_id_type& ret_vid,
             message_type& ret_msg) = 0;


    /**
     * Get a message for a specific vertex.
     */
    virtual sched_status::status_enum
    get_specific(vertex_id_type vid, message_type& ret_msg) {
      return sched_status::EMPTY;
    }


    /**
     * Inserts a message for vertex vid to be maintained, but do not
     * update the schedule. Must be implemented for the distributed
     * engines.
     *
     * TODO: do we need this functionality
     */
    virtual void place(vertex_id_type vid, const message_type& msg) {  }


    /**
     * Schedules vertex vid using the stored message that was
     * previously placed using place. Must be implemented for the
     * distributed engines.
     *
     * TODO: do we need this functionality
     */
    virtual void
    schedule_from_execution_thread(const size_t cpuid, vertex_id_type vid) {
      schedule(vid);
    }
    
    /**
     * Schedules vertex vid using the stored message that was previously
     * placed using place. Must be implemented for the distributed engines.
     *
     * TODO: do we need this functionality
     */
    virtual void schedule(vertex_id_type vid) { }
    
    /**
     * This is called after a message has been received.
     */
    virtual void completed(const size_t cpuid,
                           const vertex_id_type vid,
                           const message_type& message) { }


    /**
     * Get the terminator associated with this scheduler
     */
    virtual iterminator& terminator() = 0;
    
    /**
     * Optional to implement. Count the number of message combination
     * operations performed. Returns (size_t)(-1) if not available.
     */
    virtual size_t num_joins() const { return (size_t)(-1);}

    /**
     * Print a help string describing the options that this scheduler
     * accepts.
     */
    static void print_options_help(std::ostream& out) { };

  };

}
#endif

