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



/**
 * This class defines a basic priority scheduler.
 **/
#ifndef GRAPHLAB_PRIORITY_SCHEDULER_HPP
#define GRAPHLAB_PRIORITY_SCHEDULER_HPP

#include <cmath>


#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>

#include <graphlab/engine/terminator/task_count_terminator.hpp>
#include <graphlab/options/options_map.hpp>



// #include <graphlab/logger/assertions.hpp>
// #include <graphlab/graph/graph.hpp>
// #include <graphlab/scope/iscope.hpp>
// 
// #include <graphlab/tasks/update_task.hpp>
// #include <graphlab/schedulers/ischeduler.hpp>
// #include <graphlab/parallel/pthread_tools.hpp>
// #include <graphlab/schedulers/support/direct_callback.hpp>
// #include <graphlab/schedulers/support/vertex_task_set.hpp>
// #include <graphlab/util/task_count_termination.hpp>

//#include <bitmagic/bm.h>

#include <graphlab/macros_def.hpp>
namespace graphlab {

   /** \ingroup group_schedulers
    */
  template<typename Engine>
  class priority_scheduler : public ischeduler<Engine> {

  public:

    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;

    typedef mutable_queue<vertex_id_type, double> priority_queue_type;


  private:

    /** The vertex functor set */
    vertex_functor_set<engine_type> vfun_set;

    /** The lock on the priority queue */
    spinlock queue_lock;     
    /** The queue over vertices */
    priority_queue_type pqueue;
      
    /** Used to assess termination */
    task_count_terminator term;
    
    
  public:
    
    priority_scheduler(const graph_type& graph, 
                       size_t ncpus,
                       const options_map& opts) :
      vfun_set(graph.num_vertices()) { }       

    void start() { term.reset(); };


    void schedule(vertex_id_type vid, 
                  const update_functor_type& fun) {   
      queue_lock.lock();
      const bool first_add = vfun_set.add(vid, fun);
      if(first_add) {
        term.new_job();
        pqueue.push(vid, fun.priority());
      } else {
        // Update the priority queue
        const std::pair<bool, double> pair = vfun_set.priority(vid);
        ASSERT_TRUE(pair.first);
        pqueue.update(vid, pair.second);
      }
      queue_lock.unlock();
    } // end of add_task
    
    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
        schedule(vid, fun);      
    } // end of schedule_all

    void completed(size_t cpuid,
                   vertex_id_type vid,
                   const update_functor_type& fun) {
      term.completed_job();
    }

    /** Get the next element in the queue */
    sched_status::status_enum get_next(size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {         
      queue_lock.lock();
      const bool success = !pqueue.empty();
      if(success) {
        ret_vid = pqueue.pop().first;
        const bool get_success = vfun_set.test_and_get(ret_vid, ret_fun);
        ASSERT_TRUE(get_success);                
      }
      queue_lock.unlock();
      if(success) {
        return sched_status::NEW_TASK;
      } else {
        return sched_status::EMPTY;
      }
    } // end of get_next_task

    iterminator& terminator() { return term; }

    static void print_options_help(std::ostream &out) { };

  }; // end of priority_queue class


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

