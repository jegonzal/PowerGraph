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
 * This class defines a basic FIFO (First In First Out) scheduler
 **/
#ifndef GRAPHLAB_FIFO_SCHEDULER_HPP
#define GRAPHLAB_FIFO_SCHEDULER_HPP

#include <cmath>
#include <cassert>

#include <queue>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>

#include <graphlab/scheduler/terminator/count_termination.hpp>
#include <graphlab/options/options_map.hpp>





#include <graphlab/macros_def.hpp>

namespace graphlab {


  /** \ingroup group_schedulers
   */
  template<typename Engine>
  class fifo_scheduler : public ischeduler<Engine> {
  public:
    

    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;

  private:

    //! The set of functions to apply
    vertex_functor_set<engine_type> vfun_set;
    std::queue<vertex_id_type> queue; 
    spinlock queue_lock; 
    count_termination term;



  public:

    fifo_scheduler(const graph_type& graph, 
                   size_t ncpus,
                   const options_map& opts) :
      vfun_set(graph.num_vertices()) { }
    

    //! Reset the terminator
    void start() { term.reset(); }

    void schedule(const size_t cpuid,
                  const vertex_id_type vid, 
                  const update_functor_type& fun) {      
      if (vfun_set.add(vid, fun)) {
        term.new_job();
        queue_lock.lock();
        queue.push(vid);
        queue_lock.unlock();
      } 
    } // end of schedule

    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
        schedule(0, vid, fun);      
    } // end of schedule_all

    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const update_functor_type& fun) {
      term.completed_job();
    }


    /** Get the next element in the queue */
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {         
      queue_lock.lock();
      const bool success = !queue.empty();
      if(success) {
        ret_vid = queue.front();
        queue.pop();
      }
      queue_lock.unlock();
      if(success) {
        const bool get_success = vfun_set.test_and_get(ret_vid, ret_fun);
        ASSERT_TRUE(get_success);
        return sched_status::NEW_TASK;
      } else {
        return sched_status::EMPTY;
      }
    } // end of get_next_task

    iterminator& terminator() { return term; }


  }; // end of fifo_scheduler


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

