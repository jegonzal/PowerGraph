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
 * This class defines a basic FIFO (First In First Out) scheduler
 **/
#ifndef GRAPHLAB_FIFO_SCHEDULER_HPP
#define GRAPHLAB_FIFO_SCHEDULER_HPP

#include <cmath>
#include <cassert>

#include <queue>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/engine/terminator/iterminator.hpp>
#include <graphlab/options/options_map.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>




#include <graphlab/macros_def.hpp>

namespace graphlab {


  /** \ingroup group_schedulers
   */
  template<typename Engine>
  class fifo_scheduler: public ischeduler<Engine> {
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
    mutex queue_lock; 
    iterminator& terminator;



  public:

    fifo_scheduler(const graph_type& graph, 
                   iterminator& terminator, 
                   size_t ncpus,
                   const options_map& opts) :
      vertex_tasks(graph.num_vertices()), 
      terminator(terminator) {  }
    

    void start() { }

    void schedule(vertex_id_type vid, 
                  const update_functor_type& fun) {      
      if (vfun_set.add(task)) {
        terminator.new_job();
        queue_lock().lock();
        queue.push(vid);
        queue_lock.unlock();
      } 
    } // end of add_task

    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
        schedule(vid, fun);      
    } // end of add_task_to_all

    void completed(size_t cpuid,
                   vertex_id_type vid,
                   const update_functor_type& fun) {
      terminator.completed_job();
    }


    /** Get the next element in the queue */
    sched_status::status_enum get_next(size_t cpuid,
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
        const bool get_success = vfun_set.test_and_get(ret_fun);
        ASSERT_TRUE(get_success);
        return sched_status::NEWTASK;
      } else {
        return sched_status::EMPTY;
      }
    } // end of get_next_task

    void set_options(const scheduler_options &opts) { }


  }; // end of fifo_scheduler


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

