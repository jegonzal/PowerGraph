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
 * \author akyrola
 * This class defines a multiqueue FIFO scheduler, i.e each thread manages
 * one or more FIFO queues. Tasks are added to queues with an effort to balance
 * load. 
 **/
#ifndef GRAPHLAB_MULTIQUEUE_FIFO_SCHEDULER_HPP
#define GRAPHLAB_MULTIQUEUE_FIFO_SCHEDULER_HPP

#include <algorithm>
#include <queue>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>

#include <graphlab/scheduler/terminator/task_count_terminator.hpp>
#include <graphlab/options/options_map.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \ingroup group_schedulers 
   *
   * This class defines a multiple queue approximate fifo scheduler.
   * Each processor has its own in_queue which it puts new tasks in
   * and out_queue which it pulls tasks from.  Once a processors
   * in_queue gets too large, the entire queue is placed at the end of
   * the shared master queue.  Once a processors out queue is empty it
   * grabs the next out_queue from the master.
   */
  template<typename Engine>
  class multiqueue_fifo_scheduler : public ischeduler<Engine> {
  
  public:

    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;


    typedef std::deque<vertex_id_type> queue_type;
    

  private:

    vertex_functor_set<engine_type> vfun_set;
    std::deque<queue_type> master_queue;
    mutex master_lock;
    size_t sub_queue_size;
    std::vector<queue_type> in_queues;
    std::vector<queue_type> out_queues;

    // Terminator
    shared_termination term;
 


  public:

    multiqueue_fifo_scheduler(const graph_type& graph, 
                              size_t ncpus,
                              const options_map& opts) :
      vfun_set(graph.num_vertices()), 
      sub_queue_size(1000), 
      in_queues(ncpus), out_queues(ncpus), term(ncpus) {  }

    void start() { term.reset(); }
   

    void schedule(const size_t cpuid,
                  const vertex_id_type vid, 
                  const update_functor_type& fun) {      
      if (vfun_set.add(vid, fun)) {
        term.new_job(cpuid);
        queue_type& queue = in_queues[cpuid];
        queue.push_back(vid);
        if(queue.size() > sub_queue_size) {
          master_lock.lock();
          master_queue.push_back(queue);
          master_lock.unlock();
          queue.clear();
        }
      } 
    } // end of schedule

    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
        schedule(vid % in_queues.size(), vid, fun);      
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
      // if the local queue is empty try to get a queue from the master
      if(out_queues[cpuid].empty()) {
        master_lock.lock();
        if(!master_queue.empty()) {
          out_queues[cpuid] = master_queue.front();
          master_queue.pop_front();
        }
        master_lock.unlock();
      }
      // if the local queue is still empty see if there is any local
      // work left
      if(out_queues[cpuid].empty() && !in_queues[cpuid].empty()) {
        out_queues[cpuid].swap(in_queues[cpuid]);
      }
      // end of get next
      queue_type& queue = out_queues[cpuid];
      if(!queue.empty()) {
        ret_vid = queue.front();
        queue.pop_front();
        const bool get_success = vfun_set.test_and_get(ret_vid, ret_fun);
        ASSERT_TRUE(get_success);
        return sched_status::NEW_TASK;
      } else {
        return sched_status::EMPTY;
      }
    } // end of get_next_task

    iterminator& terminator() { return term; }



  }; 


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

