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
 * \author jegonzal This class defines a multiqueue version of the
 * priority scheduler.
 **/
#ifndef GRAPHLAB_PRIORITY_SCHEDULER_HPP
#define GRAPHLAB_PRIORITY_SCHEDULER_HPP

#include <queue>
#include <cmath>
#include <cassert>


#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_map.hpp>

#include <graphlab/scheduler/terminator/critical_termination.hpp>
#include <graphlab/options/options_map.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  /** \ingroup group_schedulers
   */
  template<typename Message>
  class priority_scheduler : public ischeduler<Message> {  
  public:

    typedef Message message_type;

    typedef mutable_queue<vertex_id_type, double> queue_type;

  private:
    vertex_map<message_type> messages;
    std::vector<queue_type> queues;
    std::vector<spinlock>   locks;
    size_t multi;
    std::vector<size_t>     current_queue;

    /** Max priority */
    double min_priority;

    // Terminator
    critical_termination term;
 


  public:

    priority_scheduler(size_t num_vertices,
                       size_t ncpus,
                       const options_map& opts) :
      messages(num_vertices), multi(0),
      current_queue(ncpus), 
      min_priority(-std::numeric_limits<double>::max()),
      term(ncpus) {     
      const bool is_set = opts.get_option("min_priority", min_priority);
      if(is_set) {
        logstream(LOG_INFO) << "The minimum scheduling priority was set to " 
                            << min_priority << std::endl;
      }
      opts.get_option("multi", multi);
      const size_t nqueues = std::max(multi*ncpus, size_t(1));
      if(multi > 0) {
        logstream(LOG_INFO) << "Using " << multi 
                            << " queues per thread." << std::endl;
      }
      queues.resize(nqueues);
      locks.resize(nqueues);
    }

    void start() { term.reset(); }
   

    void schedule(const vertex_id_type vid, 
                  const message_type& msg) {      
      const size_t idx = vid % queues.size();
      double priority = 0;
      locks[idx].lock(); 
      if (messages.add(vid, msg, priority)) {
        queues[idx].push(vid, priority); 
      } else { queues[idx].update(vid, priority); }
      locks[idx].unlock();
      // signal thread to weak up if they have not already
      term.new_job();     
    } // end of schedule

    void schedule_all(const message_type& msg,
                      const std::string& order) {
      if(order == "shuffle") {
        std::vector<vertex_id_type> permutation = 
          random::permutation<vertex_id_type>(messages.size());       
        foreach(vertex_id_type vid, permutation)  schedule(vid, msg);
      } else {
        for (vertex_id_type vid = 0; vid < messages.size(); ++vid)
          schedule(vid, msg);      
      }
    } // end of schedule_all

    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const message_type& msg) { term.completed_job(); }


    /** Get the next element in the queue */
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       message_type& ret_msg) {
      /* Check all of my queues for a task */
      for(size_t i = 0; i < multi; ++i) {
        const size_t idx = (++current_queue[cpuid] % multi) + cpuid * multi;
        locks[idx].lock();
        if(!queues[idx].empty() && 
           queues[idx].top().second >= min_priority) {
          ret_vid = queues[idx].pop().first;
          const bool get_success = messages.test_and_get(ret_vid, ret_msg);
          locks[idx].unlock();
          ASSERT_TRUE(get_success);
          return sched_status::NEW_TASK;          
        }
        locks[idx].unlock();
      }
      /* Check all the queues */
      for(size_t i = 0; i < queues.size(); ++i) {
        const size_t idx = ++current_queue[cpuid] % queues.size();
        if(!queues[idx].empty()) { // quick pretest
          locks[idx].lock();
          if(!queues[idx].empty() && 
             queues[idx].top().second >= min_priority) {
            ret_vid = queues[idx].pop().first;
            const bool get_success = messages.test_and_get(ret_vid, ret_msg);
            locks[idx].unlock();
            ASSERT_TRUE(get_success);
            return sched_status::NEW_TASK;          
          }
          locks[idx].unlock();
        }
      }
      return sched_status::EMPTY;     
    } // end of get_next_task

    iterminator& terminator() { return term; }

    size_t num_joins() const {
      return messages.num_joins();
    }


    static void print_options_help(std::ostream& out) { 
      out << "\t mult=1: number of queues per thread.\n" 
          << "\t min_priority=-infty Minimum priority required "
          << "\t    to receive the message." << std::endl;
    }

  }; // end of class priority scheduler


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

