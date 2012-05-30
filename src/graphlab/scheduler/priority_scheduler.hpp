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
#include <graphlab/parallel/atomic_add_vector2.hpp>
#include <graphlab/scheduler/get_message_priority.hpp>
#include <graphlab/options/graphlab_options.hpp>


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
    atomic_add_vector2<message_type> messages;
    std::vector<queue_type> queues;
    std::vector<spinlock>   locks;
    size_t multi;
    std::vector<size_t>     current_queue;

    double min_priority;

    // Terminator
 


  public:

    priority_scheduler(size_t num_vertices,
                       const graphlab_options& opts) :
      messages(num_vertices), multi(3),
      current_queue(opts.get_ncpus()), 
      min_priority(-std::numeric_limits<double>::max()) {     
      set_options(opts);
    }

    void set_options(const graphlab_options& opts) {
      size_t new_ncpus = opts.get_ncpus();
      // check if ncpus changed
      if (new_ncpus != current_queue.size()) {
        logstream(LOG_INFO) << "Changing ncpus from " << current_queue.size()
                            << " to " << new_ncpus << std::endl;
        ASSERT_GE(new_ncpus, 1);
        current_queue.resize(new_ncpus);
      }
      
      std::vector<std::string> keys = opts.get_scheduler_args().get_option_keys();
      foreach(std::string opt, keys) {
        if (opt == "multi") {
          opts.get_scheduler_args().get_option("multi", multi);
        } else if (opt == "min_priority") {
          opts.get_scheduler_args().get_option("min_priority", min_priority);
        } else {
          logstream(LOG_ERROR) << "Unexpected Scheduler Option: " << opt << std::endl;
        }
      }
      
      const size_t nqueues = std::max(multi*current_queue.size(), size_t(1));
      // changing the number of queues.
      // reinsert everything
      if (nqueues != queues.size()) {
        std::vector<queue_type> old_queues;
        std::swap(old_queues, queues);
        queues.resize(nqueues);
        locks.resize(nqueues);
        
        size_t idx = 0;
        for (size_t i = 0;i < old_queues.size(); ++i) {
          while (!old_queues[i].empty()) {
            queues[idx].push(old_queues[i].top().first, 
                            old_queues[i].top().second);
            old_queues[i].pop();
            ++idx;
          }
        }
      }
    }

    void start() {  }
   

    void schedule(const vertex_id_type vid, 
                  const message_type& msg) {      
      const size_t idx = vid % queues.size();
      message_type combined_message;
      
      messages.add(vid, msg, combined_message);
      // If the new priority will is above priority, put it in the queue
      if (scheduler_impl::get_message_priority(combined_message) >= min_priority) {
        locks[idx].lock(); 
        queues[idx].push_or_update(vid, 
                                    scheduler_impl::get_message_priority(combined_message)); 
        locks[idx].unlock();
      }
    } // end of schedule



    void schedule_from_execution_thread(const size_t cpuid,
                                        const vertex_id_type vid) {
      const size_t idx = vid % queues.size();
      message_type combined_message;
      
      messages.peek(vid, combined_message);
      // If the new priority will is above priority, put it in the queue
      if (scheduler_impl::get_message_priority(combined_message) >= min_priority) {
        locks[idx].lock(); 
        queues[idx].push_or_update(vid, 
                                    scheduler_impl::get_message_priority(combined_message)); 
        locks[idx].unlock();
      }
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
                   const message_type& msg) { }


    /** Get the next element in the queue */
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       message_type& ret_msg) {
      while(1) {
        /* Check all of my queues for a task */
        for(size_t i = 0; i < multi; ++i) {
          const size_t idx = (++current_queue[cpuid] % multi) + cpuid * multi;
          locks[idx].lock();
          if(!queues[idx].empty() && 
            queues[idx].top().second >= min_priority) {
            ret_vid = queues[idx].pop().first;
            const bool get_success = messages.test_and_get(ret_vid, ret_msg);
            locks[idx].unlock();
            if(get_success) return sched_status::NEW_TASK;
            else continue;
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
              if(get_success) return sched_status::NEW_TASK;
              else continue;
            }
            locks[idx].unlock();
          }
        }
        break;
      }
      return sched_status::EMPTY;     
    } // end of get_next_task


    size_t num_joins() const {
      return messages.num_joins();
    }



    sched_status::status_enum 
    get_specific(vertex_id_type vid,
                 message_type& ret_msg) {
      bool get_success = messages.test_and_get(vid, ret_msg); 
      if (get_success) return sched_status::NEW_TASK;
      else return sched_status::EMPTY;
    }

    void place(vertex_id_type vid,
                 const message_type& msg) {
      messages.add(vid, msg);
    }


    static void print_options_help(std::ostream& out) { 
      out << "\t multi = [number of queues per thread. Default = 3].\n" 
          << "min_priority = [double, minimum priority required to receive \n"
          << "\t a message, default = -inf]\n";
    }

  }; // end of class priority scheduler


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

