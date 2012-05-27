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


#ifndef GRAPHLAB_QUEUED_FIFO_SCHEDULER_HPP
#define GRAPHLAB_QUEUED_FIFO_SCHEDULER_HPP

#include <algorithm>
#include <queue>


#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/parallel/atomic_add_vector.hpp>

#include <graphlab/scheduler/get_message_priority.hpp>
#include <graphlab/options/graphlab_options.hpp>


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
  template<typename Message>
  class queued_fifo_scheduler : public ischeduler<Message> {
  
  public:

    typedef Message message_type;

    typedef std::deque<vertex_id_type> queue_type;

  private:
    atomic_add_vector<message_type> messages;
    std::deque<queue_type> master_queue;
    mutex master_lock;
    size_t sub_queue_size;
    std::vector<queue_type> in_queues;
    std::vector<mutex> in_queue_locks;
    std::vector<queue_type> out_queues;
    double min_priority;
    // Terminator
  public:

    queued_fifo_scheduler(size_t num_vertices,
                          const graphlab_options& opts) :
      messages(num_vertices), 
      sub_queue_size(100), 
      in_queues(opts.get_ncpus()), in_queue_locks(opts.get_ncpus()), 
      out_queues(opts.get_ncpus()),
      min_priority(-std::numeric_limits<double>::max()){ 
      set_options(opts);
    }

    void set_options(const graphlab_options& opts) {
      size_t new_ncpus = opts.get_ncpus();
      // check if ncpus changed
      if (new_ncpus != in_queues.size()) {
        logstream(LOG_INFO) << "Changing ncpus from " << in_queues.size()
                            << " to " << new_ncpus << std::endl;
        ASSERT_GE(new_ncpus, 1);
        // if increasing ncpus, we just resize the queues
        // push everything in in_queues to the master queue
        for (size_t i = 0;i < in_queues.size(); ++i) {
          master_queue.push_back(in_queues[i]);
          in_queues[i].clear();
        }
        // resize the queues
        in_queues.resize(new_ncpus);
        in_queue_locks.resize(new_ncpus);
        out_queues.resize(new_ncpus);
      }
      
      // read the remaining options.
      std::vector<std::string> keys = opts.get_engine_args().get_option_keys();
      foreach(std::string opt, keys) {
        if (opt == "queuesize") {
          opts.get_engine_args().get_option("queuesize", sub_queue_size);
        } else if (opt == "min_priority") {
          opts.get_engine_args().get_option("min_priority", min_priority);
        } else {
          logstream(LOG_ERROR) << "Unexpected Scheduler Option: " << opt << std::endl;
        }
      }
    }
    
    void start() { 
      master_lock.lock();
      for (size_t i = 0;i < in_queues.size(); ++i) {
        master_queue.push_back(in_queues[i]);
        in_queues[i].clear();
      }
      master_lock.unlock();
    }

    void schedule(const vertex_id_type vid, 
                  const message_type& msg) {
      // If this is a new message, schedule it
      // the min priority will be taken care of by the get_next function
      if (messages.add(vid, msg)) {
        const size_t cpuid = random::rand() % in_queues.size();
        in_queue_locks[cpuid].lock();
        queue_type& queue = in_queues[cpuid];
        queue.push_back(vid);
        if(queue.size() > sub_queue_size) {
          master_lock.lock();
          queue_type emptyq;
          master_queue.push_back(emptyq);
          master_queue.back().swap(queue);
          master_lock.unlock();
        }
        in_queue_locks[cpuid].unlock();
      } 
    } // end of schedule

    void schedule_from_execution_thread(const size_t cpuid,
                                        const vertex_id_type vid) {      
      if (!messages.empty(vid)) {
        ASSERT_LT(cpuid, in_queues.size());
        in_queue_locks[cpuid].lock();
        queue_type& queue = in_queues[cpuid];
        queue.push_back(vid);
        if(queue.size() > sub_queue_size) {
          master_lock.lock();
          queue_type emptyq;
          master_queue.push_back(emptyq);
          master_queue.back().swap(queue);
          master_lock.unlock();
        }
        in_queue_locks[cpuid].unlock();
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
                   const message_type& msg) {
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

    
    /** Get the next element in the queue */
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       message_type& ret_msg) {
      // if the local queue is empty try to get a queue from the master
      while(1) {
        if(out_queues[cpuid].empty()) {
          master_lock.lock();
          if(!master_queue.empty()) {
            out_queues[cpuid].swap(master_queue.front());
            master_queue.pop_front();
          }
          master_lock.unlock();
        }
        // if the local queue is still empty see if there is any local
        // work left
        in_queue_locks[cpuid].lock();
        if(out_queues[cpuid].empty() && !in_queues[cpuid].empty()) {
          out_queues[cpuid].swap(in_queues[cpuid]);
        }
        in_queue_locks[cpuid].unlock();
        // end of get next
        queue_type& queue = out_queues[cpuid];
        if(!queue.empty()) {
          ret_vid = queue.front();
          queue.pop_front();
          if(messages.test_and_get(ret_vid, ret_msg)) {
            if (scheduler_impl::get_message_priority(ret_msg) >= min_priority) {
              return sched_status::NEW_TASK;
            } else {
                // it is below priority. try to put it back. If putting it back
                // makes it exceed priority, reschedule it
              message_type combined_message;
              messages.add(ret_vid, ret_msg, combined_message);
              double ret_priority = scheduler_impl::get_message_priority(combined_message);
              if(ret_priority >= min_priority) {
                // aargh. we put it back and it exceeded priority
                // stick it back in the queue.
                in_queue_locks[cpuid].lock();
                in_queues[cpuid].push_back(ret_vid);
                in_queue_locks[cpuid].unlock();
              }
            }
          }
        } else {
          return sched_status::EMPTY;
        }
      }
    } // end of get_next_task


    size_t num_joins() const {
      return messages.num_joins();
    }
    /**
     * Print a help string describing the options that this scheduler
     * accepts.
     */
    static void print_options_help(std::ostream& out) { 
      out << "\t queuesize: [the size at which a subqueue is "
          << "placed in the master queue. default = 100]\n"
          << "min_priority = [double, minimum priority required to receive \n"
          << "\t a message, default = -inf]\n";
    }


  }; 


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

