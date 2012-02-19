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


#ifndef GRAPHLAB_FIFO_SCHEDULER_HPP
#define GRAPHLAB_FIFO_SCHEDULER_HPP

#include <algorithm>
#include <queue>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>

#include <graphlab/scheduler/terminator/critical_termination.hpp>
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
  class fifo_scheduler : public ischeduler<Engine> {
  
  public:

    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;

    typedef std::deque<vertex_id_type> queue_type;

  private:

    vertex_functor_set<update_functor_type> vfun_set;
    std::vector<queue_type> queues;
    std::vector<spinlock>   locks;
    size_t multi;
    std::vector<size_t>     current_queue;

    // Terminator
    critical_termination term;
 


  public:

    fifo_scheduler(const graph_type& graph, 
                   size_t ncpus,
                   const options_map& opts) :
      vfun_set(graph.num_vertices()), multi(0),
      current_queue(ncpus), term(ncpus) {     
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
                  const update_functor_type& fun) {      
      if (vfun_set.add(vid, fun)) {
        /* "Randomize" the task queue task is put in. Note that we do
           not care if this counter is corrupted in race conditions
           Find first queue that is not locked and put task there (or
           after iteration limit) Choose two random queues and use the
           one which has smaller size */
        // M.D. Mitzenmacher The Power of Two Choices in Randomized
        // Load Balancing (1991)
        // http://www.eecs.harvard.edu/~michaelm/postscripts/mythesis.
        size_t idx = 0;
        if(queues.size() > 1) {
          const size_t prod = 
            random::fast_uniform(size_t(0), queues.size() * queues.size() - 1);
          const size_t r1 = prod / queues.size();
          const size_t r2 = prod % queues.size();
          idx = (queues[r1].size() < queues[r2].size()) ? r1 : r2;  
        }
        locks[idx].lock(); queues[idx].push_back(vid); locks[idx].unlock();
        term.new_job(idx);
      }
    } // end of schedule

    void schedule_all(const update_functor_type& fun) {
      for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid) {
        if(vfun_set.add(vid,fun)) {
          term.new_job();
          const size_t idx = vid % queues.size();
          locks[idx].lock(); queues[idx].push_back(vid); locks[idx].unlock();
        }
      }
    } // end of schedule_all

    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const update_functor_type& fun) { term.completed_job(); }


    /** Get the next element in the queue */
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {
      /* Check all of my queues for a task */
      for(size_t i = 0; i < multi; ++i) {
        const size_t idx = (++current_queue[cpuid] % multi) + cpuid * multi;
        locks[idx].lock();
        if(!queues[idx].empty()) {
          ret_vid = queues[idx].front();
          queues[idx].pop_front();
          const bool get_success = vfun_set.test_and_get(ret_vid, ret_fun);
          ASSERT_TRUE(get_success);
          locks[idx].unlock();
          return sched_status::NEW_TASK;          
        }
        locks[idx].unlock();
      }
      /* Check all the queues */
      for(size_t i = 0; i < queues.size(); ++i) {
        const size_t idx = ++current_queue[cpuid] % queues.size();
        if(!queues[idx].empty()) { // quick pretest
          locks[idx].lock();
          if(!queues[idx].empty()) {
            ret_vid = queues[idx].front();
            queues[idx].pop_front();
            const bool get_success = vfun_set.test_and_get(ret_vid, ret_fun);
            ASSERT_TRUE(get_success);
            locks[idx].unlock();
            return sched_status::NEW_TASK;          
          }
          locks[idx].unlock();
        }
      }
      return sched_status::EMPTY;     
    } // end of get_next_task

    iterminator& terminator() { return term; }

    static void print_options_help(std::ostream& out) { 
      out << "\t mult=3: number of queues per thread." << std::endl;
    }


  }; // end of fifo scheduler 


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

