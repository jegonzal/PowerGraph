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

#include <cstring>
#include <queue>
#include <cmath>
#include <cassert>

#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/random.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
#include <graphlab/schedulers/support/binary_vertex_task_set.hpp>
//#include <graphlab/util/shared_termination.hpp>
#include <graphlab/util/task_count_termination.hpp>
#include <graphlab/metrics/metrics.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \ingroup group_schedulers
     This class defines a simple First-In-First_Out scheduler
   */
  template<typename Engine>
  class multiqueue_fifo_scheduler : public ischeduler<Engine> {
  
  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef ischeduler<Engine> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::engine_type engine_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;


    typedef std::queue<vertex_id_type> taskqueue_t;
    typedef task_count_termination terminator_type;

    

  private:

    vertex_functor_set<engine_type> vfun_set;

    size_t num_queues;
    size_t queues_per_cpu;
  
    std::vector<taskqueue_t> task_queues; /// The actual task queue
    std::vector<spinlock> queue_locks;
    std::vector<size_t> lastqueue;


    // Terminator
    task_count_termination terminator;
 


  public:

    multiqueue_fifo_scheduler(const graph_type& graph, 
                              size_t ncpus,
                              const options_map& opts) :
      vfun_set(graph.num_vertices()) { 
      queues_per_cpu = 2;
      num_queues = queues_per_cpu * ncpus;
       
      /* Each cpu keeps record of the queue it last 
         used to keep balance */
      lastqueue.resize(ncpus, 0);
       
      // Do this in the preconstructor
      task_queues.resize(num_queues);
      queue_locks.resize(num_queues);

    }

    void start() { term.reset(); }
   

    void schedule(vertex_id_type vid, 
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

 
    /** Get the next element in the queue */
    sched_status::status_enum get_next_task(size_t cpuid,
                                            update_task_type &ret_task) {
      bool found = false;
      /* First check my own queues. Keep track which own queue was checked
         so next time I can check next of my own queues to keep balance. */
      size_t firstown = cpuid * queues_per_cpu;
      for(size_t ownq_i = 0; ownq_i < queues_per_cpu; ++ownq_i) {
        size_t queueidx = 
          firstown + ((ownq_i + lastqueue[cpuid] + 1) % queues_per_cpu);
        taskqueue_t& queue = task_queues[queueidx];
        queue_locks[queueidx].lock();
        if (!queue.empty()) {
          ret_task = queue.front();
          queue.pop();
          found = true;
          lastqueue[cpuid] = ownq_i;
        }
        queue_locks[queueidx].unlock();
        if (found) break;
      }
  
      /* Ok, my queues were empty - now check every other queue */
      if (!found) {
        /* First check own queue - if it is empty, check others */
        for(size_t roundrobin = 0; roundrobin < num_queues; ++roundrobin) {
          size_t queueidx = 
            (firstown + queues_per_cpu + roundrobin) % num_queues;
          taskqueue_t& queue = task_queues[queueidx];
          queue_locks[queueidx].lock();
          if (!queue.empty()) {
            ret_task = queue.front();
            queue.pop();
            found = true;
          }
          queue_locks[queueidx].unlock();
          if (found)  break;
        }
      }
 
      if(!found) {
        return sched_status::EMPTY;
      }
      
      binary_vertex_tasks.remove(ret_task);
      
      if (monitor != NULL) 
        monitor->scheduler_task_scheduled(ret_task, 0.0);
      return sched_status::NEWTASK;
    } // end of get_next_task


    void add_task(update_task_type task, double priority) {
      if (binary_vertex_tasks.add(task)) {
        terminator.new_job();
        // Check if task should be pruned
        /* "Randomize" the task queue task is put in. Note that we do
           not care if this counter is corrupted in race conditions */
    
   
        /* Find first queue that is not locked and put task there (or
           after iteration limit)*/
         
        /* Choose two random queues and use the one which has smaller
           size */
        // M.D. Mitzenmacher The Power of Two Choices in Randomized
        // Load Balancing (1991)
        // http://www.eecs.harvard.edu/~michaelm/postscripts/mythesis.pdf

//         size_t r1 = random::rand_int(num_queues - 1);
//         size_t r2 = random::rand_int(num_queues - 1);
        // TODO: this can easily be done with some bit operations on the 
        const size_t prod = 
          random::fast_uniform(size_t(0), num_queues * num_queues - 1);
        const size_t r1 = prod / num_queues;
        const size_t r2 = prod % num_queues;

        size_t qidx = 
          (task_queues[r1].size() < task_queues[r2].size()) ? r1 : r2;
        
        queue_locks[qidx].lock();
        task_queues[qidx].push(task);
        queue_locks[qidx].unlock();
    
        if (monitor != NULL) 
          monitor->scheduler_task_added(task, priority);
      } else {
        prunecounter[thread::thread_id()]++;
        if (monitor != NULL) 
          monitor->scheduler_task_pruned(task);
      }
  
    }

    void add_tasks(const std::vector<vertex_id_type> &vertices,
                   update_function_type func,
                   double priority) {
      foreach(vertex_id_type vertex, vertices) {
        add_task(update_task_type(vertex, func), priority);
      }
    }


    void add_task_to_all(update_function_type func, double priority)  {
      for (vertex_id_type vertex = 0; vertex < numvertices; ++vertex){
        add_task(update_task_type(vertex, func), priority);
      }
    }



    void completed_task(size_t cpuid, const update_task_type &task) {
      terminator.completed_job();
    }


    bool is_task_scheduled(update_task_type task)  {
      return binary_vertex_tasks.get(task);
    }


    void print() {
      std::cout << "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS" << std::endl;
      std::cout << "Printing task queue sizes: " << std::endl;
      for(size_t i = 0; i < task_queues.size(); ++i) {
        std::cout << task_queues[i].size() << std::endl;
      }
      std::cout << "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS" << std::endl;
        
    }

    terminator_type& get_terminator() {
      return terminator;
    };

    void set_options(const scheduler_options &opts) { }

    metrics get_metrics() {
      for(unsigned int i=0; i<prunecounter.size(); i++) sched_metrics.add("pruned", (double)prunecounter[i], INTEGER); 
      return sched_metrics;
    }

    void reset_metrics() {
      for(unsigned int i=0; i<prunecounter.size(); i++) prunecounter[i] = 0;
      sched_metrics.clear();
    }


  }; 


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

