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


#include <iostream>
#include <string>
#include <map>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/fiber_async_consensus.hpp>
#include <graphlab/util/blocking_queue.hpp>
#include <graphlab/parallel/fiber_group.hpp>
using namespace graphlab;

#define NTHREADS 1000


class simple_engine_test {
 public:
  dc_dist_object<simple_engine_test> rmi;
  blocking_queue<size_t> queue;
  fiber_async_consensus cons;
  atomic<size_t> numactive;;

  simple_engine_test(distributed_control &dc):rmi(dc, this), cons(dc, NTHREADS) {
    numactive.value = NTHREADS; 
    dc.barrier();
  }

  void add_task_local(size_t i) {
    queue.enqueue(i);
    if (numactive.value < NTHREADS) cons.cancel();
  }  
  
  void task(size_t i) {
    if (i < 5) std::cout << "Task " << i << std::endl;
    if (i > 0) {
      if (rmi.numprocs() == 1) {
        add_task_local(i - 1);
      }
      else {
        rmi.remote_call((procid_t)((rmi.procid() + 1) % rmi.numprocs()),
                    &simple_engine_test::add_task_local,
                    i - 1);
      }
    }
  }
  
  bool try_terminate(size_t cpuid, std::pair<size_t, bool> &job) {
    job.second = false;
    
    numactive.dec();
    cons.begin_done_critical_section(cpuid);
    job = queue.try_dequeue();
    if (job.second == false) {
      bool ret = cons.end_done_critical_section(cpuid);
      numactive.inc();
      return ret;
    }
    else {
      cons.cancel_critical_section(cpuid);
      numactive.inc();
      return false;
    }
  }
  
  void thread(size_t cpuid) {
    while(1) {
       std::pair<size_t, bool> job = queue.try_dequeue();
       if (job.second == false) {
          bool ret = try_terminate(cpuid, job);
          if (ret == true) break;
          if (ret == false && job.second == false) continue;
       }
       task(job.first);
    }
  }
  
  void start_thread() {
    fiber_group thrgrp; 
    for (size_t i = 0;i < NTHREADS; ++i) {
      thrgrp.launch(boost::bind(
                            &simple_engine_test::thread,
                            this, i));
    }
    
    thrgrp.join();
    ASSERT_EQ(queue.size(), 0);
  }
};


int main(int argc, char ** argv) {
  /** Initialization */
  mpi_tools::init(argc, argv);
  global_logger().set_log_level(LOG_DEBUG);

  dc_init_param param;
  if (init_param_from_mpi(param) == false) {
    return 0;
  }
  distributed_control dc(param);
  simple_engine_test test(dc);
  test.add_task_local(300);
  test.start_thread();
  dc.barrier();
  mpi_tools::finalize();
}
