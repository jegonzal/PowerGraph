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


#ifndef GRAPHLAB_TASK_COUNT_TERMINATION_HPP
#define GRAPHLAB_TASK_COUNT_TERMINATION_HPP

#include <cassert>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


namespace graphlab {
  /**
   *  \ingroup util_internal
   *  termination checker based on task counting.
   * simple condition variable based shared termination checker.
   * When a processor finds that it is out of work, it should call
   * - begin_critical_section(cpuid),
   * - check the state of the queue.
   * - If the queue has jobs, call cancel_critical_section().
   * - If the queue has no jobs, then call end_critical_section(cpuid)
   * - If (end_critical_section() returns true, the scheduler can terminate.
   * Otherwise it must loop again.
   */
  class task_count_termination {
    atomic<size_t> newtaskcount;
    atomic<size_t> finishedtaskcount;
    bool force_termination; //signal computation is aborted

  public:
    task_count_termination() : newtaskcount(0), 
                               finishedtaskcount(0), 
                               force_termination(false) { }
    
    ~task_count_termination(){ }

    void begin_critical_section(size_t cpuid) { }
    void cancel_critical_section(size_t cpuid)  { }
    
    bool end_critical_section(size_t cpuid) {
      return newtaskcount.value == finishedtaskcount.value 
        || force_termination;
    }
    
    void abort(){
      force_termination = true;
    }

    bool is_aborted(){ return force_termination; }


    void new_job() {
      newtaskcount.inc();
    }
    
    void new_job(size_t cpuhint) {
      new_job();
    }
    
    void completed_job() {
      finishedtaskcount.inc();
      assert(finishedtaskcount.value <= newtaskcount.value);
    }
    
    void print() {
      std::cout << finishedtaskcount.value << " of "
                << newtaskcount.value << std::endl;
    }
  };

}
#endif

