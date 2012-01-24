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


#ifndef GRAPHLAB_CRITICAL_TERMINATION_HPP
#define GRAPHLAB_CRITICAL_TERMINATION_HPP

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>


namespace graphlab {
  /**
   \ingroup util_internal
     simple condition variable based shared termination checker.
     When a processor finds that it is out of work, it should call
     - begin_critical_section(cpuid),
     - check the state of the queue.
     - If the queue has jobs, call cancel_critical_section().
     - If the queue has no jobs, then call end_critical_section(cpuid)
     - If (end_critical_section() returns true, the scheduler can terminate.
     Otherwise it must loop again.
  */
  class critical_termination : public iterminator {
  public:
    critical_termination(size_t ncpus) :
      cond(ncpus),
      numactive(ncpus),
      ncpus(ncpus),
      done(false),
      trying_to_sleep(0),
      critical(ncpus, 0),
      sleeping(ncpus, 0) { }
  
  
    void begin_critical_section(size_t cpuid) {
      trying_to_sleep.inc();
      critical[cpuid] = true;
      m.lock();
    }

    void cancel_critical_section(size_t cpuid) {
      m.unlock();
      critical[cpuid] = false;
      trying_to_sleep.dec();
    }

    bool end_critical_section(size_t cpuid) {
      // if done flag is set, quit immediately
      if (done) {
        m.unlock();
        critical[cpuid] = false;
        trying_to_sleep.dec();
        return true;
      }
      /*
        Assertion: Since numactive is decremented only within 
        a critical section, and is incremented only within the same critical 
        section. Therefore numactive is a valid counter of the number of threads 
        outside of this critical section. 
      */
      --numactive;
    
      /*
        Assertion: If numactive is ever 0 at this point, the algorithm is done.
        WLOG, let the current thread which just decremented numactive be thread 0
      
        Since there is only 1 active thread (0), there must be no threads 
        performing insertions, and are no othe threads which are waking up.
        All threads must therefore be sleeping in cond.wait().
      */
      if (numactive == 0) {
        done = true;
        for (size_t i = 0;i < cond.size(); ++i) cond[i].signal();
      } else {
        sleeping[cpuid] = true;
        while(1) {
          cond[cpuid].wait(m);
          // here we are protected by the mutex again.
          
          // woken up by someone else. leave the 
          // terminator
          if (sleeping[cpuid] == false || done) {
            break;
          }
        }
      }
      m.unlock();
      critical[cpuid] = false;
      trying_to_sleep.dec();
      return done;
    }
  
    /**
      called if a new task is available, and all sleeping CPUs should
      be woken up to handle the job. This is used in the case where
      it is not known which processor is responsible for the new job, 
      or where any processor can handle the job.
    */
    void new_job() {
      /*
        Assertion: numactive > 0 if there is work to do.
        If there are threads trying to sleep, lets wake them up
      */
      if (trying_to_sleep > 0 || numactive < ncpus) {
        m.lock();
        // once I acquire this lock, all threads must be
        // in the following states
        // 1: still running and has not reached begin_critical_section()
        // 2: is sleeping in cond.wait()
        // 3: has called begin_critical_section() but has not acquired
        //    the mutex
        // In the case of 1,3: These threads will perform one more sweep
        // of their task queues. Therefore they will see any new job if available
        // in the case of 2: numactive must be < ncpus since numactive
        // is mutex protected. Then I can wake them up by
        // clearing their sleeping flags and broadcasting.
        if (numactive < ncpus) {
          // this is safe. Note that it is done from within 
          // the critical section.
          for (size_t i = 0;i < ncpus; ++i) {
            numactive += sleeping[i];
            if (sleeping[i]) cond[i].signal();
            sleeping[i] = 0;
          }
        }
        m.unlock();
      }
    }

    /**
      called if a new task is available, and the CPU meant to process
      the job is known. Only the processor [cpuhint] will be woken 
      up to process the job. Note that this is not particular efficient
      since a broadcast is used, all CPUs will actually wake up
      briefly.
    */
    void new_job(size_t cpuhint) {
      if (critical[cpuhint]) {
        m.lock();
        // see new_job() for detailed comments
        if (sleeping[cpuhint]) {
          numactive += sleeping[cpuhint];
          cond[cpuhint].signal();
        }
        m.unlock();
      }
    }

    void completed_job() { }
    
    size_t num_active() {
      return numactive;
    }

    bool is_aborted() {
      return forced_abort;
    }

    void abort() { 
      forced_abort = true;
      done = true;
    }

    void reset() {
      numactive = ncpus;
      done = false;
      forced_abort = false;
      trying_to_sleep.value = 0;
      for (size_t i = 0; i < ncpus; ++i) critical[i] = 0;
    }
  private:
    std::vector<conditional> cond;
    mutex m;
    
    /// counts the number of threads which are not sleeping
    /// protected by the mutex
    size_t numactive; 
    
    /// Total number of CPUs
    size_t ncpus;
    
    /// once flag is set, the terminator is invalid, and all threads
    /// should leave
    bool done;
    
    /// set if abort() is called
    bool forced_abort;    
    
    /// Number of threads which have called
    /// begin_critical_section(), and have not left end_critical_section()
    /// This is an atomic counter and is not protected.
    atomic<size_t> trying_to_sleep;
    
    /// critical[i] is set if thread i has called 
    /// begin_critical_section(), but has not left end_critical_section()
    /// sum of critical should be the same as trying_to_sleep
    std::vector<char> critical;
    
    /// sleeping[i] is set if threads[i] is in cond.wait()
    std::vector<char> sleeping;    
  };

}
#endif

