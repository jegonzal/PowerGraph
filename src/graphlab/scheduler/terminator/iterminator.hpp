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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#ifndef GRAPHLAB_ITERMINATOR_HPP
#define GRAPHLAB_ITERMINATOR_HPP

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
  class iterminator {
  public:
    virtual ~iterminator() { }
    virtual void begin_critical_section(size_t cpuid) = 0;
    virtual void cancel_critical_section(size_t cpuid) = 0;
    virtual bool end_critical_section(size_t cpuid) = 0;
    virtual void new_job(size_t cpuhint = size_t(-1)) = 0;
    virtual void completed_job() = 0;
    virtual bool is_aborted() = 0;
    virtual void abort() = 0;
    virtual void reset() = 0;


  };

}
#endif

