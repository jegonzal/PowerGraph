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


#ifndef CONTROLLED_TERMINATION_HPP
#define CONTROLLED_TERMINATION_HPP

#include <cassert>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


namespace graphlab {
  /**
   *  \ingroup util_internal
   * The simplest possible terminator.
   * Always fails until a flag is set
   */
  class controlled_termination {
    bool quit;
  public:
    controlled_termination():quit(false) { }

    ~controlled_termination(){ }

    void begin_critical_section(size_t cpuid) { }
    void cancel_critical_section(size_t cpuid)  { }

    bool end_critical_section(size_t cpuid) { return quit;  }

    void new_job() { }

    void new_job(size_t cpuhint) { }

    void completed_job() { }

    void complete() { quit = true; }

    void reset() { quit = false; }
  };

}
#endif

