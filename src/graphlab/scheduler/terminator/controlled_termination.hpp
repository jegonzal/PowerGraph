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


#ifndef GRAPHLAB_CONTROLLED_TERMINATION_HPP
#define GRAPHLAB_CONTROLLED_TERMINATION_HPP

#include <cassert>

#include <graphlab/scheduler/terminator/iterminator.hpp>

namespace graphlab {
  /**
   *  \ingroup util_internal
   * The simplest possible terminator.
   * Always fails until a flag is set
   */
  class controlled_termination : public iterminator {
    bool natural_quit;
    bool forced_abort;
  public:
    controlled_termination() : 
      natural_quit(false), forced_abort(false) { }


    void begin_critical_section(size_t cpuid) { }
    void cancel_critical_section(size_t cpuid)  { }
    bool end_critical_section(size_t cpuid) { 
      return natural_quit || forced_abort;  }
    void new_job(size_t cpuhint) { }
    void completed_job() { }
    bool is_aborted() { return forced_abort; }
    void abort() { forced_abort = true; }
    void reset() { 
      natural_quit = false; 
      forced_abort = false;
    }
    void complete() { natural_quit = true; }

  };

}
#endif

