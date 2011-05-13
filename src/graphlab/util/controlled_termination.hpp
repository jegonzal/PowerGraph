/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
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

