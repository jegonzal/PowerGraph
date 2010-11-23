#ifndef CONTROLLED_TERMINATION_HPP
#define CONTROLLED_TERMINATION_HPP

#include <cassert>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


namespace graphlab {
  /**
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

    bool end_critical_section(size_t cpuid) {
      return quit;
    }

    void new_job() { }

    void new_job(size_t cpuhint) { }

    void completed_job() { }

    void complete() {
      quit = true;
    }
  };

}
#endif
