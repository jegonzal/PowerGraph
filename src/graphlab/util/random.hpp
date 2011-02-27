#ifndef GRAPHLAB_RANDOM_HPP
#define GRAPHLAB_RANDOM_HPP

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>

namespace graphlab {

  /**
   * A collection of thread safe random number routines.  Each thread
   * is assigned its own generator.
   */
  struct random {
    
    /**
     * Generate a thread specific random number between 0 and 1
     */
    static double rand01() {
      return thread::get_thread_specific_data().rand01();
    }

    /**
     * Generate a gamma distribution random variable 
     */
    static double rand_gamma(double alpha = 1) {
      return thread::get_thread_specific_data().rand_gamma(alpha);
    }


    /**
     * Generate a gaussian random variable with zero mean and unit
     * variance.
     */
    static double rand_gaussian() {
      return thread::get_thread_specific_data().rand_gaussian();
    }

    /**
     * This function has been deprecated in favor of the rand_gaussian
     * name
     */
     __attribute__((__deprecated__)) static double gaussian_rand() {
      return thread::get_thread_specific_data().rand_gaussian();
    }

    /**
     * Generate a thread specific random integer between 0 and max inclusive
     */
    static size_t rand_int(size_t max) {
      return thread::get_thread_specific_data().rand_int(max);
    }
    
    /**
     * Seed the random number generator.  Note that the engine
     * automatically seeds the number generators.
     */
    static void seed(size_t value) {
      thread::get_thread_specific_data().seed(value);
    }

    /** Seed the random number generator with the default seed */
    static void seed() {
      thread::get_thread_specific_data().seed(timer::usec_of_day());
    }
  }; // end of random 
} // end of graphlab


#endif
