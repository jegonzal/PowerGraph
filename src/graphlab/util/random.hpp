#ifndef GRAPHLAB_RANDOM_HPP
#define GRAPHLAB_RANDOM_HPP

#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>

#include <graphlab/macros_def.hpp>
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
     * Generate a draw from a multinomial.  This function
     * automatically normalizes as well.
     */
    static size_t rand_multi(const std::vector<double>& prb) {
      ASSERT_GT(prb.size(),0);
      if (prb.size() == 1) {
	return 0;
      }
      double sum(0);
      for(size_t i = 0; i < prb.size(); ++i) {
        ASSERT_GE(prb[i], 0);
        sum += prb[i];
      }
      ASSERT_GT(sum, 0);
      const double rnd(rand01());
      size_t ind = 0;
      for(double cumsum(prb[ind]/sum); 
          rnd > cumsum && (ind+1) < prb.size(); 
          cumsum += (prb[++ind]/sum));
      return ind;
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
#include <graphlab/macros_undef.hpp>

#endif
