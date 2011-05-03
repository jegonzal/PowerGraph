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

#ifndef GRAPHLAB_RANDOM_HPP
#define GRAPHLAB_RANDOM_HPP

#include <cstdlib>
#include <stdint.h>


#include <vector>
#include <limits>
#include <algorithm>

#include <boost/random.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {

  /**
   * \ingroup random
   * A collection of thread safe random number routines.  Each thread
   * is assigned its own generator however assigning a seed affects
   * all current and future generators.
   */
  namespace random {        


    ///////////////////////////////////////////////////////////////////////
    //// Underlying generator definition



    namespace distributions {
      /**
       * The uniform distribution struct is used for partial function
       * specialization. Generating uniform random real numbers is
       * accomplished slightly differently than for integers.
       * Therefore the base case is for integers and we then
       * specialize the two real number types (floats and doubles).
       */
      template<typename IntType>
      struct uniform {
        typedef boost::uniform_int<IntType> distribution_type;
        template<typename RealRNG, typename DiscreteRNG>
        static inline IntType sample(RealRNG& real_rng, 
                                     DiscreteRNG& discrete_rng, 
                                     const IntType& min, const IntType& max) {
          return distribution_type(min, max)(discrete_rng);
        }
      };
      template<>
      struct uniform<double> {
        typedef boost::uniform_real<double> distribution_type;
        template<typename RealRNG, typename DiscreteRNG>
        static inline double sample(RealRNG& real_rng, 
                                    DiscreteRNG& discrete_rng, 
                                    const double& min, const double& max) {
          return distribution_type(min, max)(real_rng);
        }
      };
      template<>
      struct uniform<float> {
        typedef boost::uniform_real<float> distribution_type;
        template<typename RealRNG, typename DiscreteRNG>
        static inline float sample(RealRNG& real_rng, 
                                  DiscreteRNG& discrete_rng, 
                                  const float& min, const float& max) {
          return distribution_type(min, max)(real_rng);
        }
      };

    };

    /**
     * The generator class is the base underlying type used to
     * generate random numbers.  User threads should use the functions
     * provided in the random namespace.
     */
    class generator {
    public:
      // base Generator types
      typedef boost::lagged_fibonacci607 real_rng_type;
      typedef boost::mt11213b            discrete_rng_type;  
      typedef boost::rand48              fast_discrete_rng_type;       
    
      //! Seed the generator using the default seed
      inline void seed() {
        mut.lock();
        real_rng.seed();
        discrete_rng.seed();
        fast_discrete_rng.seed();
        mut.unlock();
      }

      //! Seed the generator nondeterministically
      void nondet_seed();


      //! Seed the generator using the current time in microseconds
      inline void time_seed() {
        seed( graphlab::timer::usec_of_day() );
      }

      //! Seed the random number generator based on a number
      void seed(size_t number) {
        mut.lock();
        fast_discrete_rng.seed(number);
        real_rng.seed(fast_discrete_rng);
        discrete_rng.seed(fast_discrete_rng);
        mut.unlock();
      }
      
      //! Seed the generator using another generator
      void seed(generator& other){
        mut.lock();
        real_rng.seed(other.real_rng);
        discrete_rng.seed(other.discrete_rng);
        fast_discrete_rng.seed(other.fast_discrete_rng());
        mut.unlock();
      } 
   
      /**
       * Generate a random number in the uniform real with range [min,
       * max) or [min, max] if the number type is discrete.
       */
      template<typename NumType>
      inline NumType uniform(const NumType min, const NumType max) { 
        mut.lock();
        const NumType result = distributions::uniform<NumType>::
          sample(real_rng, discrete_rng, min, max);
        mut.unlock();
        return result;
      } // end of uniform

      /**
       * Generate a random number in the uniform real with range [min,
       * max) or [min, max] if the number type is discrete.
       */
      template<typename NumType>
      inline NumType fast_uniform(const NumType min, const NumType max) { 
        mut.lock();
        const NumType result = distributions::uniform<NumType>::
          sample(real_rng, fast_discrete_rng, min, max);
        mut.unlock();
        return result;
      } // end of fast_uniform


      /**
       * Generate a random number in the uniform real with range [min,
       * max);
       */
      inline double gamma(const double alpha = double(1)) {
        boost::gamma_distribution<double> gamma_dist(alpha);
        mut.lock();
        const double result = gamma_dist(real_rng);
        mut.unlock();
        return result;
      } // end of gamma

      /**
       * Generate a gaussian random variable with zero mean and unit
       * variance.
       */
      inline double gaussian(const double mean = double(0), 
                             const double var = double(1)) {
        boost::normal_distribution<double> normal_dist(mean,var);
        mut.lock();
        const double result = normal_dist(real_rng);
        mut.unlock();
        return result;
      } // end of gaussian

      inline bool bernoulli(const double p = double(0.5)) {
        boost::bernoulli_distribution<double> dist(p);
        mut.lock();
        const double result(dist(discrete_rng));
        mut.unlock();
        return result;
      } // end of bernoulli

      inline bool fast_bernoulli(const double p = double(0.5)) {
        boost::bernoulli_distribution<double> dist(p);
        mut.lock();
        const double result(dist(fast_discrete_rng));
        mut.unlock();
        return result;
      } // end of bernoulli


      /**
       * Draw a random number from a multinomial
       */
      size_t multinomial(const std::vector<double>& prb) {
        ASSERT_GT(prb.size(),0);
        if (prb.size() == 1) { return 0; }
        double sum(0);
        for(size_t i = 0; i < prb.size(); ++i) {
          ASSERT_GE(prb[i], 0); // Each entry must be P[i] >= 0
          sum += prb[i];
        }
        ASSERT_GT(sum, 0); // Normalizer must be positive
        // actually draw the random number
        const double rnd(uniform<double>(0,1));
        size_t ind = 0;
        for(double cumsum(prb[ind]/sum); 
            rnd > cumsum && (ind+1) < prb.size(); 
            cumsum += (prb[++ind]/sum));
        return ind;
      } // end of multinomial

      
      /** 
       * Shuffle a standard vector
       */ 
      template<typename T>
      void shuffle(std::vector<T>& vec) { shuffle(vec.begin(), vec.end()); }

      /** 
       * Shuffle a range using the begin and end iterators
       */ 
      template<typename Iterator>
      void shuffle(Iterator begin, Iterator end) {
        mut.lock();
        shuffle_functor functor(*this);
        std::random_shuffle(begin, end, functor);
        mut.unlock();
      }

    private:
      //////////////////////////////////////////////////////
      /// Data members
      struct shuffle_functor {
        generator& gen;
        inline shuffle_functor(generator& gen) : gen(gen) { }
        inline std::ptrdiff_t operator()(std::ptrdiff_t end) {
          return distributions::uniform<ptrdiff_t>::
            sample(gen.real_rng, gen.fast_discrete_rng, 0, end-1);
        }
      };

      
      //! The real random number generator
      real_rng_type real_rng;
      //! The discrete random number generator
      discrete_rng_type discrete_rng;
      //! The fast discrete random number generator
      fast_discrete_rng_type fast_discrete_rng;
      //! lock used to access local members
      mutex mut;      
    }; // end of class generator















    /**
     * \ingroup random
     * Seed all generators using the default seed
     */
    void seed();

    /**
     * \ingroup random
     * Seed all generators using an integer
     */
    void seed(size_t seed_value);

    /**
     * \ingroup random
     * Seed all generators using a nondeterministic source
     */
    void nondet_seed();

    /**
     * \ingroup random
     * Seed all generators using the current time in microseconds
     */
    void time_seed();
    

    /**
     * \ingroup random
     * Get the local generator
     */
    generator& get_source();

    /**
     * \ingroup random
     * Generate a random number in the uniform real with range [min,
     * max) or [min, max] if the number type is discrete.
     */
    template<typename NumType>
    inline NumType uniform(const NumType min, const NumType max) { 
      return get_source().uniform<NumType>(min, max);
    } // end of uniform
    
    /**
     * \ingroup random
     * Generate a random number in the uniform real with range [min,
     * max) or [min, max] if the number type is discrete.
     */
    template<typename NumType>
    inline NumType fast_uniform(const NumType min, const NumType max) { 
      return get_source().fast_uniform<NumType>(min, max);
    } // end of fast_uniform
    
    /**
     * \ingroup random
     * Generate a random number between 0 and 1
     */
    inline double rand01() { return uniform<double>(0, 1); }

    /**
     * \ingroup random
     * Simulates the standard rand function as defined in cstdlib
     */
    inline int rand() { return fast_uniform(0, RAND_MAX); }


    /**
     * \ingroup random
     * Generate a random number from a gamma distribution.
     */
    inline double gamma(const double alpha = double(1)) {
      return get_source().gamma(alpha);
    }



    /**
     * \ingroup random
     * Generate a gaussian random variable with zero mean and unit
     * variance.
     */
    inline double gaussian(const double mean = double(0), 
                    const double var = double(1)) {
      return get_source().gaussian(mean, var);
    }

    /**
     * \ingroup random
     * Draw a sample from a bernoulli distribution
     */
    inline bool bernoulli(const double p = double(0.5)) {
      return get_source().bernoulli(p);
    }

    /**
     * \ingroup random
     * Draw a sample form a bernoulli distribution using the faster generator
     */
    inline bool fast_bernoulli(const double p = double(0.5)) {
      return get_source().fast_bernoulli(p);
    }

    /**
     * \ingroup random
     * Generate a draw from a multinomial.  This function
     * automatically normalizes as well.
     */
    inline size_t multinomial(const std::vector<double>& prb) {
      return get_source().multinomial(prb);
    }


    /** 
     * \ingroup random
     * Shuffle a standard vector
     */ 
    template<typename T>
    void shuffle(std::vector<T>& vec) { 
      get_source().shuffle(vec); 
    }

    /** 
     * \ingroup random
     * Shuffle a range using the begin and end iterators
     */ 
    template<typename Iterator>
    void shuffle(Iterator begin, Iterator end) {
      get_source().shuffle(begin, end);
    }

  }; // end of random 
}; // end of graphlab


#endif
