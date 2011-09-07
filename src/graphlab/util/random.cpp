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


#include <pthread.h>

#include <set>
#include <iostream>
#include <fstream>

#include <boost/random.hpp>
#include <boost/integer_traits.hpp>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/timer.hpp>

#include <graphlab/util/random.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {
  namespace random {

    /**
     * A truely nondeterministic generator
     */
    class nondet_generator {
    public:
      typedef size_t result_type;
      BOOST_STATIC_CONSTANT(result_type, min_value = 
                            boost::integer_traits<result_type>::const_min);
      BOOST_STATIC_CONSTANT(result_type, max_value = 
                            boost::integer_traits<result_type>::const_max);
      result_type min BOOST_PREVENT_MACRO_SUBSTITUTION () const { return min_value; }
      result_type max BOOST_PREVENT_MACRO_SUBSTITUTION () const { return max_value; }
      
      nondet_generator() {
        rnd_dev.open("/dev/urandom", std::ios::binary | std::ios::in);
        ASSERT_TRUE(rnd_dev.good());
      }
      // Close the random number generator
      ~nondet_generator() { rnd_dev.close(); }
      // read a size_t from the source
      result_type operator()() {
        // read a machine word into result
        result_type result(0);
        mut.lock();
        ASSERT_TRUE(rnd_dev.good());
        rnd_dev.read(reinterpret_cast<char*>(&result), sizeof(result_type));
        ASSERT_TRUE(rnd_dev.good());
        mut.unlock();
        //        std::cout << result << std::endl;
        return result;
      }      
    private:
      std::ifstream rnd_dev;
      mutex mut;
    };
    nondet_generator global_nondet_rng;






    /**
     * This class represents a master registery of all active random
     * number generators
     */
    struct source_registry {
      std::set<generator*> generators;
      generator master;
      mutex mut;
      
      /**
       * Seed all threads using the default seed
       */
      void seed() {
        mut.lock();
        master.seed();
        foreach(generator* generator, generators) {
          ASSERT_TRUE(generator != NULL);
          generator->seed(master);
        }
        mut.unlock();
      }

      /**
       * Seed all threads using the default seed
       */
      void nondet_seed() {
        mut.lock();
        master.nondet_seed();
        foreach(generator* generator, generators) {
          ASSERT_TRUE(generator != NULL);
          generator->seed(master);
        }
        mut.unlock();
      }


      /**
       * Seed all threads using the default seed
       */
      void time_seed() {
        mut.lock();
        master.time_seed();
        foreach(generator* generator, generators) {
          ASSERT_TRUE(generator != NULL);
          generator->seed(master);
        }
        mut.unlock();
      }

      
      /**
       *  Seed all threads with a fixed number
       */     
      void seed(const size_t number) {
        mut.lock();
        master.seed(number);
        foreach(generator* generator, generators) {
          ASSERT_TRUE(generator != NULL);
          generator->seed(master);
        }
        mut.unlock();
      }
      
      /**
       * Register a source with the registry and seed it based on the
       * master.
       */
      void register_generator(generator* tls_ptr) {
        ASSERT_TRUE(tls_ptr != NULL);
        mut.lock();
        generators.insert(tls_ptr);
        tls_ptr->seed(master);
        mut.unlock();
      }
      
      /**
       * Unregister a source from the registry
       */
      void unregister_source(generator* tls_ptr) {
        mut.lock();
        generators.erase(tls_ptr);
        mut.unlock();
      }
    };
    source_registry registry;








    //////////////////////////////////////////////////////////////
    /// Pthread TLS code

    /**
     * this function is responsible for destroying the random number
     * generators
     */
    void destroy_tls_data(void* ptr) {
      generator* tls_rnd_ptr = 
        reinterpret_cast<generator*>(ptr);
      if(tls_rnd_ptr != NULL) { 
        registry.unregister_source(tls_rnd_ptr);
        delete tls_rnd_ptr; 
      }
    }


    /**
     * Simple struct used to construct the thread local storage at
     * startup.
     */
    struct tls_key_creator {
      pthread_key_t TLS_RANDOM_SOURCE_KEY;
      tls_key_creator() : TLS_RANDOM_SOURCE_KEY(0) {
        pthread_key_create(&TLS_RANDOM_SOURCE_KEY,
                           destroy_tls_data);
      }
    }; 
    /**
     * This static constant instantiates forces the pthread key
     * allocation.
     */
    const tls_key_creator key;     

    /////////////////////////////////////////////////////////////
    //// Implementation of header functions
    
       
 
    generator& get_source() {
      // get the thread local storage
      generator* tls_rnd_ptr = 
        reinterpret_cast<generator*>
        (pthread_getspecific(key.TLS_RANDOM_SOURCE_KEY));
      // Create a tls_random_source if none was provided
      if(tls_rnd_ptr == NULL) {
        tls_rnd_ptr = new generator();      
        assert(tls_rnd_ptr != NULL);
        // This will seed it with the master rng
        registry.register_generator(tls_rnd_ptr);
        pthread_setspecific(key.TLS_RANDOM_SOURCE_KEY, 
                            tls_rnd_ptr);      
      }
      // assert(tls_rnd_ptr != NULL);
      return *tls_rnd_ptr;
    } // end of get local random source



    void seed() { registry.seed();  } 

    void nondet_seed() { registry.nondet_seed(); } 

    void time_seed() { registry.time_seed(); } 

    void seed(const size_t seed_value) { registry.seed(seed_value);  } 


    void generator::nondet_seed() {
      mut.lock();
      // std::cout << "initializing real rng" << std::endl;
      real_rng.seed(global_nondet_rng());
      // std::cout << "initializing discrete rng" << std::endl;
      discrete_rng.seed(global_nondet_rng());
      // std::cout << "initializing fast discrete rng" << std::endl;
      fast_discrete_rng.seed(global_nondet_rng());
      mut.unlock();
    }



  
  }; // end of namespace random

};// end of namespace graphlab

