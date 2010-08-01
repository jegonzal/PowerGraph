

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <boost/random.hpp>
#include <util/timer.hpp>
#include <parallel/pthread_tools.hpp>

using namespace graphlab;

class RandomTestSuite: public CxxTest::TestSuite {
  size_t iterations;
 
  
  public:

  RandomTestSuite() : iterations(1E8) { }
  
  void test_laggedfib(void) {
    typedef boost::lagged_fibonacci607 generator_type;
    boost::lagged_fibonacci607 rng;
    boost::uniform_01<boost::lagged_fibonacci607>  unifrand(rng);
    
    double sum = 0.0;
    timer tmr;
    tmr.start();
    for(size_t i = 0; i < iterations; ++i) {
      sum += unifrand();
    }
    double elapsed_time = tmr.current_time();
    std::cout << std::endl;
    std::cout << "-----------------------------------------------"
              << std::endl;
    std::cout << "Fib Runtime: " << elapsed_time <<  std::endl;
    std::cout << "Total Value: " << sum <<  std::endl;    
  }

  void test_rand(void) {
    double sum = 0.0;
    timer tmr;
    tmr.start();
    for(size_t i = 0; i < iterations; ++i) {
      sum += double(rand()) / RAND_MAX;
    }
    double elapsed_time = tmr.current_time();
    std::cout << std::endl;
    std::cout << "-----------------------------------------------"
              << std::endl;
    std::cout << "Rand Runtime: " << elapsed_time << std::endl;
    std::cout << "Total Value: " << sum <<  std::endl;    
  }

  // Threaded random tests
  struct worker_laggedfib : public runnable {
    double sum;
    size_t iterations;
    int seed;
    void run() {
      boost::lagged_fibonacci607 fib(seed);
      boost::uniform_01<boost::lagged_fibonacci607>  unifrand(fib);
      for(size_t i = 0; i < iterations; ++i) {
        sum += unifrand();
      }
    }
  };

    // Threaded random tests
  struct worker_rand : public runnable {
    double sum;
    size_t iterations;    
    void run() {
      for(size_t i = 0; i < iterations; ++i) {
        sum += rand();
      }
    }
  };

  void test_threaded_randomness_lagged(void) {
    timer tmr;
    tmr.start();
    std::vector<worker_laggedfib> workers(4);
    thread_group group;
    for(size_t i = 0; i < workers.size(); ++i) {
      workers[i].iterations = iterations;
      workers[i].seed = i + 33;
      group.launch(&(workers[i]));      
    }
    group.join();
    double sum = 0.0;
    for(size_t i = 0; i < workers.size(); ++i) {
      sum += workers[i].sum;
    }
    double elapsed_time = tmr.current_time();
    std::cout << std::endl;
    std::cout << "------------ threaded -------------------------"
              << std::endl;
    std::cout << "Fib Runtime: " << elapsed_time << std::endl;
    std::cout << "Total Value: " << sum <<  std::endl;    
  }

  void test_threaded_randomness_rand(void) {
    timer tmr;
    tmr.start();
    std::vector<worker_rand> workers(4);
    thread_group group;
    for(size_t i = 0; i < workers.size(); ++i) {
      workers[i].iterations = iterations;
      group.launch(&(workers[i]));      
    }
    group.join();
    double sum = 0.0;
    for(size_t i = 0; i < workers.size(); ++i) {
      sum += workers[i].sum;
    }
    double elapsed_time = tmr.current_time();
    std::cout << std::endl;
    std::cout << "------------- threaded ------------------------"
              << std::endl;
    std::cout << "Rand Runtime: " << elapsed_time << std::endl;
    std::cout << "Total Value:  " << sum <<  std::endl;    
  }


  
  
};
