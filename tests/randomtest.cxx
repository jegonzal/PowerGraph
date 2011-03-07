

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <boost/random.hpp>

#include <graphlab.hpp>


using namespace graphlab;

typedef size_t vertex_data_type;
typedef size_t edge_data_type;
typedef graphlab::core<vertex_data_type, edge_data_type> core_type;
static void update_function(core_type::types::iscope& scope,
                            core_type::types::icallback& callback,
                            core_type::types::ishared_data* shared_data) {
  
  std::cout << scope.vertex() << ": " << random::rand01() << std::endl;
}


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


 
  
  
  void test_engine() {
    core_type core;
    core.graph().add_vertex(vertex_data_type(1));
    core.graph().add_vertex(vertex_data_type(2));
    core.graph().add_edge(0,1, edge_data_type(1));
    core.graph().add_edge(1,0, edge_data_type(2));
    core.graph().compute_coloring();
    core.set_scheduler_type("chromatic");
    core.set_scope_type("none");
    core.set_ncpus(4);
    core.sched_options().add_option("update_function", update_function);
    for(size_t i = 0; i < 3; ++i) {
      core.sched_options().add_option("max_iterations", 2*(i+1));
      std::cout << "--------------------------------" << std::endl;
      core.start();
    }

  }

  
  
};
