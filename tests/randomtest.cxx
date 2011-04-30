

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <graphlab.hpp>


using namespace graphlab;

typedef double vertex_data_type;
typedef double edge_data_type;
typedef graphlab::core<vertex_data_type, edge_data_type> core_type;

template<typename NumType>
void uniform_speed(const size_t max_iter) {
  NumType sum(0);
  timer ti;
  ti.start();
  for(size_t i = 0; i < max_iter; ++i) {
    sum += graphlab::random::uniform<NumType>(0, 10);
  }
  double slow_time = ti.current_time();
  ti.start();
  for(size_t i = 0; i < max_iter; ++i) {
    sum += graphlab::random::fast_uniform<NumType>(0, 10);
  }
  double fast_time = ti.current_time();
  std::cout << slow_time << ", " << fast_time << std::endl; 
}


static void update_function(core_type::types::iscope& scope,
                            core_type::types::icallback& callback) {
  namespace random = graphlab::random;

  scope.vertex_data() += 
    random::uniform<int>(0,9) +
    random::fast_uniform<int>(0,9) +
    random::uniform<char>(0, 2) +
    random::fast_uniform<char>(0, 2) +
    random::uniform<uint32_t>(0, 5) +
    random::fast_uniform<uint32_t>(0, 5) +
    random::rand() % 10 +
    random::uniform<size_t>(0,9) +
    random::fast_uniform<size_t>(0, 10) +
    random::uniform<double>(0,1) +
    random::fast_uniform<double>(0,1) +
    random::uniform<float>(0,1) +
    random::fast_uniform<float>(0,1) +
    random::gamma() +
    random::gaussian();
  std::vector<double> weights(5);
  for(size_t i = 0; i < weights.size(); ++i) {
    weights[i] = random::uniform<double>(0,1);
  }
  scope.vertex_data() += random::multinomial(weights);                 

}




class thread_worker {
public:
  std::vector<int> values;
  void run() {
    namespace random = graphlab::random;
    for(size_t i = 0; i < values.size(); ++i) {
      values[i] = random::uniform<int>(0,3);
    }
  }
};

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& values) {
  out << "{";
  for(size_t i = 0; i < values.size(); ++i) {
    out << values[i];
    if(i + 1 < values.size()) out << ", ";
  }
  return out << "}";
}


std::vector<int> operator+(const std::vector<int>& v1, 
                           const std::vector<int>& v2) {
  assert(v1.size() == v2.size());
  std::vector<int> result(v1.size());
  for(size_t i = 0; i < result.size(); ++i) {
    result[i] = v1[i] + v2[i];
  }
  return result;
}




class RandomTestSuite: public CxxTest::TestSuite {
  size_t iterations;
 
  
  public:

  RandomTestSuite() : iterations(1E8) { }
  
 

  void test_random_number_generators() {
    std::cout << std::endl;
    std::cout << "beginning seed" << std::endl;
    namespace random = graphlab::random;
    graphlab::random::seed();
    graphlab::random::time_seed();
    graphlab::random::nondet_seed();
    graphlab::random::seed(12345);
    std::cout << "finished" << std::endl;

    const size_t num_iterations(20);
    std::vector<thread_worker> workers(10);
    for(size_t i = 0; i < workers.size(); ++i) 
      workers[i].values.resize(num_iterations);
    thread_group threads;
    for(size_t i = 0; i < workers.size(); ++i) {
      threads.launch(boost::bind(&thread_worker::run, &(workers[i])));
    }
    threads.join();
    for(size_t i = 0; i < workers.size(); ++i) {
      std::cout << workers[i].values << std::endl;
    }
    std::vector<int> sum(workers[0].values.size());
    for(size_t i = 0; i < workers.size(); ++i) {
      sum = sum + workers[i].values;
    }
    std::cout << "Result: " << sum << std::endl;
  }




  void test_randomness_in_engine() {
    core_type core;
    for(size_t i = 0; i < 32; ++i) 
      core.graph().add_vertex(vertex_data_type(0));
    for(size_t i = 0; i+1 < core.graph().num_vertices(); ++i) {
      core.graph().add_edge(i, i+1, edge_data_type(1));
      core.graph().add_edge(i+1, i, edge_data_type(2));
    }
    core.graph().compute_coloring();
    core.set_scheduler_type("chromatic");
    core.set_scope_type("none");
    core.set_ncpus(4);
    core.sched_options().add_option("update_function", update_function);
    for(size_t i = 0; i < 3; ++i) {
      core.sched_options().add_option("max_iterations", 2);
      std::cout << "--------------------------------" << std::endl;
      core.start();
      for(size_t i = 0; i < core.graph().num_vertices(); ++i) {
        std::cout << core.graph().vertex_data(i) << "\t";
      }
      std::cout << std::endl;
    }
  }

  void test_shuffle() {
    namespace random = graphlab::random;
    random::nondet_seed();
    std::vector<int> numbers(6);
    for(size_t i = 0; i < numbers.size(); ++i) numbers[i] = i + 1;
    for(size_t j = 0; j < 10; ++j) {
      // shuffle the numbers
      random::shuffle(numbers);
      std::cout << numbers << std::endl;
    }

  }



  void test_speed() {
    namespace random = graphlab::random;
    std::cout << "speed test run: " << std::endl;
    const size_t MAX_ITER(10000000);
    std::cout << "size_t:   "; 
    uniform_speed<size_t>(MAX_ITER);
    std::cout << "int:      "; 
    uniform_speed<int>(MAX_ITER);
    std::cout << "uint32_t: "; 
    uniform_speed<uint32_t>(MAX_ITER);
    std::cout << "uint16_t: "; 
    uniform_speed<uint16_t>(MAX_ITER);
    std::cout << "char:     "; 
    uniform_speed<char>(MAX_ITER);
    std::cout << "float:    "; 
    uniform_speed<float>(MAX_ITER);
    std::cout << "double:   "; 
    uniform_speed<double>(MAX_ITER);
    
    std::cout << "gaussian: ";
    double sum = 0;
    timer time;
    time.start();
    for(size_t i = 0; i < MAX_ITER; ++i) 
      sum += random::gaussian();
    std::cout << time.current_time() << std::endl;
    
    std::cout << "shuffle:  "; 
    std::vector<int> numbers(6);
    for(size_t i = 0; i < numbers.size(); ++i) numbers[i] = i + 1;
    time.start();
    for(size_t j = 0; j < MAX_ITER/numbers.size(); ++j) {
      // shuffle the numbers
      random::shuffle(numbers);
    }
    std::cout << time.current_time() << ", ";
    time.start();
    for(size_t j = 0; j < MAX_ITER/numbers.size(); ++j) {
      // shuffle the numbers
      std::random_shuffle(numbers.begin(), numbers.end());
    }
    std::cout << time.current_time() << std::endl;
    
    
    
    
  }


  
  
};
