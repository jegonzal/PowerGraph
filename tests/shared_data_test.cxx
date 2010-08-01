#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include <cxxtest/TestSuite.h>

#include <graphlab.hpp>

using namespace graphlab;

typedef graph<double, double> graph_type;
typedef thread_shared_data<graph_type> thread_shared_data_type;
typedef ishared_data<graph_type> ishared_data_type;
typedef ishared_data_type::iscope_type iscope_type;
typedef general_scope_factory<graph_type> scope_factory_type;

struct accumulator_type : public unsupported_serialize  {
  double sum;
  size_t count;
  accumulator_type() : sum(0), count(0) { }
};
  
void sync_sum_fun(size_t index,
                  const ishared_data_type& sdm,
                  iscope_type& iscope,
                  any& acc) {
  acc.as<accumulator_type>().sum += iscope.vertex_data();
  acc.as<accumulator_type>().count++;
}

void apply_fun(size_t index,
               const ishared_data_type& sdm,
               any& current_data,
               const any& acc) {
  current_data.as<double>() =
    acc.as<accumulator_type>().sum /
    acc.as<accumulator_type>().count;
}


class SharedDataTest : public CxxTest::TestSuite {
public:
  
  
  void test_thread_shared_data(void) {
    std::cout << "Testing thread sdm" << std::endl;
    thread_shared_data_type tsdm;
    ishared_data_type& sdm(tsdm); 
    tsdm.create_atomic(3, any(std::string("hello world")));
    std::cout << sdm.get(3).as<std::string>() << std::endl;

    graph_type graph(6);
    std::cout << "Graph Before: " << std::endl;
    std::cout << graph << std::endl;
    for(size_t i = 0; i < 6; ++i) {
      graph.vertex_data(i) = i;
      if(i > 0) graph.add_edge(i-1, i, i - 0.5);
    }
    std::cout << "Graph After: " << std::endl;
    std::cout << graph << std::endl;


    


    
    accumulator_type zero;    
    tsdm.set_sync(1,
                  sync_sum_fun,
                  apply_fun,
                  any(zero),
                  1000);
    tsdm.atomic_set(1, any(double(0)));
    std::cout << "Before: " << sdm.get(1).as<double>()
              << std::endl;
    tsdm.sync(graph, 1);
    std::cout << "After: " << sdm.get(1).as<double>()
              << std::endl;


    

    scope_factory_type scope_factory(graph, 1, scope_range::FULL_CONSISTENCY);
    tsdm.set_scope_factory(&scope_factory);

    double current_value = sdm.get(1).as<double>();

    graph.vertex_data(2) = 7.3;
    std::cout << "Waiting on signal" << std::endl;
    while(sdm.get(1).as<double>() == current_value) {
      tsdm.signal(1);
      usleep(1000);
    }
    std::cout << "After: " << sdm.get(1).as<double>()
              << std::endl;    
    
    
  }


  
};

