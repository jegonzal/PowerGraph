#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include <cxxtest/TestSuite.h>

#include <graphlab.hpp>

using namespace graphlab;

typedef graph<double, double> graph_type;
typedef types<graph_type> gl;


glshared<std::string> helloworld;
glshared<double> graphaverage;

struct accumulator_type : public unsupported_serialize  {
  double sum;
  size_t count;
  accumulator_type() : sum(0), count(0) { }
};



void sync_sum_fun(gl::iscope& iscope,
                  any& acc) {
  acc.as<accumulator_type>().sum += iscope.vertex_data();
  acc.as<accumulator_type>().count++;
}

void apply_fun(any& current_data,
               const any& acc) {
  current_data.as<double>() =
    acc.as<accumulator_type>().sum /
    acc.as<accumulator_type>().count;
}


class SharedVariableTest : public CxxTest::TestSuite {
public:
  
  
  void test_shared_variable(void) {
    gl::core core;
    // set initial values
    helloworld.set("hello world");
    ASSERT_TRUE("hello world" == helloworld.get_val());
    graphaverage.set(0);
    ASSERT_EQ(graphaverage.get_val(), 0);
    

    for(size_t i = 0; i < 6; ++i) {
      core.graph().add_vertex(double(i));
    }
    for(size_t i = 1; i < 6; ++i) {
      core.graph().add_edge(i-1, i, i - 0.5);
    }
    
    std::cout << "Graph: " << std::endl;
    std::cout << core.graph() << std::endl;

    accumulator_type zero;    
    core.engine().set_sync(graphaverage,
                           sync_sum_fun,
                           apply_fun,
                           any(zero),
                           1);

    core.engine().sync_now(graphaverage);
    ASSERT_EQ(graphaverage.get_val(), 2.5);


    
  }


  
};

