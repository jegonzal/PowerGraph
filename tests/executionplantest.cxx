
// Test the graph class

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <cxxtest/TestSuite.h>

#include <graphlab.hpp>


#include <graphlab/macros_def.hpp>

using namespace graphlab;


typedef double vertex_type;
typedef double edge_type;

typedef graphlab::graph<vertex_type, edge_type> graph_type;
typedef graphlab::types<graph_type> gl_types;


bool is_even(vertex_id_t v, const vertex_type& vdata) {
  return (v % 2 == 0);
}

class ExecutionPlanTest: public CxxTest::TestSuite {
public:
  
  static const size_t N = 10;
  graph_type g;
  
    
  // construct the graph at startup 
  ExecutionPlanTest() {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
  }
 
  void test_graph() {
    // Add the vertices
    TS_TRACE("Checking Add Vertex");
    ss_set_type allv;
    for(size_t i = 0; i < N; ++i) {     
      g.add_vertex();
      ss_insert(allv, i);
    }
    // Make a ring
    TS_TRACE("Checking Add Edge");
    for(size_t i = 0; i < N; ++i) {
      size_t j = (i+1) % N;
      g.add_edge(i, j);
    }   
    
    gl_types::vset v1;
    gl_types::rvset v2(is_even);
    v1.init(&g, NULL, NULL, 1);
    v2.init(&g, &v1, NULL, 1);

    v1.resolve_event_handlers();
    // set v1 to everything
    v1.rebuild(NULL, allv);
    std::cout << "v1 size = " << v1.size() << "\n";
    std::cout << "v2 size = " << v2.size() << "\n";
    gl_types::execution_plan eplan;
    eplan.execute(v1, NULL);
    eplan.execute(v2, NULL);
    eplan.generate_plan(g, 2);
    std::cout <<" PLAN 0 ============================ \n";
    eplan.print_plan(0);
    std::cout <<" PLAN 1 ============================ \n";
    eplan.print_plan(1);
    
  }

};



#include <graphlab/macros_undef.hpp>
