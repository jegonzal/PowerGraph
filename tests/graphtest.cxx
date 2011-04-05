
// Test the graph class

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <cxxtest/TestSuite.h>


#include <graphlab.hpp>


#include <graphlab/macros_def.hpp>

using namespace graphlab;

struct vertex_data {
  size_t bias;
  size_t sum;
};

struct edge_data {
  size_t weight;
  size_t sum;
};



class GraphTestSuite: public CxxTest::TestSuite {
public:
  
  static const size_t N = 1000;
  graphlab::graph<vertex_data, edge_data> g;
  vertex_data verts[N];
  edge_data edges[N];
  
    
  // construct the graph at startup 
  GraphTestSuite() {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
  }
 
  // void test_graph() {
  //   // Add the vertices
  //   TS_TRACE("Checking Add Vertex");
  //   for(size_t i = 0; i < N; ++i) {     
  //     verts[i].bias = i; 
  //     verts[i].sum = 0;
  //     g.add_vertex(verts[i]);
  //     vertex_data& vdata = g.vertex_data(i);
  //     TS_ASSERT_EQUALS(vdata.bias, verts[i].bias);
  //     TS_ASSERT_EQUALS(vdata.sum, verts[i].sum);
  //     vdata.sum = 3;
  //     TS_ASSERT_EQUALS(vdata.sum, g.vertex_data(i).sum);
  //     TS_ASSERT_DIFFERS(vdata.sum, verts[i].sum);      
  //   }
  //   // Make a ring
  //   TS_TRACE("Checking Add Edge");
  //   for(size_t i = 0; i < N; ++i) {
  //     edges[i].weight = i * i;
  //     edges[i].sum = 0;
  //     size_t j = (i+1) % N;
  //     g.add_edge(i, j, edges[i]);
  //     g.finalize();
  //     edge_data& edata = g.edge_data(i,j);
  //     TS_ASSERT_EQUALS(edata.weight, i * i);
  //     TS_ASSERT_EQUALS(edata.sum, (size_t)0);
  //     edata.sum = 3;
  //     TS_ASSERT_EQUALS(edata.sum, g.edge_data(i,j).sum);
  //     TS_ASSERT(edata.sum != edges[i].sum);
      
  //   }   

  //   TS_TRACE("Checking Num vertices");
  //   TS_ASSERT_EQUALS(g.num_vertices(),  N);

  //   // Make a ring
  //   TS_TRACE("Checking Add Edge Again with bi directed edges");
  //   for(size_t i = 0; i < N; ++i) {
  //     edges[i].weight = i * i;
  //     edges[i].sum = 0;
  //     size_t j = (i+1) % N;
  //     g.add_edge(i, j, edges[i]);
  //     g.add_edge(j, i, edges[i]);
  //   }   
    

  // }

  // void test_vertex_colors() {
  //   for(size_t i = 0; i < N; ++i) {     
  //     g.set_vertex_color(i, i%10); 
  //   }
  //   for(size_t i = 0; i < N; ++i) {     
  //     TS_ASSERT_EQUALS(g.vertex_color(i), i%10);
  //   }          
  // }
  
                               
                               

  void test_finalize() {
    namespace random = graphlab::random;

    typedef graph<char, char> graph_type;
    typedef types< graph_type > gl;
    size_t num_verts = 100000;
    size_t degree = 1000;
    
    std::cout << "Creating graph" << std::endl;
    gl::graph graph(num_verts);

    { timer time;
      std::cout << "Adding random edges" << std::endl;
      time.start();
      for(size_t i = 0; i < num_verts; ++i) {
        if(i % (num_verts / 100) == 0) {
          std::cout << "\b\b\b\b\b" << size_t(100.0 * i / num_verts) << "%";
          std::cout.flush();
        }

        std::set<gl::vertex_id_t> neighbors;
        
        for(size_t j = 0; j < degree; ++j) {
          size_t neighbor = random::uniform<size_t>(0, num_verts - 1);
          if(neighbors.insert(neighbor).second)
            graph.add_edge(i, neighbor);
        }
      }
      std::cout << "Finished in " << time.current_time()
                << " seconds." << std::endl;
    }
    
    { timer time;
      std::cout << "Finalizing graph" << std::endl;
      time.start();
      graph.finalize();
      std::cout << "Finished in " << time.current_time()
                << " seconds." << std::endl;
    }
  }
};



#include <graphlab/macros_undef.hpp>
