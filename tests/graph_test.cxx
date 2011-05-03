
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
  

  
    
  // construct the graph at startup 
  GraphTestSuite() {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
  }
 
  void test_graph() {
    static const size_t N = 1000;
    graphlab::graph<vertex_data, edge_data> g;
    vertex_data verts[N];
    edge_data edges[N];
    g.clear();
    // Add the vertices
    TS_TRACE("Checking Add Vertex");
    for(vertex_id_t i = 0; i < N; ++i) {     
      verts[i].bias = i; 
      verts[i].sum = 0;
      g.add_vertex(verts[i]);
      vertex_data& vdata = g.vertex_data(i);
      TS_ASSERT_EQUALS(vdata.bias, verts[i].bias);
      TS_ASSERT_EQUALS(vdata.sum, verts[i].sum);
      vdata.sum = 3;
      TS_ASSERT_EQUALS(vdata.sum, g.vertex_data(i).sum);
      TS_ASSERT_DIFFERS(vdata.sum, verts[i].sum);      
    }
    // Make a ring
    TS_TRACE("Checking Add Edge");
    for(vertex_id_t i = 0; i < N; ++i) {
      edges[i].weight = i * i;
      edges[i].sum = 0;
      vertex_id_t j = (i+1) % N;
      g.add_edge(i, j, edges[i]);
      g.finalize();
      edge_data& edata = g.edge_data(i,j);
      TS_ASSERT_EQUALS(edata.weight, i * i);
      TS_ASSERT_EQUALS(edata.sum, (size_t)0);
      edata.sum = 3;
      TS_ASSERT_EQUALS(edata.sum, g.edge_data(i,j).sum);
      TS_ASSERT(edata.sum != edges[i].sum);
      
    }   

    TS_TRACE("Checking Num vertices");
    TS_ASSERT_EQUALS(g.num_vertices(),  N);

    // Make a ring
    TS_TRACE("Checking Add Edge Again with bi directed edges");
    for(vertex_id_t i = 0; i < N; ++i) {
      edges[i].weight = i * i;
      edges[i].sum = 0;
      vertex_id_t j = (i+1) % N;
      g.add_edge(i, j, edges[i]);
      g.add_edge(j, i, edges[i]);
    }   
  }



  void test_finalize_and_colors() {
    typedef graph<char, char> graph_type;
    typedef types< graph_type > gl;
    size_t num_verts = 10000;
    size_t degree = 100;
    TS_TRACE("Constructing random graph");
    gl::graph graph(num_verts);
    // create a random graph
    for(vertex_id_t i = 0; i < num_verts; ++i) {
      std::set<gl::vertex_id_t> neighbors;
      
      for(size_t j = 0; j < degree; ++j) {
        vertex_id_t neighbor = (vertex_id_t)graphlab::random::uniform<size_t>(0, num_verts - 1);
        if(neighbor != i && neighbors.insert(neighbor).second) graph.add_edge(i, neighbor);
      }
    }
    graph.finalize();
    TS_TRACE("Testing Coloring");
    graph.compute_coloring();
    
    for(vertex_id_t i = 0; i < num_verts; ++i) {
      foreach(edge_id_t e, graph.in_edge_ids(i)) {
        TS_ASSERT_DIFFERS(graph.color(graph.source(e)), graph.color(i));
      }
      foreach(edge_id_t e, graph.out_edge_ids(i)) {
        TS_ASSERT_DIFFERS(graph.color(graph.target(e)), graph.color(i));
      }
    }
  }
                               
  void test_partition() {
    typedef graph<char, char> graph_type;
    typedef types< graph_type > gl;
    // make a 100x100 grid
    size_t dim = 100;
    gl::graph g(dim * dim);
    
    for (size_t i = 0;i < dim; ++i) {
      for (size_t j = 0;j < dim - 1; ++j) {
        g.add_edge((vertex_id_t)(dim * i + j), (vertex_id_t)(dim * i + j + 1), char(1));
        g.add_edge((vertex_id_t)(dim * i + j + 1), (vertex_id_t)(dim * i + j), char(1));
        g.add_edge((vertex_id_t)(dim * j + i), (vertex_id_t)(dim * (j + 1) + i), char(1));
        g.add_edge((vertex_id_t)(dim * (j + 1) + i), (vertex_id_t)(dim * j + i), char(1));
      }
    }
    g.finalize();    
    TS_TRACE("Random Partitioning");
    try_partition_method(g, partition_method::PARTITION_RANDOM);
    TS_TRACE("Metis Partitioning");
    try_partition_method(g, partition_method::PARTITION_METIS);
    TS_TRACE("BFS Partitioning");
    try_partition_method(g, partition_method::PARTITION_BFS);
    TS_TRACE("Edge Number Partitioning");
    try_partition_method(g, partition_method::PARTITION_EDGE_NUM);
  }
  
private:
  void try_partition_method(graph<char, char> &g, 
                            partition_method::partition_method_enum partmethod) {
    std::vector<vertex_id_t> parts;
    g.partition(partmethod, 4, parts);
    TS_ASSERT_EQUALS(parts.size(), g.num_vertices());  

    std::vector<size_t> vcount(4);
    for (size_t i = 0;i < parts.size(); ++i) {
      TS_ASSERT(parts[i] >= 0);
      TS_ASSERT(parts[i] < 4);
      vcount[parts[i]]++;
    }
    //make sure we actually do have 4 partitions
    for (size_t i = 0;i < vcount.size(); ++i) {
      TS_ASSERT(vcount[i] > 0);
    }
  }
};



#include <graphlab/macros_undef.hpp>
