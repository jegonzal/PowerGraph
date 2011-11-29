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

class update_functor { };


class GraphTestSuite: public CxxTest::TestSuite {
public:
  

  
    
  // construct the graph at startup 
  GraphTestSuite() {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
  }
 
  void test_graph() {
    static const size_t N = 1000;
    typedef graphlab::graph<vertex_data, edge_data> graph_type;
    typedef graph_type::vertex_id_type vertex_id_type;
    typedef graph_type::vertex_id_type vertex_id_type;
    typedef graph_type::edge_id_type edge_id_type;

    graph_type g;
    vertex_data verts[N];
    edge_data edges[N];
    g.clear();
    // Add the vertices
    TS_TRACE("Checking Add Vertex");
    for(vertex_id_type i = 0; i < N; ++i) {     
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
    for(vertex_id_type i = 0; i < N; ++i) {
      edges[i].weight = i * i;
      edges[i].sum = 0;
      vertex_id_type j = (i+1) % N;
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
    for(vertex_id_type i = 0; i < N; ++i) {
      edges[i].weight = i * i;
      edges[i].sum = 0;
      vertex_id_type j = (i+1) % N;
      g.add_edge(i, j, edges[i]);
      g.add_edge(j, i, edges[i]);
    }   
  }



  void test_finalize_and_colors() {
    timer ti;
    ti.start();
    typedef graph<char, char> graph_type;
    typedef graph_type::vertex_id_type vertex_id_type;
    typedef graph_type::edge_type edge_type;
    size_t num_verts = 10000;
    size_t degree = 100;
    TS_TRACE("Constructing random graph");
    graph_type graph(num_verts);

    // create a random graph
    for(vertex_id_type i = 0; i < num_verts; ++i) {
      std::set<vertex_id_type> neighbors;
      
      for(size_t j = 0; j < degree; ++j) {
        vertex_id_type neighbor = 
          vertex_id_type(graphlab::random::uniform<vertex_id_type>(0, num_verts - 1));
        if(neighbor != i && 
           neighbors.insert(neighbor).second) graph.add_edge(i, neighbor);
      }
    }
    graph.finalize();
    std::cerr << num_verts << " * " << degree 
              << " edges created in " << ti.current_time() << " s" << std::endl;
    TS_TRACE("Testing Coloring");
    graph.compute_coloring();
    
    for(vertex_id_type i = 0; i < num_verts; ++i) {
      foreach(edge_type e, graph.in_edges(i)) {
        TS_ASSERT_DIFFERS(graph.color(e.source()), graph.color(i));
      }
      foreach(edge_type e, graph.out_edges(i)) {
        TS_ASSERT_DIFFERS(graph.color(e.target()), graph.color(i));
      }
    }
  }
                               
  void test_partition() {
    typedef graph<char, char> graph_type;
    typedef graph_type::vertex_id_type vertex_id_type;   
    // make a 100x100 grid
    size_t dim = 100;
    graph_type g(dim * dim);
    
    for (size_t i = 0;i < dim; ++i) {
      for (size_t j = 0;j < dim - 1; ++j) {
        g.add_edge(vertex_id_type(dim * i + j), 
                   vertex_id_type(dim * i + j + 1), char(1));
        g.add_edge(vertex_id_type(dim * i + j + 1), 
                   vertex_id_type(dim * i + j), char(1));
        g.add_edge(vertex_id_type(dim * j + i), 
                   vertex_id_type(dim * (j + 1) + i), char(1));
        g.add_edge(vertex_id_type(dim * (j + 1) + i), 
                   (vertex_id_type)(dim * j + i), char(1));
      }
    }
    g.finalize();    
    TS_TRACE("Random Partitioning");
    try_partition_method(g, "random");
    TS_TRACE("Metis Partitioning");
    try_partition_method(g, "metis");
    TS_TRACE("BFS Partitioning");
    try_partition_method(g, "bfs");
    TS_TRACE("Edge Number Partitioning");
    try_partition_method(g, "edge_num");
  }
  
private:
  void try_partition_method(graph<char, char> &g, 
                            const std::string& partmethod) {
    std::vector<graph_partitioner::part_id_type> parts;
    graph_partitioner::partition(partmethod, g, 4, parts);
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
