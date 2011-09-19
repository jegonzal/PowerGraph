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
#include <graphlab/graph/disk_graph.hpp>

#include <graphlab/macros_def.hpp>

using namespace graphlab;

struct vertex_data {
  size_t bias;
  size_t sum;
};
SERIALIZABLE_POD(vertex_data);

struct edge_data {
  size_t weight;
  size_t sum;
};
SERIALIZABLE_POD(edge_data);

typedef graph<vertex_data, edge_data> graph_type;
//typedef graph_type::vertex_id_type vertex_id_t;

vertex_data vertexdata_generator(vertex_id_t idx) {
  return vertex_data();
}

template <size_t degree, size_t num_verts>
void edge_generator(vertex_id_t vid,
                    const vertex_data& vdata,
                    std::vector<vertex_id_t>& inv, 
                    std::vector<edge_data>& inedata,
                    std::vector<vertex_id_t>& outv, 
                    std::vector<edge_data>& outedata) {
  std::set<vertex_id_t> neighbors;
  for(size_t j = 0; j < degree; ++j) {
    vertex_id_t neighbor = (vertex_id_t)graphlab::random::uniform<size_t>(0, num_verts - 1);
    if(neighbor != vid && neighbors.insert(neighbor).second) {
      outv.push_back(neighbor);
      outedata.push_back(edge_data());
    }
  }
  
}

class GraphTestSuite: public CxxTest::TestSuite {
public:
  

  
    
  // construct the graph at startup 
  GraphTestSuite() {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
  }
 
  void test_diskgraph_construction() {
    static const size_t N = 1000;
    graphlab::disk_graph<vertex_data, edge_data> g("dg0", 10);
    vertex_data verts[N];
    edge_data edges[N];
    g.clear();
    // Add the vertices
    TS_TRACE("Checking Add Vertex");
    for(vertex_id_t i = 0; i < N; ++i) {     
      verts[i].bias = i; 
      verts[i].sum = 0;
      g.add_vertex(verts[i]);
      vertex_data vdata = g.get_vertex_data(i);
      TS_ASSERT_EQUALS(vdata.bias, verts[i].bias);
      TS_ASSERT_EQUALS(vdata.sum, verts[i].sum);
      vdata.sum = 3;
      g.set_vertex_data(i, vdata);
      TS_ASSERT_EQUALS(vdata.sum, g.get_vertex_data(i).sum);
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
      edge_data edata = g.get_edge_data(i,j); 
      TS_ASSERT_EQUALS(edata.weight, i * i);
      TS_ASSERT_EQUALS(edata.sum, (size_t)0);
      edata.sum = 3;
      g.set_edge_data(i,j,edata);
      TS_ASSERT_EQUALS(edata.sum, g.get_edge_data(i,j).sum);
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
    g.finalize();
  }

  void test_diskgraph_structure() {
    const size_t num_verts = 10000;
    const size_t degree = 10;
    graphlab::graph<vertex_data, edge_data> memgraph;
    
    {
      TS_TRACE("Constructing random graph");
            
      graphlab::disk_graph<vertex_data, edge_data> graph("dg1", 10);
      graph.clear();

      graph.add_vertex_collection("coll", num_verts, vertexdata_generator);
      for(vertex_id_t i = 0; i < num_verts; ++i) memgraph.add_vertex(vertex_data());
      
      // create a random graph
      for(vertex_id_t i = 0; i < num_verts; ++i) {
        std::set<vertex_id_t> neighbors;
        
        for(size_t j = 0; j < degree; ++j) {
          vertex_id_t neighbor = (vertex_id_t)graphlab::random::uniform<size_t>(0, num_verts - 1);
          if(neighbor != i && neighbors.insert(neighbor).second) {
            edge_data ed;
            ed.weight = i;
            ed.sum = neighbor;
            memgraph.add_edge(i, neighbor, ed);
            graph.add_edge(i, neighbor, ed);
          }
        }
      }
      TS_ASSERT_EQUALS(graph.num_vertices(), memgraph.num_vertices());
      TS_ASSERT_EQUALS(graph.num_edges(), memgraph.num_edges());
      graph.finalize();
    }
    // check the graph
    {
      TS_TRACE("Reloading graph");
      graphlab::disk_graph<vertex_data, edge_data> graph("dg1.idx");
      TS_ASSERT_EQUALS(graph.num_vertices(), memgraph.num_vertices());
      TS_ASSERT_EQUALS(graph.num_edges(), memgraph.num_edges());
      TS_TRACE("Checking graph");
      for(vertex_id_t i = 0; i < num_verts; ++i) {
        // get the outvertices for each vertex
        std::vector<vertex_id_t> outv = graph.out_vertices(i);
        std::vector<vertex_id_t> outvmem = memgraph.out_vertices(i);
        // test if they are the same size
        TS_ASSERT_EQUALS(outv.size(), outvmem.size());
        std::sort(outv.begin(), outv.end());
        std::sort(outvmem.begin(), outvmem.end());
        // compare if the vector contents are identical
        for (size_t j = 0;j < outv.size(); ++j) {
          TS_ASSERT_EQUALS(outv[j], outvmem[j]);
          // compare if the edge data is identical
          edge_data ed1 = graph.get_edge_data(i, outv[j]);
          edge_data ed2 = memgraph.edge_data(i, outvmem[j]);
          TS_ASSERT_EQUALS(ed1.weight, ed2.weight);
          TS_ASSERT_EQUALS(ed1.sum, ed2.sum);
        }
        
        // repeat for in vertices
        std::vector<vertex_id_t> inv = graph.in_vertices(i);
        std::vector<vertex_id_t> invmem = memgraph.in_vertices(i);
        // test if they are the same size
        TS_ASSERT_EQUALS(inv.size(), invmem.size());
        std::sort(inv.begin(), inv.end());
        std::sort(invmem.begin(), invmem.end());
        // compare if the vector contents are identical
        for (size_t j = 0;j < inv.size(); ++j) {
          TS_ASSERT_EQUALS(inv[j], invmem[j]);
        }
      }
    }
  }

  
  
   void test_graph_perf2() {
    const size_t num_verts = 10000;
    const size_t degree = 100;
    TS_TRACE("Constructing random graph");
    graphlab::disk_graph<vertex_data, edge_data> graph("dg2", 10);
    graph.clear();
    timer ti;
    ti.start();
    graph.add_vertex_collection("coll", num_verts, vertexdata_generator);
    // create a random graph
    graph.add_edge_indirect("coll", edge_generator<degree, num_verts>);
    std::cerr << "add_edge_collection: " << num_verts << " * " << degree << " edges created in " << ti.current_time() << " s" << std::endl;
    
    TS_TRACE("Testing Coloring");
    ti.start();
    graph.compute_coloring();
    std::cout << "Coloring took " << ti.current_time() << " s" << std::endl;
    TS_TRACE("Verifying Coloring");
    for(vertex_id_t i = 0; i < num_verts; ++i) {
      size_t colori = graph.get_color(i);
      foreach(vertex_id_t nbr, graph.in_vertices(i)) {
        TS_ASSERT_DIFFERS(graph.get_color(nbr),colori);
      }
      foreach(vertex_id_t nbr, graph.out_vertices(i)) {
        TS_ASSERT_DIFFERS(graph.get_color(nbr),colori);
      }
    }
    graph.finalize();
  }
  
  
  
  
  void test_diskgraph_from_memgraph() {
    const size_t num_verts = 10000;
    const size_t degree = 100;
    graphlab::graph<vertex_data, edge_data> memgraph;
    
    TS_TRACE("Constructing disk graph from mem graph");
          
    for(vertex_id_t i = 0; i < num_verts; ++i) memgraph.add_vertex(vertex_data());
    
    // create a random graph
    for(vertex_id_t i = 0; i < num_verts; ++i) {
      std::set<vertex_id_t> neighbors;
      
      for(size_t j = 0; j < degree; ++j) {
        vertex_id_t neighbor = (vertex_id_t)graphlab::random::uniform<size_t>(0, num_verts - 1);
        if(neighbor != i && neighbors.insert(neighbor).second) {
          edge_data ed;
          ed.weight = i;
          ed.sum = neighbor;
          memgraph.add_edge(i, neighbor, ed);
        }
      }
    }
    {
      timer ti;
      ti.start();
      
      graphlab::disk_graph<vertex_data, edge_data> graph("dg3", 10);
      graph = memgraph;
      graph.finalize();
            
      std::cerr << "disk from mem: " << num_verts << " * " << degree << " edges created in " << ti.current_time() << " s" << std::endl;
      
      TS_TRACE("Checking graph");
      TS_ASSERT_EQUALS(graph.num_vertices(), memgraph.num_vertices());
      TS_ASSERT_EQUALS(graph.num_edges(), memgraph.num_edges());
      for(vertex_id_t i = 0; i < num_verts; ++i) {
        // get the outvertices for each vertex
        std::vector<vertex_id_t> outv = graph.out_vertices(i);
        std::vector<vertex_id_t> outvmem = memgraph.out_vertices(i);
        // test if they are the same size
        TS_ASSERT_EQUALS(outv.size(), outvmem.size());
        std::sort(outv.begin(), outv.end());
        std::sort(outvmem.begin(), outvmem.end());
        // compare if the vector contents are identical
        for (size_t j = 0;j < outv.size(); ++j) {
          TS_ASSERT_EQUALS(outv[j], outvmem[j]);
          // compare if the edge data is identical
          edge_data ed1 = graph.get_edge_data(i, outv[j]);
          edge_data ed2 = memgraph.edge_data(i, outvmem[j]);
          TS_ASSERT_EQUALS(ed1.weight, ed2.weight);
          TS_ASSERT_EQUALS(ed1.sum, ed2.sum);
        }
        
        // repeat for in vertices
        std::vector<vertex_id_t> inv = graph.in_vertices(i);
        std::vector<vertex_id_t> invmem = memgraph.in_vertices(i);
        // test if they are the same size
        TS_ASSERT_EQUALS(inv.size(), invmem.size());
        std::sort(inv.begin(), inv.end());
        std::sort(invmem.begin(), invmem.end());
        // compare if the vector contents are identical
        for (size_t j = 0;j < inv.size(); ++j) {
          TS_ASSERT_EQUALS(inv[j], invmem[j]);
        }
      }
      TS_TRACE("Making Constant Atoms");
      graph.make_const_atoms();
    }   
    
    {
      timer ti;
      ti.start();
      
      graphlab::disk_graph<vertex_data, edge_data> graph("dg3", 10, true);
      
      TS_TRACE("Checking constant graph");
      TS_ASSERT_EQUALS(graph.num_vertices(), memgraph.num_vertices());
      TS_ASSERT_EQUALS(graph.num_edges(), memgraph.num_edges());
      for(vertex_id_t i = 0; i < num_verts; ++i) {
        // get the outvertices for each vertex
        std::vector<vertex_id_t> outv = graph.out_vertices(i);
        std::vector<vertex_id_t> outvmem = memgraph.out_vertices(i);
        // test if they are the same size
        TS_ASSERT_EQUALS(outv.size(), outvmem.size());
        std::sort(outv.begin(), outv.end());
        std::sort(outvmem.begin(), outvmem.end());
        // compare if the vector contents are identical
        for (size_t j = 0;j < outv.size(); ++j) {
          TS_ASSERT_EQUALS(outv[j], outvmem[j]);
          // compare if the edge data is identical
          edge_data ed1 = graph.get_edge_data(i, outv[j]);
          edge_data ed2 = memgraph.edge_data(i, outvmem[j]);
          TS_ASSERT_EQUALS(ed1.weight, ed2.weight);
          TS_ASSERT_EQUALS(ed1.sum, ed2.sum);
        }
        
        // repeat for in vertices
        std::vector<vertex_id_t> inv = graph.in_vertices(i);
        std::vector<vertex_id_t> invmem = memgraph.in_vertices(i);
        // test if they are the same size
        TS_ASSERT_EQUALS(inv.size(), invmem.size());
        std::sort(inv.begin(), inv.end());
        std::sort(invmem.begin(), invmem.end());
        // compare if the vector contents are identical
        for (size_t j = 0;j < inv.size(); ++j) {
          TS_ASSERT_EQUALS(inv[j], invmem[j]);
        }
      }
    }
  }

};



#include <graphlab/macros_undef.hpp>
