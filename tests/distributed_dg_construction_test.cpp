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


#include <graphlab/graph/disk_graph.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/graph/mr_disk_graph_construction.hpp>

using namespace graphlab;

/*
 * Simple example of a graph constructor. This graph constructor constructs
 * a distributed graph from an in-memory graph. Of course, this is quite 
 * unecessary as a disk graph can be easily constructed from an im memory graph directly,
 * but this provides a simple test case / example for using the distributed graph construction
 * scheme
 */
class graph_constructor: public igraph_constructor<float, double>{

public: 
  typedef igraph_constructor<float, double> base_type;
  typedef base_type::vertex_id_type vertex_id_type;

private:

  graph<float, double>& mg;
  size_t i;
  size_t max;
  
  size_t viterator;
  size_t eiterator;
public:
  graph_constructor(graph<float, double>& mg):mg(mg) { }
  graph_constructor(graph_constructor& g):mg(g.mg) { }
  
  /**
     Begin iteration. 
     Instance "i" will store vertex IDs i, i+N, i+2N, etc.
     and edge IDs i, i+N, i+2N, etc.
  */
  iteration_method begin(size_t i_, size_t max_) {
    std::cout << "begin: " << i_ << " " << max_ << std::endl;
    i = i_;
    max = max_;
    viterator = i;
    eiterator = i;
    return ExternalIteration;
  }
  
  /**
     We use a simple modular atom partitioning
  */
  uint16_t vertex_to_atomid(vertex_id_type vid, uint16_t numatoms)  {
    return vid % numatoms;
  }
  
  /**
     This function will be called repeatedly until NoMoreData is returned.
  */
  iterate_return_type iterate(vertex_id_type& vtx, 
                              float& vdata,
                              uint32_t& color,
                              std::pair<vertex_id_type, vertex_id_type>& edge, 
                              double& edata) {
    /// while I still have vertices to insert, return a vertex
    while (viterator < mg.num_vertices()) {
      vtx = viterator;
      vdata = mg.vertex_data(viterator);
      color = mg.color(viterator);
      viterator += max;
      return Vertex;
    }

    /// while I still have edges to insert, return an edge
    while (eiterator < mg.num_edges()) {
      edge.first = mg.source(eiterator);
      edge.second = mg.target(eiterator);
      edata = mg.edge_data(eiterator);
      eiterator += max;
      return Edge;
    }
    return NoMoreData;
  }
};

  
int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  mpi_tools::init(argc, argv);
  
  
  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  // create distributed control
  distributed_control dc(param);
  
  
  const size_t num_verts = 10000;
  const size_t degree = 100;
  
  // here we first create an in memory graph. 
  // Every machine creates exactly the same graph.
  
  srand(42);

  std::cout << "Testing Distributed Disk Graph Construction" << std::endl;
  std::cout << "Creating a graph" << std::endl;
  graph<float, double> memgraph;
  typedef graph<float, double>::vertex_id_type vertex_id_type;

  for(vertex_id_type i = 0; i < num_verts; ++i) memgraph.add_vertex(rand());
  for(vertex_id_type i = 0; i < num_verts; ++i) {
    for(size_t j = 0; j < degree; ++j) memgraph.add_edge(i, (i + j + 1) % num_verts, rand());
  }
  memgraph.compute_coloring();
  //-----------------------------------------------------------------------
  // create an instance of a graph cosntructor
  graph_constructor gc(memgraph);
  std::cout << "Starting MR construction..." << std::endl;
  timer ti;
  ti.start();
  // call the mr_disk_graph_construction function
  mr_disk_graph_construction<graph_constructor,float, double>(dc, gc, 2, "dg", 16, 
                                                              disk_graph_atom_type::WRITE_ONLY_ATOM, "/tmp", "./");
  // thats all!
  std::cout << "Completed in " << ti.current_time() << " s" << std::endl;
  //-----------------------------------------------------------------------
  
  // everything after this is just for verification
  if (dc.procid() == 0) {
    std::cout << "Checking constructed graph ... " << std::endl;
    disk_graph<float,double> graph(disk_graph_atom_type::MEMORY_ATOM, "dg.idx");
    ASSERT_EQ(graph.num_vertices(), memgraph.num_vertices());
    ASSERT_EQ(graph.num_edges(), memgraph.num_edges());
    for(vertex_id_type i = 0; i < num_verts; ++i) {
      // get the outvertices for each vertex
      std::vector<vertex_id_type> outv = graph.out_vertices(i);
      std::vector<vertex_id_type> outvmem = memgraph.out_vertices(i);
      // test if they are the same size
      ASSERT_EQ(outv.size(), outvmem.size());
      std::sort(outv.begin(), outv.end());
      std::sort(outvmem.begin(), outvmem.end());
      // compare if the vector contents are identical
      float v1 = graph.get_vertex_data(i);
      float v2 = memgraph.vertex_data(i);
      ASSERT_EQ(v1, v2);
      for (size_t j = 0;j < outv.size(); ++j) {
        ASSERT_EQ(outv[j], outvmem[j]);
        // compare if the edge data is identical
        double ed1 = graph.get_edge_data(i, outv[j]);
        double ed2 = memgraph.edge_data(i, outvmem[j]);
        ASSERT_EQ(ed1, ed2);
      }
      
      // repeat for in vertices
      std::vector<vertex_id_type> inv = graph.in_vertices(i);
      std::vector<vertex_id_type> invmem = memgraph.in_vertices(i);
      // test if they are the same size
      ASSERT_EQ(inv.size(), invmem.size());
      std::sort(inv.begin(), inv.end());
      std::sort(invmem.begin(), invmem.end());
      // compare if the vector contents are identical
      for (size_t j = 0;j < inv.size(); ++j) {
        ASSERT_EQ(inv[j], invmem[j]);
      }
    }
    graph.make_memory_atoms();
  }
  
  if (dc.procid() == 0) {
    std::cout << "Checking constructed fast graph ... " << std::endl;
    disk_graph<float,double> graph(disk_graph_atom_type::MEMORY_ATOM, "dg.idx");

    ASSERT_EQ(graph.num_vertices(), memgraph.num_vertices());
    ASSERT_EQ(graph.num_edges(), memgraph.num_edges());
    for(vertex_id_type i = 0; i < num_verts; ++i) {
      // get the outvertices for each vertex
      std::vector<vertex_id_type> outv = graph.out_vertices(i);
      std::vector<vertex_id_type> outvmem = memgraph.out_vertices(i);
      // test if they are the same size
      ASSERT_EQ(outv.size(), outvmem.size());
      std::sort(outv.begin(), outv.end());
      std::sort(outvmem.begin(), outvmem.end());
      // compare if the vector contents are identical
      float v1 = graph.get_vertex_data(i);
      float v2 = memgraph.vertex_data(i);
      ASSERT_EQ(v1, v2);
      for (size_t j = 0;j < outv.size(); ++j) {
        ASSERT_EQ(outv[j], outvmem[j]);
        // compare if the edge data is identical
        double ed1 = graph.get_edge_data(i, outv[j]);
        double ed2 = memgraph.edge_data(i, outvmem[j]);
        ASSERT_EQ(ed1, ed2);
      }
      
      // repeat for in vertices
      std::vector<vertex_id_type> inv = graph.in_vertices(i);
      std::vector<vertex_id_type> invmem = memgraph.in_vertices(i);
      // test if they are the same size
      ASSERT_EQ(inv.size(), invmem.size());
      std::sort(inv.begin(), inv.end());
      std::sort(invmem.begin(), invmem.end());
      // compare if the vector contents are identical
      for (size_t j = 0;j < inv.size(); ++j) {
        ASSERT_EQ(inv[j], invmem[j]);
      }
    }
  }
  mpi_tools::finalize();
}
