/*
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


// standard C++ headers
#include <iostream>
#include <cxxtest/TestSuite.h>

// includes the entire graphlab framework
#include <graphlab/graph/dynamic_local_graph.hpp>
#include <graphlab/macros_def.hpp>


class graph_test : public CxxTest::TestSuite {
public:

  struct vertex_data {
    size_t num_flips;
    vertex_data() : num_flips(0) { }
  };

  struct edge_data {
    int from;
    int to;
    edge_data (int f = 0, int t = 0) : from(f), to(t) {}
  };

  struct edge_data_empty { };

  typedef graphlab::dynamic_local_graph<vertex_data, edge_data> graph_type;
  typedef graph_type::edge_list_type edge_list_type;
  typedef graph_type::edge_type edge_type;
  typedef graph_type::vertex_type vertex_type;

  typedef uint32_t vertex_id_type;

  graph_type g;

  void test_sparse_graph () {
    size_t num_v = 10;
    size_t num_e = 6;

    g.finalize();

    for (size_t i = 0; i < num_v; ++i) {
      vertex_data vdata;
      g.add_vertex(vertex_id_type(i), vdata);
    }

    g.add_edge(1,3,edge_data(1,3));
    g.add_edge(2,3,edge_data(2,3));
    g.add_edge(4,3,edge_data(4,3));
    g.add_edge(5,3,edge_data(5,3));

    g.add_edge(3,2, edge_data(3,2));
    g.add_edge(3,5, edge_data(3,5));

    g.finalize();

    ASSERT_EQ(g.num_vertices(), num_v);
    ASSERT_EQ(g.num_edges(), num_e);

    for (vertex_id_type i = 0; i < 6; ++i) {
      std::cout << i << std::endl;
      edge_list_type inedges = g.in_edges(i);
      edge_list_type outedges = g.out_edges(i);
      size_t arr_insize[] = {0,0,1,4,0,1};
      size_t arr_outsize[] = {0,1,1,2,1,1};
      if (i != 3) {
        ASSERT_EQ(inedges.size(), arr_insize[i]);
        ASSERT_EQ(outedges.size(), arr_outsize[i]);
        if (outedges.size() > 0)
          {
            ASSERT_EQ(outedges[0].source().id(), i);
            ASSERT_EQ(outedges[0].target().id(), 3);

            edge_data data = (outedges[0]).data();
            ASSERT_EQ(data.from, i);
            ASSERT_EQ(data.to, 3);
          }
      } else {
        ASSERT_EQ(outedges.size(), 2);
        size_t arr_out[] = {5,2};
        std::set<size_t> set_out(arr_out, arr_out+2);
        for (size_t j = 0; j < 2; ++j) {
          edge_data data = (outedges[j]).data();
          ASSERT_EQ(data.from, 3);
          ASSERT_TRUE(set_out.find(data.to) != set_out.end());
          set_out.erase(data.to);
        }

        size_t arr_in[] = {5,4,2,1};
        std::set<size_t> set_in(arr_in, arr_in+4);
        ASSERT_EQ(inedges.size(), 4);
        for (size_t j = 0; j < 4; ++j) {
          edge_data data = (inedges[j]).data();
          ASSERT_EQ(data.to, 3);
          ASSERT_TRUE(set_in.find(data.from) != set_in.end());
          set_in.erase(data.from);
        }
      }
    }

    for (vertex_id_type i = 6; i < num_v; ++i) {
      edge_list_type inedges = g.in_edges(i);
      edge_list_type outedges = g.out_edges(i);
      ASSERT_EQ(0, inedges.size());
      ASSERT_EQ(0, outedges.size());
    }
  }

  /**
     In this function, we construct the 3 by 3 grid graph.
  */
  void test_grid_graph() {
    g.clear();
    std::cout << "-----------Begin Grid Test: --------------------" << std::endl;
    size_t dim = 3;
    size_t num_vertices = 0;
    size_t num_edge = 0;
    typedef uint32_t vertex_id_type;


    // here we create dim * dim vertices.
    for (size_t i = 0; i < dim * dim; ++i) {
      // create the vertex data, randomizing the color
      vertex_data vdata;
      vdata.num_flips = 0;
      // create the vertex
      g.add_vertex(vertex_id_type(i), vdata);
      ++num_vertices;
    }

    // create the edges. The add_edge(i,j,edgedata) function creates
    // an edge from i->j. with the edgedata attached.   edge_data edata;

    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim - 1; ++j) {
        // add the horizontal edges in both directions
        //
        g.add_edge(dim * i + j, dim * i + j + 1, edge_data(dim*i+j, dim*i+j+1));
        g.add_edge(dim * i + j + 1, dim * i + j, edge_data(dim*i+j+1, dim*i+j));

        // add the vertical edges in both directions
        g.add_edge(dim * j + i, dim * (j + 1) + i, edge_data(dim*j+i, dim*(j+1)+i));
        g.add_edge(dim * (j + 1) + i, dim * j + i, edge_data(dim*(j+1)+i, dim*j+i));
        num_edge += 4;
      }
    }

    // the graph is now constructed
    // we need to call finalize.
    g.finalize();

    printf("Test num_vertices()...\n");
    ASSERT_EQ(g.num_vertices(), num_vertices);
    printf("+ Pass test: num_vertices :)\n\n");

    printf("Test num_edges()...\n");
    ASSERT_EQ(g.num_edges(), num_edge);
    printf("+ Pass test: num_edges :)\n\n");

    // Symmetric graph: #inneighbor == outneighbor
    printf("Test num_in_neighbors() == num_out_neighbors() ...\n");
    for (size_t i = 0; i < num_vertices; ++i) {
      ASSERT_EQ(g.in_edges(i).size(), g.vertex(i).num_in_edges());
      ASSERT_EQ(g.out_edges(i).size(), g.vertex(i).num_out_edges());
      ASSERT_EQ(g.in_edges(i).size(), g.out_edges(i).size());
    }
    ASSERT_EQ(g.in_edges(4).size(), 4);
    ASSERT_EQ(g.in_edges(0).size(), 2);
    printf("+ Pass test: #in = #out...\n\n");

    printf("Test iterate over in/out_edges and get edge data. Access through local_graph: \n");
    for (vertex_id_type i = 0; i < num_vertices; ++i) {
      const edge_list_type& out_edges = g.out_edges(i);
      const edge_list_type& in_edges = g.in_edges(i);
      check_adjacency(i, in_edges, out_edges);
    }
    printf("+ Pass test: iterate edgelist and get data. :) \n");

    printf("Test iterate over in/out_edges and get edge data. Access through vertex: \n");
    for (vertex_id_type i = 0; i < num_vertices; ++i) {
      vertex_type v = g.vertex(i);
      const edge_list_type& out_edges = v.out_edges();
      const edge_list_type& in_edges = v.in_edges();
      check_adjacency(i, in_edges, out_edges);
    }
    printf("+ Pass test: iterate edgelist and get data. :) \n");
    std::cout << "-----------End Grid Test--------------------" << std::endl;
  }

  void test_dynamic_grid_graph() {
    g.clear();
    std::cout << "-----------Begin Dynamic Insertion Test: --------------------" << std::endl;
    size_t dim = 3;
    size_t num_vertices = 0;
    size_t num_edge = 0;
    typedef uint32_t vertex_id_type;

    // here we create dim * dim vertices.
    for (size_t i = 0; i < dim * dim; ++i) {
      // create the vertex data, randomizing the color
      vertex_data vdata;
      vdata.num_flips = 0;
      // create the vertex
      g.add_vertex(vertex_id_type(i), vdata);
      ++num_vertices;
    }

    g.finalize();

    // create the edges. The add_edge(i,j,edgedata) function creates
    // an edge from i->j. with the edgedata attached.   edge_data edata;

    for (size_t i = 0;i < dim; ++i) {
      for (size_t j = 0;j < dim - 1; ++j) {
        // add the horizontal edges in both directions
        g.add_edge(dim * i + j, dim * i + j + 1, edge_data(dim*i+j, dim*i+j+1));
        g.add_edge(dim * i + j + 1, dim * i + j, edge_data(dim*i+j+1, dim*i+j));
        num_edge += 2;
      }
    }

    g.finalize();

    for (size_t i = 0;i < dim; ++i) {
      for (size_t j = 0;j < dim - 1; ++j) {
        // add the vertical edges in both directions
        g.add_edge(dim * j + i, dim * (j + 1) + i, edge_data(dim*j+i, dim*(j+1)+i));
        g.add_edge(dim * (j + 1) + i, dim * j + i, edge_data(dim*(j+1)+i, dim*j+i));
        num_edge += 2;
      }
    }
    g.finalize();

    printf("Test num_vertices()...\n");
    ASSERT_EQ(g.num_vertices(), num_vertices);
    printf("+ Pass test: num_vertices :)\n\n");

    printf("Test num_edges()...\n");
    ASSERT_EQ(g.num_edges(), num_edge);
    printf("+ Pass test: num_edges :)\n\n");

    // Symmetric graph: #inneighbor == outneighbor
    printf("Test num_in_neighbors() == num_out_neighbors() ...\n");
    for (size_t i = 0; i < num_vertices; ++i) {
      ASSERT_EQ(g.in_edges(i).size(), g.vertex(i).num_in_edges());
      ASSERT_EQ(g.out_edges(i).size(), g.vertex(i).num_out_edges());
      ASSERT_EQ(g.in_edges(i).size(), g.out_edges(i).size());
    }
    ASSERT_EQ(g.in_edges(4).size(), 4);
    ASSERT_EQ(g.in_edges(0).size(), 2);
    printf("+ Pass test: #in = #out...\n\n");


    printf("Test iterate over in/out_edges and get edge data: \n");
    for (vertex_id_type i = 0; i < num_vertices; ++i) {
      const edge_list_type& out_edges = g.out_edges(i);
      const edge_list_type& in_edges = g.in_edges(i);

      printf("Test v: %u\n", i);
      printf("In edge ids: ");
      foreach(edge_type edge, in_edges)
        std::cout << "(" << edge.data().from << ","
                  << edge.data().to << ") ";
      std::cout <<std::endl;

      printf("Out edge ids: ");
      foreach(edge_type edge, out_edges)
        std::cout << "(" << edge.data().from << ","
                  << edge.data().to << ") ";
      std::cout <<std::endl;

      foreach(edge_type edge, out_edges) {
        edge_data edata = edge.data();
        ASSERT_EQ(edge.source().id(), i);
        ASSERT_EQ(edata.from, edge.source().id());
        ASSERT_EQ(edata.to, edge.target().id());
      }

      foreach(edge_type edge, in_edges) {
        edge_data edata = edge.data();
        ASSERT_EQ(edge.target().id(), i);
        ASSERT_EQ(edata.from, edge.source().id());
        ASSERT_EQ(edata.to, edge.target().id());
      }
    }
    printf("+ Pass test: iterate edgelist and get data. :) \n");

    std::cout << "-----------End Grid Test--------------------" << std::endl;
  }

private:
  void check_adjacency(size_t i,
                       const edge_list_type& in_edges,
                       const edge_list_type& out_edges) {
      std::cout << "Test v: " << i << std::endl;
      printf("In edge ids: ");
      foreach(edge_type edge, in_edges)
        std::cout << "(" << edge.data().from << ","
                  << edge.data().to << ") ";
      std::cout <<std::endl;

      printf("Out edge ids: ");
      foreach(edge_type edge, out_edges)
        std::cout << "(" << edge.data().from << ","
                  << edge.data().to << ") ";
      std::cout <<std::endl;

      foreach(edge_type edge, out_edges) {
        edge_data edata = edge.data();
        ASSERT_EQ(edge.source().id(), i);
        ASSERT_EQ(edata.from, edge.source().id());
        ASSERT_EQ(edata.to, edge.target().id());
      }

      foreach(edge_type edge, in_edges) {
        edge_data edata = edge.data();
        ASSERT_EQ(edge.target().id(), i);
        ASSERT_EQ(edata.from, edge.source().id());
        ASSERT_EQ(edata.to, edge.target().id());
      }
  }
};

#include <graphlab/macros_undef.hpp>
  // void test_add_vertex() {
  //   std::cout << "Building graph" << std::endl;
  //   graphlab::graph<std::map<int,int>, edge_data> graph;
  //   std::map<int, int> data;
  //   for(size_t i = 0; i < 50; ++i) data[i] = i;
  //   graphlab::timer time;
  //   time.start();
  //   for(vertex_id_type vid = 0; vid < 1000000; ++vid) {
  //     graph.add_vertex(vid, data);
  //   }
  //   std::cout << "add vertex runtime: " << time.current_time() << std::endl;
  // } // end of test add vertex


