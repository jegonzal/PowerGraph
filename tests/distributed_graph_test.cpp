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

#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/graph_vertex_join.hpp>
#include <graphlab/macros_def.hpp>



struct vertex_data: public graphlab::IS_POD_TYPE  {
  size_t i;
  vertex_data() : i(0) { }
};

struct edge_data: public graphlab::IS_POD_TYPE  {
  int from;
  int to;
  edge_data (int f = 0, int t = 0) : from(f), to(t) {}
};

struct edge_data_empty { };

typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graph_type::local_edge_list_type local_edge_list_type;
typedef graph_type::local_edge_type local_edge_type;
typedef graph_type::local_vertex_type local_vertex_type;


//graph 2 used for testing joins
struct vertex_data2: public graphlab::IS_POD_TYPE  {
  size_t i;
  size_t j;
  vertex_data2() : i(0), j((size_t)(-1)) { }
};

typedef graphlab::distributed_graph<vertex_data2, edge_data> graph_type2;

template <typename VType>
size_t get_vid(const VType& t) {
  return t.id();
}
template <typename VType>
size_t get_vid_half(const VType& t) {
  return t.id() / 2;
}

template <typename VType, typename U>
void assign_data_i_to_i(VType& vtype, const U& u) {
  vtype.data().i = u.i;
}

template <typename VType, typename U>
void assign_data_j_to_i(VType& vtype, const U& u) {
  vtype.data().i = u.j;
}


int main(int argc, char** argv) {
  graphlab::mpi_tools::init(argc, argv);
  global_logger().set_log_level(LOG_INFO);

  graphlab::distributed_control dc;
  graph_type g(dc);
  graph_type2 g2(dc);

  std::cout << "-----------Begin Grid Test: Object Accessors--------------------" << std::endl;
  size_t dim = 3;
  size_t num_vertices = 0;
  size_t num_edge = 0;
  typedef uint32_t vertex_id_type;

  if (dc.procid() == 0) {
    // here we create dim * dim vertices.
    for (size_t i = 0; i < dim * dim; ++i) {
      // create the vertex data, randomizing the color
      vertex_data vdata;
      vdata.i = i;
      // create the vertex
      g.add_vertex(vertex_id_type(i), vdata);
      g2.add_vertex(vertex_id_type(2 * i));
      ++num_vertices;
    }

    // create the edges. The add_edge(i,j,edgedata) function creates
    // an edge from i->j. with the edgedata attached.   edge_data edata;

    for (size_t i = 0;i < dim; ++i) {
      for (size_t j = 0;j < dim - 1; ++j) {
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
  }
  dc.all_reduce(num_vertices);
  dc.all_reduce(num_edge);
  // the graph is now constructed
  // we need to call finalize.
  g.finalize();
  g2.finalize();

  ASSERT_EQ(g.num_vertices(), g2.num_vertices());

  dc.cout() << "Test num_vertices()...\n";
  ASSERT_EQ(g.num_vertices(), num_vertices);
  dc.cout() << "+ Pass test: num_vertices :)\n\n";

  dc.cout() << "Test num_edges()...\n";
  ASSERT_EQ(g.num_edges(), num_edge);
  dc.cout() << "+ Pass test: num_edges :)\n\n";

  // Symmetric graph: #inneighbor == outneighbor
  dc.cout() << "Test num_in_neighbors() == num_out_neighbors() ...\n";
  for (size_t i = 0; i < g.num_local_vertices(); ++i) {
    graph_type::vertex_type v(g.l_vertex(i));
    ASSERT_EQ(v.num_in_edges(), v.num_out_edges());
    if (v.id() == 4) {
      ASSERT_EQ(v.num_in_edges(), 4);
    }
    if (v.id() == 0) {
      ASSERT_EQ(v.num_in_edges(), 2);
    }
  }
  dc.cout() << "+ Pass test: #in = #out...\n\n";


  dc.cout() << "Test iterate over in/out_edges and get edge data: \n";
  for (graphlab::lvid_type i = 0; i < g.num_local_vertices(); ++i) {
    local_vertex_type v = local_vertex_type(g.l_vertex(i));
    const local_edge_list_type& out_edges = v.out_edges();
    const local_edge_list_type& in_edges = v.in_edges();

    printf("Test v: %u\n", v.global_id());
    printf("In edge ids: ");
    foreach(local_edge_type edge, in_edges)
      std::cout << "(" << edge.data().from << ","
                << edge.data().to << ") ";
    std::cout <<std::endl;

    printf("Out edge ids: ");
    foreach(local_edge_type edge, out_edges)
      std::cout << "(" << edge.data().from << ","
                << edge.data().to << ") ";
    std::cout <<std::endl;

    foreach(local_edge_type edge, out_edges) {
      edge_data edata = edge.data();
      ASSERT_EQ(edata.from, edge.source().global_id());
      ASSERT_EQ(edata.to, edge.target().global_id());
    }

    foreach(local_edge_type edge, in_edges) {
      edge_data edata = edge.data();
      ASSERT_EQ(edata.from, edge.source().global_id());
      ASSERT_EQ(edata.to, edge.target().global_id());
    }
  }
  dc.cout() << "+ Pass test: iterate edgelist and get data. :) \n";
  std::cout << "-----------End Grid Test--------------------" << std::endl;

  
  dc.cout() << "Testing Injective join\n";
  graphlab::graph_vertex_join<graph_type, graph_type2> join(dc, g, g2);
  join.prepare_injective_join(get_vid<graph_type::vertex_type>,
                              get_vid_half<graph_type2::vertex_type>);
  join.right_injective_join(assign_data_i_to_i<graph_type2::vertex_type,
                                               graph_type::vertex_data_type>);
  // check g2
  for (graphlab::lvid_type i = 0; i < g2.num_local_vertices(); ++i) {
    ASSERT_EQ(g2.l_vertex(i).data().i, g2.l_vertex(i).global_id() / 2);
  }

  // now fold back this time assigning j from graph2 to i from graph 1
  // this should make all i's in graph 1 take the value -1
  join.left_injective_join(assign_data_j_to_i<graph_type::vertex_type,
                                              graph_type2::vertex_data_type>);
  for (graphlab::lvid_type i = 0; i < g.num_local_vertices(); ++i) {
    ASSERT_EQ(g.l_vertex(i).data().i, (size_t)(-1));
  }
dc.barrier();
  dc.cout() << "Injective join pass\n";
  graphlab::mpi_tools::finalize();
}

#include <graphlab/macros_undef.hpp>
