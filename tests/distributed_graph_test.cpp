// standard C++ headers
#include <iostream>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/macros_def.hpp>



struct vertex_data: public graphlab::IS_POD_TYPE  {
  size_t num_flips;
  vertex_data() : num_flips(0) { }
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



int main(int argc, char** argv) {
  graphlab::mpi_tools::init(argc, argv);
  global_logger().set_log_level(LOG_INFO);

  graphlab::distributed_control dc;
  graph_type g(dc);

  std::cout << "-----------Begin Grid Test: Object Accessors--------------------" << std::endl;
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
    ASSERT_EQ(g.vertex(i).num_in_edges(), g.vertex(i).num_out_edges());
  }
  ASSERT_EQ(g.vertex(4).num_in_edges(), 4);
  ASSERT_EQ(g.vertex(0).num_in_edges(), 2);
  printf("+ Pass test: #in = #out...\n\n");


  printf("Test iterate over in/out_edges and get edge data: \n");
  for (vertex_id_type i = 0; i < num_vertices; ++i) {
    local_vertex_type v = local_vertex_type(g.vertex(i));
    const local_edge_list_type& out_edges = v.out_edges();
    const local_edge_list_type& in_edges = v.in_edges();

    printf("Test v: %u\n", i);
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
      ASSERT_EQ(edge.source().global_id(), i);
      ASSERT_EQ(edata.from, edge.source().global_id());
      ASSERT_EQ(edata.to, edge.target().global_id());
    }

    foreach(local_edge_type edge, in_edges) {
      edge_data edata = edge.data();
      ASSERT_EQ(edge.target().global_id(), i);
      ASSERT_EQ(edata.from, edge.source().global_id());
      ASSERT_EQ(edata.to, edge.target().global_id());
    }
  }
  printf("+ Pass test: iterate edgelist and get data. :) \n");
  std::cout << "-----------End Grid Test--------------------" << std::endl;
}

#include <graphlab/macros_undef.hpp>
