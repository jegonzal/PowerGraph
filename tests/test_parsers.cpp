#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/distributed_graph_ops.hpp>
#include <graphlab/macros_def.hpp>

typedef graphlab::distributed_graph<size_t, size_t> graph_type;

void check_structure(graph_type &graph) {
  ASSERT_EQ(graph.num_vertices(), 5);
  ASSERT_EQ(graph.num_edges(), 7);
  // check vertex 0 
  {
    graph_type::vertex_type vtype = graph.vertex(0);
    graph_type::local_edge_list_type v0_out = graph_type::local_vertex_type(vtype).out_edges();
    ASSERT_EQ(v0_out.size(), 1);
    ASSERT_EQ(v0_out[0].target().global_id(), 5);
  }
  // vertex 1
  {
    graph_type::vertex_type vtype = graph.vertex(1);
    graph_type::local_edge_list_type v0_out = graph_type::local_vertex_type(vtype).out_edges();
    ASSERT_EQ(v0_out.size(), 2);
    ASSERT_EQ(v0_out[0].target().global_id(), 0);
    ASSERT_EQ(v0_out[1].target().global_id(), 5);
  }
  
  // vertex 2
  {
    graph_type::vertex_type vtype = graph.vertex(2);
    graph_type::local_edge_list_type v0_out = graph_type::local_vertex_type(vtype).out_edges();
    ASSERT_EQ(v0_out.size(), 2);
    ASSERT_EQ(v0_out[0].target().global_id(), 0);
    ASSERT_EQ(v0_out[1].target().global_id(), 5);
  }
  // vertex 3
  {
    graph_type::vertex_type vtype = graph.vertex(3);
    graph_type::local_edge_list_type v0_out = graph_type::local_vertex_type(vtype).out_edges();
    ASSERT_EQ(v0_out.size(), 2);
    ASSERT_EQ(v0_out[0].target().global_id(), 0);
    ASSERT_EQ(v0_out[1].target().global_id(), 5);
  }
}



void test_adj(graphlab::distributed_control& dc) {
  graphlab::distributed_graph<size_t, size_t> graph(dc);
  graphlab::graph_ops::load(graph, "data/test_adj", "adj");
  graph.finalize();
  check_structure(graph);  
}

void test_snap(graphlab::distributed_control& dc) {
  graphlab::distributed_graph<size_t, size_t> graph(dc);
  graphlab::graph_ops::load(graph, "data/test_snap", "snap");
  graph.finalize();
  check_structure(graph);  
}

void test_tsv(graphlab::distributed_control& dc) {
  graphlab::distributed_graph<size_t, size_t> graph(dc);
  graphlab::graph_ops::load(graph, "data/test_tsv", "tsv");
  graph.finalize();
  check_structure(graph);  
}

void test_powerlaw(graphlab::distributed_control& dc) {
  graphlab::distributed_graph<size_t, size_t> graph(dc);
  graphlab::graph_ops::load_synthetic_powerlaw(graph, 1000);
  graph.finalize();
  ASSERT_EQ(graph.num_vertices(), 1000);
  std::cout << graph.num_edges() << " Edges\n";
}

int main(int argc, char** argv) {
  graphlab::distributed_control dc;
  test_adj(dc);
  test_snap(dc);
  test_tsv(dc);
  test_powerlaw(dc);
};

