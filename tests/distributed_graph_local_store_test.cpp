
#include <graphlab/distributed2/graph/graph_local_store.hpp>

int main(int argc, char** argv) {
  graphlab::dist_graph_impl::graph_local_store<size_t, double> graph;
  graph.create_store(10000,10000,"v.txt","e.txt");
  for (size_t i = 0;i < 9999; ++i) {
    graph.add_edge(i, i, i+1);
  }
  graph.add_edge(9999,0,9999);
  graph.compute_minimal_prefetch();
  
  for (size_t i = 0;i < 10; ++i) {
    graph.print_prefetch_list(i);
  }
  for (size_t i = 9990;i < 10000; ++i) {
    graph.print_prefetch_list(i);
  }

  graph.edge_data(5) = 10.1;
  graph.vertex_data(10) = 11;
  graph.flush();
}
