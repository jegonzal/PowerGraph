#include <iostream>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>

void generate_atoms() {
  graphlab::graph<size_t, double> testgraph;
  for (size_t v = 0; v < 10000; ++v) testgraph.add_vertex(65);
  for (size_t i = 0;i < 9999; ++i) {
    testgraph.add_edge(i, i+1, 66);
  }
  testgraph.add_edge(9999,0,66);
  std::vector<uint32_t> parts;
  testgraph.partition(graphlab::partition_method::PARTITION_METIS, 4, parts);
  for (size_t i = 0;i < parts.size(); ++i) std::cout << parts[i];
  std::cout << "\n";
  testgraph.compute_coloring();
  graph_partition_to_atomindex(testgraph, parts, "atomidx.txt", "atom");
}

using namespace graphlab;
int main(int argc, char** argv) {
//  generate_atoms();
  size_t machineid = atoi(argv[1]);
  global_logger().set_log_level(LOG_INFO);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");
  distributed_control dc(machines,"", machineid, 8, TCP_COMM);
  dc.services().barrier();
  distributed_graph<size_t, double> dg(dc, "atomidx.txt");
  if (dc.procid() == 0) {
    std::cout << "----------- P0 -----------\n";
    std::cout << dg;
  }
  dc.services().barrier();
  if (dc.procid() == 1) {
    std::cout << "----------- P1 -----------\n";
    std::cout << dg;
  }
  dc.services().barrier();
}
