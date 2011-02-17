#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/distributed_chromatic_engine.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
void generate_atoms() {
  graphlab::graph<size_t, double> testgraph;
  for (size_t v = 0; v < 10000; ++v) testgraph.add_vertex(v);
  for (size_t i = 0;i < 9999; ++i) {
    testgraph.add_edge(i, i+1, i);
  }
  testgraph.add_edge(9999,0,9999);
  std::vector<uint32_t> parts;
  testgraph.partition(graphlab::partition_method::PARTITION_METIS, 4, parts);
  for (size_t i = 0;i < parts.size(); ++i) std::cout << parts[i];
  std::cout << "\n";
  testgraph.compute_coloring();
  graph_partition_to_atomindex(testgraph, parts, "atomidx_ne.txt", "atom_ne", true);
}



using namespace graphlab;

int main(int argc, char** argv) {
  dc_init_param param;
  if (init_param_from_env(param) == false) {
    generate_atoms(); return 0;
  }
  
  global_logger().set_log_level(LOG_DEBUG);

  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
  distributed_control dc(param);

  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
//  distributed_control dc(machines, "", machineid, 1, TCP_COMM);
  dc.barrier();
  distributed_graph<size_t, double> dg(dc, "atomidx_ne.txt");

  dc.barrier();
  std::cout << "Graph Constructed!" << std::endl;
  // now we make an engine
  distributed_chromatic_engine<distributed_graph<size_t, double> > engine(dc, dg, 1);
}
