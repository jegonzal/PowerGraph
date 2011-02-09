#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/graph/graph_lock.hpp>
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

volatile bool firstlockacquired;
void lock_acquired(vertex_id_t vid) {
  std::cout << "lock on " << vid << std::endl;
  firstlockacquired = true;
}


int main(int argc, char** argv) {
  firstlockacquired = false;
  dc_init_param param;
  if (init_param_from_env(param) == false) {
    generate_atoms(); return 0;
  }
  
  global_logger().set_log_level(LOG_DEBUG);

  param.numhandlerthreads = 1;
  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
  distributed_control dc(param);
  dc.services().barrier();
  distributed_graph<size_t, double> dg(dc, "atomidx_ne.txt");
  graph_lock<size_t, double> graphlock(dc, dg);
  dc.services().barrier();
  std::cout << "Constructed!" << std::endl;
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  assert(localvertices.size() > 0);
  graphlock.scope_request(localvertices[0], lock_acquired, scope_range::EDGE_CONSISTENCY);
  graphlock.scope_request(localvertices[1], lock_acquired, scope_range::EDGE_CONSISTENCY);
  while(firstlockacquired == false) sched_yield();
  graphlock.scope_unlock(localvertices[0], scope_range::EDGE_CONSISTENCY);
  sleep(2);
  dc.services().barrier();
}

