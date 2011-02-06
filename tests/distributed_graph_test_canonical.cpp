#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/logger/assertions.hpp>
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
  if (argc == 1) {
    generate_atoms(); return 0;
  }
  assert(argc == 3);
  size_t machineid = atoi(argv[1]);
  size_t nummachines = atoi(argv[2]);
  global_logger().set_log_level(LOG_DEBUG);
  std::vector<std::string> machines;
  for (size_t i = 0;i < nummachines; ++i) {
    std::stringstream strm;
    strm << "127.0.0.1:" << 10000+i;
    machines.push_back(strm.str());
  }
  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
  distributed_control dc(machines, "", machineid, 8, TCP_COMM);
  dc.services().barrier();
  distributed_graph<size_t, double> dg(dc, "atomidx_ne.txt");

  dc.services().barrier();
  std::cout << "Constructed!" << std::endl;
  
  ASSERT_EQ(dg.num_vertices(), 10000);
  ASSERT_EQ(dg.num_edges(), 10000);
  // everyone do a check on the local graph structure
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  std::set<edge_id_t> eids;
  for (size_t v_ = 0; v_ < localvertices.size(); ++v_) {
    vertex_id_t v = localvertices[v_];
    // everything has one in and one out
    ASSERT_EQ(dg.num_in_neighbors(v), 1);
    ASSERT_EQ(dg.num_out_neighbors(v), 1);
    // check find
    std::pair<bool, edge_id_t> ret = dg.find(v, (v+1) % 10000);
    ASSERT_TRUE(ret.first);

    // check source and target
    ASSERT_EQ(dg.source(ret.second), v);
    ASSERT_EQ(dg.target(ret.second), (v+1) % 10000);

    // check in edgeids
    dgraph_edge_list delist = dg.in_edge_ids(v);
    ASSERT_EQ(delist.size(), 1);
    eids.insert(*delist.begin());
    delist = dg.out_edge_ids(v);
    ASSERT_EQ(delist.size(), 1);
    eids.insert(*delist.begin());
  }

  // check for data
  for (size_t v_ = 0; v_ < localvertices.size(); ++v_) {
    vertex_id_t v = localvertices[v_];
    ASSERT_EQ(dg.get_vertex_data(v), v);
    ASSERT_EQ(dg.get_edge_data(v, (v+1) % 10000), v);
    if (dg.vertex_is_local(v)) {
      ASSERT_EQ(dg.vertex_data(v), v);
    }
  }
  
  std::vector<edge_id_t> alledges;
  std::copy(eids.begin(), eids.end(), std::back_inserter(alledges));
  for (size_t e_ = 0; e_ < alledges.size(); ++e_) {
    edge_id_t e = alledges[e_];
    
    ASSERT_EQ(dg.get_edge_data(e), dg.source(e));
    if (dg.edge_is_local(e)) {
      ASSERT_EQ(dg.edge_data(e), dg.source(e));
    }
  }
  dc.services().barrier();
}
