#include <iostream>
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
  graph_partition_to_atomindex(testgraph, parts, "atomidx.txt", "atom");
}

using namespace graphlab;
int main(int argc, char** argv) {
  //generate_atoms(); return 0;
  size_t machineid = atoi(argv[1]);
  global_logger().set_log_level(LOG_INFO);
  std::vector<std::string> machines;
  machines.push_back("127.0.0.1:10000");
  machines.push_back("127.0.0.1:10001");
  machines.push_back("127.0.0.1:10002");
  machines.push_back("127.0.0.1:10003");
  //distributed_control dc(machines,"buffered_send=yes,buffered_recv=yes", machineid, 8, SCTP_COMM);
  distributed_control dc(machines, "", machineid, 8, TCP_COMM);
  dc.services().barrier();
  distributed_graph<size_t, double> dg(dc, "atomidx.txt");

  dc.services().barrier();
  std::cout << "Constructed!" << std::endl;
  
  ASSERT_EQ(dg.num_vertices(), 10000);
  ASSERT_EQ(dg.num_edges(), 10000);
  // everyone do a check on the global graph structure
  for (size_t v = 0; v < 10000; ++v) {
    // everything has one in and one out
    ASSERT_EQ(dg.num_in_neighbors(v), 1);
    ASSERT_EQ(dg.num_out_neighbors(v), 1);
    // check find
    std::pair<bool, edge_id_t> ret = dg.find(v, (v+1) % 10000);
    ASSERT_TRUE(ret.first);
    ASSERT_EQ(ret.second, v);

    // check source and target
    ASSERT_EQ(dg.source(ret.second), v);
    ASSERT_EQ(dg.target(ret.second), (v+1) % 10000);

    // check in edgeids
    dgraph_edge_list delist = dg.in_edge_ids(v);
    ASSERT_EQ(delist.size(), 1);
    ASSERT_EQ(*delist.begin(), (v + 10000 - 1) % 10000);

    delist = dg.out_edge_ids(v);
    ASSERT_EQ(delist.size(), 1);
    ASSERT_EQ(*delist.begin(), v);
  }

  // check for a few nonexistant edges
  ASSERT_FALSE(dg.find(0,5000).first);
  ASSERT_FALSE(dg.find(5000,6000).first);
  ASSERT_FALSE(dg.find(9999,1).first);

  // check for data
  for (size_t v = 0; v < 10000; ++v) {
    ASSERT_EQ(dg.get_vertex_data(v), v);
    ASSERT_EQ(dg.get_edge_data(v, (v+1) % 10000), v);
    if (dg.vertex_is_local(v)) {
      ASSERT_EQ(dg.vertex_data(v), v);
    }
  }
  
  for (size_t e = 0; e < 10000; ++e) {
    ASSERT_EQ(dg.get_edge_data(e), e);
    if (dg.edge_is_local(e)) {
      ASSERT_EQ(dg.edge_data(e), e);
    }
  }
  dc.services().barrier();
}
