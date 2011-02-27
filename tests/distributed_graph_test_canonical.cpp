#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
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


void check_vertex_values(distributed_graph<size_t, double> &dg, size_t value) {
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  const std::vector<vertex_id_t>& ghostvertices = dg.ghost_vertices();
  for (size_t i = 0;i < localvertices.size(); ++i) {
    ASSERT_EQ(dg.vertex_data(localvertices[i]), value);
  }
  for (size_t i = 0;i < ghostvertices.size(); ++i) {
    ASSERT_EQ(dg.vertex_data(ghostvertices[i]), value);
  }
}

void check_edge_values(distributed_graph<size_t, double> &dg, double value) {
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  for (size_t i = 0;i < localvertices.size(); ++i) {
    foreach(edge_id_t eid, dg.in_edge_ids(localvertices[i])) {
      ASSERT_EQ(dg.edge_data(eid), value);
    }
    foreach(edge_id_t eid, dg.out_edge_ids(localvertices[i])) {
      ASSERT_EQ(dg.edge_data(eid), value);
    }
  }
}

void set_all_vertices_to_value(distributed_graph<size_t, double> &dg, size_t value) {
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  // ok. now everyone set to zero
  for (size_t i = 0;i < localvertices.size(); ++i) {
    dg.vertex_data(localvertices[i]) = value;
    dg.increment_vertex_version(localvertices[i]);
  }
}

void set_all_edges_to_value(distributed_graph<size_t, double> &dg, double value) {
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  for (size_t i = 0;i < localvertices.size(); ++i) {
    foreach(edge_id_t eid, dg.in_edge_ids(localvertices[i])) {
      dg.edge_data(eid) = value;
      dg.increment_edge_version(eid);
    }
  }
}

void set_all_in_boundary(distributed_graph<size_t, double> &dg, size_t vvalue, double evalue) {
  const std::vector<vertex_id_t>& boundvertices = dg.boundary_scopes();
  for (size_t i = 0;i < boundvertices.size(); ++i) {
    dg.vertex_data(boundvertices[i]) = vvalue;
    dg.vertex_is_modified(boundvertices[i]);
    dg.increment_vertex_version(boundvertices[i]);

    foreach(edge_id_t eid, dg.in_edge_ids(boundvertices[i])) {
      dg.edge_data(eid) = evalue;
      dg.edge_is_modified(eid);
      dg.increment_edge_version(eid);

      vertex_id_t sourcevid = dg.source(eid);
      dg.vertex_data(sourcevid) = vvalue;
      dg.vertex_is_modified(sourcevid);
      dg.increment_vertex_version(sourcevid);
    }

    foreach(edge_id_t eid, dg.out_edge_ids(boundvertices[i])) {
      dg.edge_data(eid) = evalue;
      dg.edge_is_modified(eid);
      dg.increment_edge_version(eid);

      vertex_id_t targetvid = dg.target(eid);
      dg.vertex_data(targetvid) = vvalue;
      dg.vertex_is_modified(targetvid);
      dg.increment_vertex_version(targetvid);
    }
  }
}


void boundary_has_at_least_one_match(distributed_graph<size_t, double> &dg, size_t vvalue, double evalue) {
  // test only owned data
  bool vmatch = false;
  bool ematch = false;
  const std::vector<vertex_id_t>& boundvertices = dg.boundary_scopes();
  for (size_t i = 0;i < boundvertices.size(); ++i) {
    if (dg.vertex_data(boundvertices[i]) == vvalue) vmatch = true;

    foreach(edge_id_t eid, dg.in_edge_ids(boundvertices[i])) {
      if (dg.edge_data(eid) == evalue) ematch = true;
    }
  }
  ASSERT_TRUE(vmatch);
  ASSERT_TRUE(ematch);
}

void sync_test(distributed_graph<size_t, double> &dg, distributed_control &dc) {
  size_t VVAL = 0;
  double EVAL = 0;
  std::cout << "Owned data modified. Test if Ghost data are updated" << std::endl;
  std::cout << "===================================================" << std::endl;
  std::cout << "Testing Synchronous Vertex sync. " << std::endl;
  set_all_vertices_to_value(dg, VVAL);
  dc.barrier();
  dg.synchronize_all_vertices();
  check_vertex_values(dg, VVAL);
  dc.barrier();

  ++VVAL;

  std::cout << "Testing Asynchronous Vertex sync. " << std::endl;
  set_all_vertices_to_value(dg, VVAL);
  dc.barrier();
  dg.synchronize_all_vertices(true);
  dg.wait_for_all_async_syncs();
  check_vertex_values(dg, VVAL);
  dc.barrier();

  ++VVAL;

  std::cout << "Testing Synchronous Edge sync. " << std::endl;
  set_all_edges_to_value(dg, EVAL);
  dc.barrier();
  dg.synchronize_all_edges();
  check_edge_values(dg, EVAL);
  dc.barrier();

  ++EVAL;

  std::cout << "Testing Asynchronous Edge sync. " << std::endl;
  set_all_edges_to_value(dg, EVAL);
  dc.barrier();
  dg.synchronize_all_edges(true);
  dg.wait_for_all_async_syncs();
  check_edge_values(dg, EVAL);
  dc.barrier();

  ++EVAL;
  ++VVAL;

  std::cout << "Testing Synchronous Scope sync. " << std::endl;
  set_all_vertices_to_value(dg, VVAL);
  set_all_edges_to_value(dg, EVAL);
  dc.barrier();
  dg.synchronize_all_scopes();
  check_vertex_values(dg, VVAL);
  check_edge_values(dg, EVAL);
  dc.barrier();

  ++EVAL;
  ++VVAL;

  std::cout << "Testing Asynchronous Scope sync. " << std::endl;
  set_all_vertices_to_value(dg, VVAL);
  set_all_edges_to_value(dg, EVAL);
  dc.barrier();
  dg.synchronize_all_scopes(true);
  dg.wait_for_all_async_syncs();
  check_vertex_values(dg, VVAL);
  check_edge_values(dg, EVAL);
  dc.barrier();



  ++VVAL;
  
  std::cout << "Testing synchronous vertex pushing" << std::endl;
  set_all_vertices_to_value(dg, VVAL);
  dg.push_all_owned_vertices_to_replicas();
  dc.barrier();
  check_vertex_values(dg, VVAL);
  dc.barrier();
  
  ++VVAL;  
  
  std::cout << "Testing asynchronous vertex pushing" << std::endl;
  set_all_vertices_to_value(dg, VVAL);
  dg.push_all_owned_vertices_to_replicas();
  dg.wait_for_all_async_pushes();
  dc.barrier();
  check_vertex_values(dg, VVAL);
  dc.barrier();    
 
 
 
  ++EVAL;
  
  std::cout << "Testing synchronous edge pushing" << std::endl;
  set_all_edges_to_value(dg, EVAL);
  dg.push_all_owned_edges_to_replicas();
  dc.barrier();
  check_edge_values(dg, EVAL);
  dc.barrier();
  
  ++EVAL;  
  
  std::cout << "Testing asynchronous edge pushing" << std::endl;
  set_all_edges_to_value(dg, EVAL);
  dg.push_all_owned_edges_to_replicas();
  dg.wait_for_all_async_pushes();
  dc.barrier();
  check_edge_values(dg, EVAL);
  dc.barrier();    

  std::cout << "Ghost data modified. Test if remote owned data are updated" << std::endl;
  std::cout << "==========================================================" << std::endl;

  std::cout << "Testing Synchronous Scope sync with ghost changes. " << std::endl;
  ++EVAL;
  ++VVAL;
  if (dc.procid() == 0) {
    set_all_in_boundary(dg, VVAL, EVAL);
    dg.synchronize_all_scopes();
  }
  dc.barrier();
  if (dc.procid() == 1) {
    boundary_has_at_least_one_match(dg, VVAL, EVAL);
    dg.synchronize_all_scopes();
  }
  dc.barrier();


  std::cout << "Testing Synchronous Vertex/Edge sync with ghost changes. " << std::endl;
  ++EVAL;
  ++VVAL;
  if (dc.procid() == 1) {
    set_all_in_boundary(dg, VVAL, EVAL);
    dg.synchronize_all_vertices();
    dg.synchronize_all_edges();
  }
  dc.barrier();
  if (dc.procid() == 0) {
    boundary_has_at_least_one_match(dg, VVAL, EVAL);
    dg.synchronize_all_scopes();
  }
  dc.barrier();



  std::cout << "Testing Asynchronous Scope sync with ghost changes. " << std::endl;
  ++EVAL;
  ++VVAL;
  if (dc.procid() == 0) {
    set_all_in_boundary(dg, VVAL, EVAL);
    dg.synchronize_all_scopes(true);
    dg.wait_for_all_async_syncs();
  }
  dc.barrier();
  if (dc.procid() == 1) {
    boundary_has_at_least_one_match(dg, VVAL, EVAL);
    dg.synchronize_all_scopes();
  }
  dc.barrier();


  std::cout << "Testing Asynchronous Vertex/Edge sync with ghost changes. " << std::endl;
  ++EVAL;
  ++VVAL;
  if (dc.procid() == 1) {
    set_all_in_boundary(dg, VVAL, EVAL);
    dg.synchronize_all_vertices(true);
    dg.synchronize_all_edges(true);
    dg.wait_for_all_async_syncs();
  }
  dc.barrier();
  if (dc.procid() == 0) {
    boundary_has_at_least_one_match(dg, VVAL, EVAL);
    dg.synchronize_all_scopes();
  }
  dc.barrier();
}

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
  dc.full_barrier();
  sync_test(dg, dc);
}
