/**  
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


#include <pthread.h>


#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/graph_partitioner.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/graph/disk_graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>


using namespace graphlab;


void generate_atoms() {
  graphlab::graph<size_t, double> testgraph;
  for (size_t v = 0; v < 10000; ++v) testgraph.add_vertex(v);
  for (size_t i = 0;i < 9999; ++i) {
    testgraph.add_edge(i, i+1, i);
  }
  testgraph.add_edge(9999,0,9999);
  std::vector<graph_partitioner::part_id_type> parts;
  graph_partitioner::partition("metis", testgraph, 4, parts);
  testgraph.compute_coloring();
   
  graphlab::disk_graph<size_t, double> dg("atom_ne", 4);
  dg.create_from_graph(testgraph, parts);
  dg.finalize();
}



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
    dg.vertex_is_modified(localvertices[i]);
  }
}

void set_all_edges_to_value(distributed_graph<size_t, double> &dg, double value) {
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  for (size_t i = 0;i < localvertices.size(); ++i) {
    foreach(edge_id_t eid, dg.in_edge_ids(localvertices[i])) {
      dg.edge_data(eid) = value;
      dg.edge_is_modified(eid);
    }
  }
}

void set_all_in_boundary(distributed_graph<size_t, double> &dg, size_t vvalue, double evalue) {
  const std::vector<vertex_id_t>& boundvertices = dg.boundary_scopes();
  for (size_t i = 0;i < boundvertices.size(); ++i) {
    dg.vertex_data(boundvertices[i]) = vvalue;
    dg.vertex_is_modified(boundvertices[i]);

    foreach(edge_id_t eid, dg.in_edge_ids(boundvertices[i])) {
      dg.edge_data(eid) = evalue;
      dg.edge_is_modified(eid);

      vertex_id_t sourcevid = dg.source(eid);
      dg.vertex_data(sourcevid) = vvalue;
      dg.vertex_is_modified(sourcevid);
    }

    foreach(edge_id_t eid, dg.out_edge_ids(boundvertices[i])) {
      dg.edge_data(eid) = evalue;
      dg.edge_is_modified(eid);

      vertex_id_t targetvid = dg.target(eid);
      dg.vertex_data(targetvid) = vvalue;
      dg.vertex_is_modified(targetvid);
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

void print_usage() {
  std::cout << "Tests distributed graph\n";
  std::cout << "First run ./distributed_graph_test -g to generate the test graph\n";
  std::cout << "Then run distributed_graph_test -b with exactly 2 MPI nodes to test\n";
  std::cout << "     (ex: mpiexec -n 2 ./distributed_graph_test -b)\n";
}

int main(int argc, char** argv) {
  if (argc != 2) {
    print_usage();
    return 0;
  }
  
  if (std::string(argv[1]) == "-g") {
    generate_atoms(); 
    return 0;
  }
  else if (std::string(argv[1]) != "-b") {
    print_usage();
    return 0;
  }
  graphlab::mpi_tools::init(argc, argv);
  
  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  ASSERT_EQ(param.machines.size(), 2);
  global_logger().set_log_level(LOG_DEBUG);

  distributed_control dc(param);

  dc.barrier();
  typedef distributed_graph<size_t, double> graph_type;
  graph_type dg(dc, "atom_ne.idx");

  dc.barrier();
  std::cout << "Constructed!" << std::endl;

  ASSERT_EQ(dg.num_vertices(), 10000);
  ASSERT_EQ(dg.num_edges(), 10000);
  // everyone do a check on the local graph structure
  const std::vector<vertex_id_t>& localvertices = dg.owned_vertices();
  std::set<edge_id_t> eids;
  std::cout << "Checking Graph Structure..." << std::endl;
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
    graph_type::edge_list_type delist = dg.in_edge_ids(v);
    ASSERT_EQ(delist.size(), 1);
    eids.insert(*delist.begin());
    delist = dg.out_edge_ids(v);
    ASSERT_EQ(delist.size(), 1);
    eids.insert(*delist.begin());
  }

  std::cout << "Checking Graph Data..." << std::endl;
  // check for data
  for (size_t v_ = 0; v_ < localvertices.size(); ++v_) {
    vertex_id_t v = localvertices[v_];
    ASSERT_EQ(dg.get_vertex_data(v), v);
    ASSERT_EQ(dg.get_edge_data(v, (v+1) % 10000), v);
    if (dg.vertex_is_local(v)) {
      ASSERT_EQ(dg.vertex_data(v), v);
    }
  }
  
  std::cout << "Testing one way vertex collection..." << std::endl;
  // check vertex collection routines
  // each machine collects a different random subset
  std::vector<vertex_id_t> vsubset_collect;
  for (size_t i = 0;i < 100; ++i) {
    vsubset_collect.push_back(rand() % 10000);
  }
  std::map<vertex_id_t, size_t> retv = dg.collect_vertex_subset_one_way(vsubset_collect);
  // check the map
  for (size_t i = 0;i < 100; ++i) {
    ASSERT_TRUE(retv.find(vsubset_collect[i]) != retv.end());
    ASSERT_EQ(retv[vsubset_collect[i]], vsubset_collect[i]);
  }
  
  dc.barrier();
  
  std::cout << "Testing all vertex collection..." << std::endl;
  // try collect all vertices
  std::vector<size_t> allvertexdata = dg.collect_vertices(0);
  if (dc.procid() == 0) {
    ASSERT_EQ(allvertexdata.size(), 10000);
    for (size_t i = 0;i < 10000; ++i) {
      ASSERT_EQ(allvertexdata[i], i);
    }
  }
  else {
    ASSERT_EQ(allvertexdata.size(), 0);
  }
  

  
  std::vector<edge_id_t> alledges;
  std::copy(eids.begin(), eids.end(), std::back_inserter(alledges));
  for (size_t e_ = 0; e_ < alledges.size(); ++e_) {
    edge_id_t e = alledges[e_];
    
    ASSERT_EQ(dg.get_edge_data(e), dg.source(e));
    ASSERT_EQ(dg.edge_data(e), dg.source(e));
  }
  dc.full_barrier();
  sync_test(dg, dc);
  graphlab::mpi_tools::finalize();
}
