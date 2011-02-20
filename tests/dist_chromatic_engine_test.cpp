#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include <graphlab/metrics/reporters/basic_reporter.hpp>
#include <graphlab/metrics/reporters/file_reporter.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/distributed_chromatic_engine.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
using namespace graphlab;
  
typedef distributed_graph<size_t, double> graph_type;
typedef distributed_chromatic_engine<distributed_graph<size_t, double> > engine_type;
typedef engine_type::iscope_type iscope_type;
typedef engine_type::icallback_type icallback_type;
typedef ishared_data<graph_type> ishared_data_type;
typedef engine_type::icallback_type icallback_type;
typedef engine_type::update_task_type update_task_type;


const size_t NUMV = 10000;
const size_t NUMITERATIONS = 1000;

void generate_atoms() {
  graph<size_t, double> testgraph;
  for (size_t v = 0; v < NUMV; ++v) testgraph.add_vertex(0);
  for (size_t i = 0;i < NUMV - 1; ++i) {
    testgraph.add_edge(i, i+1, 0);
  }
  testgraph.add_edge(NUMV - 1,0, 0);
  std::vector<uint32_t> parts;
  testgraph.partition(graphlab::partition_method::PARTITION_METIS, 4, parts);
  for (size_t i = 0;i < parts.size(); ++i) std::cout << parts[i];
  std::cout << "\n";
  testgraph.compute_coloring();
  graph_partition_to_atomindex(testgraph, parts, "atomidx_ne_chrtest.txt", "atom_ne_chrtest", true);
}

/*
 * The update function here is as such:
 *  - Each vertex will look at all of its neighbors (which should all be
 *                                          the same value in this graph)
 *    and set itself to be one higher than its neighbors value
 *  - Each edge is updated to be the max of its source or sink
 */

void add_one_static(iscope_type& scope,
                    icallback_type& scheduler,
                    ishared_data_type* data_manager) {
  size_t& vdata = scope.vertex_data();
  //logger(LOG_INFO, "eval on %d", scope.vertex());
  
  size_t srcvdata , destvdata;
  
  // try to read in data
  foreach(edge_id_t edge, scope.in_edge_ids()) {
    srcvdata = scope.const_neighbor_vertex_data(scope.source(edge));
    assert(scope.source(edge) == (scope.vertex() + NUMV - 1) % NUMV);
    // sanity check
    double edata = scope.const_edge_data(edge);
    assert(edata == std::max<size_t>(vdata, srcvdata));
  }
 
  //try to read out data
  foreach(edge_id_t edge, scope.out_edge_ids()) {
    destvdata = scope.const_neighbor_vertex_data(scope.target(edge));
    assert(scope.target(edge) == (scope.vertex() + 1) % NUMV);
    // sanity check
    double edata = scope.const_edge_data(edge);
    assert(edata == std::max<size_t>(vdata, srcvdata));
  }
  
  // consistency check
  ASSERT_EQ(srcvdata, destvdata);
  
  vdata = (srcvdata + 1);
  
  
  // write in data
  foreach(edge_id_t edge, scope.in_edge_ids()) {
    scope.edge_data(edge) = std::max<size_t>(vdata, srcvdata);
  }
 
  //write out data
  foreach(edge_id_t edge, scope.out_edge_ids()) {
    scope.edge_data(edge) = std::max<size_t>(vdata, srcvdata);
  }
}


void add_one_dynamic(iscope_type& scope,
                    icallback_type& scheduler,
                    ishared_data_type* data_manager) {
  size_t& vdata = scope.vertex_data();
  //logger(LOG_INFO, "eval on %d", scope.vertex());
  
  size_t srcvdata , destvdata;
  
  // try to read in data
  foreach(edge_id_t edge, scope.in_edge_ids()) {
    srcvdata = scope.const_neighbor_vertex_data(scope.source(edge));
    assert(scope.source(edge) == (scope.vertex() + NUMV - 1) % NUMV);
    // sanity check
    double edata = scope.const_edge_data(edge);
    assert(edata == std::max<size_t>(vdata, srcvdata));
  }
 
  //try to read out data
  foreach(edge_id_t edge, scope.out_edge_ids()) {
    destvdata = scope.const_neighbor_vertex_data(scope.target(edge));
    assert(scope.target(edge) == (scope.vertex() + 1) % NUMV);
    // sanity check
    double edata = scope.const_edge_data(edge);
    assert(edata == std::max<size_t>(vdata, srcvdata));
  }
  
  // consistency check
  ASSERT_EQ(srcvdata, destvdata);
  
  vdata = (srcvdata + 1);
  
  
  // write in data
  foreach(edge_id_t edge, scope.in_edge_ids()) {
    scope.edge_data(edge) = std::max<size_t>(vdata, srcvdata);
  }
 
  //write out data
  foreach(edge_id_t edge, scope.out_edge_ids()) {
    scope.edge_data(edge) = std::max<size_t>(vdata, srcvdata);
  }
  
  if (vdata < 2 * NUMITERATIONS) {
    scheduler.add_task(update_task_type((scope.vertex() + NUMV - 1) % NUMV, add_one_dynamic),
                      1.0);
    scheduler.add_task(update_task_type((scope.vertex() + 1) % NUMV, add_one_dynamic),
                      1.0);
  }
}


int main(int argc, char** argv) {
  dc_init_param param;

  // if not running in DC environment, make atoms
  if (init_param_from_env(param) == false) {
    generate_atoms(); return 0;
  }
  
  global_logger().set_log_level(LOG_DEBUG);
  distributed_control dc(param);

  graph_type dg(dc, "atomidx_ne_chrtest.txt");

  std::cout << "Graph Constructed!" << std::endl;
  std::cout << "Testing Static: " << std::endl;
  // now we make an engine
  distributed_chromatic_engine<distributed_graph<size_t, double> > engine(dc, dg, 2);
  scheduler_options schedopts;
  schedopts.add_option("update_function", add_one_static);
  schedopts.add_option("max_iterations", NUMITERATIONS);
  engine.set_scheduler_options(schedopts);
  engine.start();
  //return 0;
  //now go through graph again and do a consistency check.
  // should be alternating between 2 * NUMITERATIONS - 1 and 2 * NUMITERATIONS
  bool ismatch = dg.get_vertex_data(0) == 2 * NUMITERATIONS;
  for (size_t i = 0;i < dg.num_vertices(); ++i) {
    if (ismatch) {
      ASSERT_EQ(dg.get_vertex_data(i), 2 * NUMITERATIONS);
    }
    else {
      ASSERT_EQ(dg.get_vertex_data(i), 2 * NUMITERATIONS - 1);      
    }
    // reset
    ismatch = !ismatch;
    foreach(vertex_id_t source, dg.in_vertices(i)) {
      size_t mval = std::max<size_t>(dg.get_vertex_data(source),
                                     dg.get_vertex_data(i));
      
      ASSERT_EQ(dg.get_edge_data(source, i), double(mval));
    }
    foreach(vertex_id_t target, dg.out_vertices(i)) {
      size_t mval = std::max<size_t>(dg.get_vertex_data(target),
                                     dg.get_vertex_data(i));
      
      ASSERT_EQ(dg.get_edge_data(i, target), double(mval));
    }
  }
  
  
  dc.barrier();
  
  /*******************************************************************/
  
  std::cout << "Testing Dynamic: " << std::endl;
  // reset graph
  for (size_t i = 0;i < dg.num_vertices(); ++i) {
    dg.set_vertex_data(i, 0);
    foreach(vertex_id_t source, dg.in_vertices(i)) {
      dg.set_edge_data(source, i, 0);
    }
    foreach(vertex_id_t target, dg.out_vertices(i)) {
      dg.set_edge_data(i, target, 0);
    }
  }
 
  schedopts.add_option("update_function", add_one_dynamic);
  schedopts.add_option("max_iterations", 0);
  engine.set_scheduler_options(schedopts);
  
  engine.add_task_to_all(add_one_dynamic, 1.0);
  engine.start();
  std::cout << "Done!" << std::endl;
  ismatch = dg.get_vertex_data(0) == 2 * NUMITERATIONS;
  for (size_t i = 0;i < dg.num_vertices(); ++i) {
    if (ismatch) {
      ASSERT_EQ(dg.get_vertex_data(i), 2 * NUMITERATIONS);
    }
    else {
      ASSERT_EQ(dg.get_vertex_data(i), 2 * NUMITERATIONS - 1);      
    }
    // reset
    ismatch = !ismatch;
  }
  
  dc.fill_metrics();
  dg.fill_metrics();
  
  if (dc.procid() == 0) {
    basic_reporter reporter;
    metrics::report_all(reporter);
    file_reporter freporter("graphlab_metrics.txt");
    metrics::report_all(freporter);
  }
}
