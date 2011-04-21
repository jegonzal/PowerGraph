#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

#include <graphlab/metrics/reporters/basic_reporter.hpp>
#include <graphlab/metrics/reporters/file_reporter.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/distributed_locking_engine.hpp>
#include <graphlab/distributed2/distributed_glshared.hpp>
#include <graphlab/distributed2/distributed_glshared_manager.hpp>
#include <graphlab/schedulers/multiqueue_fifo_scheduler.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
using namespace graphlab;
  
typedef distributed_graph<size_t, double> graph_type;
typedef multiqueue_fifo_scheduler<graph_type> scheduler_type;
typedef distributed_locking_engine<distributed_graph<size_t, double>, scheduler_type > engine_type;
typedef engine_type::iscope_type iscope_type;
typedef engine_type::icallback_type icallback_type;
typedef ishared_data<graph_type> ishared_data_type;
typedef engine_type::icallback_type icallback_type;
typedef engine_type::update_task_type update_task_type;


const size_t NUMV = 1000;
const size_t NUMITERATIONS = 1000;
const size_t SYNC_INTERVAL = 100;
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
  graph_partition_to_atomindex(testgraph, parts, "atomidx_ne_locktest.txt", "atom_ne_locktest", true);
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
    assert(edata == std::max<size_t>(vdata, destvdata));
  }
  
  vdata = (srcvdata + 1);
  
  
  // write in data
  foreach(edge_id_t edge, scope.in_edge_ids()) {
    scope.edge_data(edge) = std::max<size_t>(vdata, srcvdata);
  }
 
  //write out data
  foreach(edge_id_t edge, scope.out_edge_ids()) {
    scope.edge_data(edge) = std::max<size_t>(vdata, destvdata);
  }
  
  if (vdata < 2 * NUMITERATIONS) {
    scheduler.add_task(update_task_type((scope.vertex() + NUMV - 1) % NUMV, add_one_dynamic),
                      1.0);
    scheduler.add_task(update_task_type((scope.vertex() + 1) % NUMV, add_one_dynamic),
                      1.0);
  }
}







struct accumulator_type {
  double sum;
  size_t count;
  accumulator_type() : sum(0), count(0) { } 
};
SERIALIZABLE_POD(accumulator_type);


void sync_sum_fun(iscope_type& iscope,
                  any& acc) {
  acc.as<accumulator_type>().sum += iscope.vertex_data();
  acc.as<accumulator_type>().count++;
}

void apply_fun(any& current_data,
               const any& acc) {
  current_data.as<double>() =
    acc.as<accumulator_type>().sum /
    acc.as<accumulator_type>().count;
}


void merge_fun(any& merge_dest,
              const any& merge_src) {
  merge_dest.as<accumulator_type>().sum += merge_src.as<accumulator_type>().sum;
  merge_dest.as<accumulator_type>().count += merge_src.as<accumulator_type>().count;
}



distributed_glshared<double> averagevalue;












int main(int argc, char** argv) {

  dc_init_param param;

  // if not running in DC environment, make atoms
  if (init_param_from_env(param) == false) {
    generate_atoms(); return 0;
  }
  
  param.initstring = "";
  param.numhandlerthreads = 8;
  global_logger().set_log_level(LOG_DEBUG);
  distributed_control dc(param);
  
  graph_type dg(dc, "atomidx_ne_locktest.txt");

  std::cout << "Graph Constructed!" << std::endl;
  std::cout << "Testing Static: " << std::endl;
  // now we make an engine
  engine_type engine(dc, dg, 2);
  averagevalue.set(0.0);
  engine.set_sync(averagevalue,
                sync_sum_fun,
                apply_fun,
                any(accumulator_type()),
                SYNC_INTERVAL,
                merge_fun);
  
  /*******************************************************************/
  
  std::cout << "Testing Dynamic: " << std::endl;
  // reset graph
  averagevalue.set(0.0);
  std::cout << "Resetting Graph..." << std::endl;
  for (size_t i = 0;i < dg.num_vertices(); ++i) {
    dg.set_vertex_data(i, 0);
    foreach(vertex_id_t source, dg.in_vertices(i)) {
      dg.set_edge_data(source, i, 0);
    }
    foreach(vertex_id_t target, dg.out_vertices(i)) {
      dg.set_edge_data(i, target, 0);
    }
  }
  dc.full_barrier();
  dc.full_barrier();
  std::cout << "Starting Dynamic..." << std::endl;
  
  engine.add_task_to_all(add_one_dynamic, 1.0);
  engine.start();
  std::cout << "Done!" << std::endl;
  std::cout << "Synced value: " << averagevalue.get_val() << std::endl;

  
  dc.fill_metrics();
  dg.fill_metrics();
  
  if (dc.procid() == 0) {
    basic_reporter reporter;
    dc.report_metrics(reporter);
    dg.report_metrics(reporter);
    engine.report_metrics(reporter);
    file_reporter freporter("graphlab_metrics.txt");
    dc.report_metrics(freporter);
    dg.report_metrics(freporter);
    engine.report_metrics(freporter);
  }
}
