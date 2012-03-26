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


#include <string>

#include <distributed_graphlab.hpp>

#include "data_structures.hpp"


#include <graphlab/macros_def.hpp>


// Constants for the algorithm. Better way would be to
// pass them in the shared_data to the update function, but for
// the simplicity of this example, we simply define them here.

double termination_bound = 1e-5;
double random_reset_prob = 0.15;   // PageRank random reset probability

/**
 * The Page rank update function
 */
void pagerank_update(gl_types::iscope &scope,
                     gl_types::icallback &scheduler) {                     
  // Get the data associated with the vertex
  vertex_data& vdata = scope.vertex_data();
  
  // Sum the incoming weights; start by adding the 
  // contribution from a self-link.
  float sum = 0;
  foreach(graphlab::edge_id_t eid, scope.in_edge_ids()) {
    // Get the neighobr vertex value
    const vertex_data& neighbor_vdata =
      scope.const_neighbor_vertex_data(scope.source(eid));    
    // Add the contribution to the sum
    sum += neighbor_vdata.value /  float(scope.out_edge_ids(scope.source(eid)).size());
  }
  // compute the jumpweight
  const float old_value = vdata.value;
  vdata.value = random_reset_prob + (1-random_reset_prob)*sum;
   
  // Schedule the neighbors as needed
  const float weight = 1.0/float(scope.out_edge_ids().size());
  const float residual =  std::fabs(old_value - vdata.value)*weight;
  if(residual > termination_bound) {  
    foreach(graphlab::edge_id_t eid, scope.out_edge_ids()) {
      gl_types::update_task task(scope.target(eid), pagerank_update);
      scheduler.add_task(task, residual);
    }
  }
} // end of pagerank update function






int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");

  // Setup the parser
  graphlab::command_line_options clopts("Run the PageRank algorithm.");
  clopts.use_distributed_options();


  // Add some command line options
  std::string graph_file;
  clopts.attach_option("graph",
                       &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.attach_option("bound",
                       &termination_bound, termination_bound,
                       "residual termination threshold");
  clopts.attach_option("resetprob",
                       &random_reset_prob, random_reset_prob,
                       "Random reset probability");

  // Parse the command line input
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }


  graphlab::mpi_tools::init(argc, argv);  
  graphlab::dc_init_param param;
  param.initstring += ", buffered_multiqueue_send=1";
  if(graphlab::init_param_from_env(param) == false) {
    graphlab::init_param_from_mpi(param); 
  }
  // create distributed control
  graphlab::distributed_control dc(param);
  // Create the distributed_graph --
  gl_types::distributed_core core(dc, graph_file, 
                                  graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM);
  // Set the engine options
  core.set_engine_options(clopts);
  core.build_engine();

  std::cout << "Scheduling update functions." << std::endl;
  core.sched_options().add_option("update_function", pagerank_update);
  core.add_task_to_all(pagerank_update, 100.0);

  // Starte the engine
  double runtime = core.start();
  if (dc.procid() == 0) {
    size_t update_count = core.last_update_count();
    std::cout << "Finished Running engine in " << runtime 
              << " seconds." << std::endl
              << "Total updates: " << update_count << std::endl
              << "Efficiency: " << (double(update_count) / runtime)
              << " updates per second "
              << std::endl;  
  }

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;

} // End of main

