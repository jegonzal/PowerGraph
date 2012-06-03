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



#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>

#include "als_vertex_program.hpp"

#include <graphlab/macros_def.hpp>


typedef graphlab::synchronous_engine<als_vertex_program> engine_type;

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir, output_dir;
  clopts.attach_option("matrix", &input_dir, input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D",
                       &(vertex_data::NLATENT), vertex_data::NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("maxupdates",
                       &(als_vertex_program::MAX_UPDATES), 
                       als_vertex_program::MAX_UPDATES,
                       "The maxumum number of udpates allowed for a vertex");
  clopts.attach_option("lambda", 
                       &(als_vertex_program::LAMBDA), 
                       als_vertex_program::LAMBDA, 
                       "ALS regularization weight"); 
  clopts.attach_option("tol",
                       &(als_vertex_program::TOLERANCE), 
                       als_vertex_program::TOLERANCE,
                       "residual termination threshold");
  clopts.attach_option("output", &output_dir, output_dir,
                       "Output results");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
  
  std::cout << dc.procid() << ": Loading graph." << std::endl;
  graphlab::timer timer; timer.start();
  graph_type graph(dc, clopts);  
  graphlab::graph_ops::load(graph, input_dir, graph_loader); 
  std::cout << dc.procid() << ": Loading graph. Finished in " 
            << timer.current_time() << std::endl;
  std::cout << dc.procid() << ": Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  std::cout << dc.procid() << ": Finalizing graph. Finished in " 
            << timer.current_time() << std::endl;
  
  if(dc.procid() == 0){
    std::cout
      << "========== Graph statistics on proc " << dc.procid() 
      << " ==============="
      << "\n Num vertices: " << graph.num_vertices()
      << "\n Num edges: " << graph.num_edges()
      << "\n Num replica: " << graph.num_replicas()
      << "\n Replica to vertex ratio: " 
      << float(graph.num_replicas())/graph.num_vertices()
      << "\n --------------------------------------------" 
      << "\n Num local own vertices: " << graph.num_local_own_vertices()
      << "\n Num local vertices: " << graph.num_local_vertices()
      << "\n Replica to own ratio: " 
      << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
      << "\n Num local edges: " << graph.num_local_edges()
      //<< "\n Begin edge id: " << graph.global_eid(0)
      << "\n Edge balance ratio: " 
      << float(graph.num_local_edges())/graph.num_edges()
      << std::endl;
  }
 
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts);
 
  // Run the PageRank ---------------------------------------------------------
  std::cout << "Running ALS" << std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  if(dc.procid() == 0) {
    std::cout << "----------------------------------------------------------"
              << std::endl;
    std::cout << "Final Runtime (seconds):   " << runtime 
              << std::endl
              << "Updates executed: " << engine.last_update_count() << std::endl
              << "Update Rate (updates/second): " 
              << engine.last_update_count() / runtime << std::endl;
  }



  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main



