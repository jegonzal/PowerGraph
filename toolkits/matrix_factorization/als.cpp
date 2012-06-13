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


/**
 * \file
 * 
 * \brief The main file for the ALS matrix factorization algorithm.
 *
 * This file contains the main body of the ALS matrix factorization
 * algorithm. 
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>

#include "als_vertex_program.hpp"

#include <graphlab/macros_def.hpp>



size_t vertex_data::NLATENT = 20;
double als_vertex_program::TOLERANCE = 1e-3;
double als_vertex_program::LAMBDA = 0.01;
size_t als_vertex_program::MAX_UPDATES = -1;


/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef graphlab::omni_engine<als_vertex_program> engine_type;

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir, output_dir;
  std::string predictions;
  size_t interval = 10;
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
  clopts.attach_option("interval",  &interval, interval, 
                       "The time in seconds between error reports");
  clopts.attach_option("predictions", &predictions, predictions,
                       "The prefix (folder and filename) to save predictions.");
  clopts.attach_option("output", &output_dir, output_dir,
                       "Output results");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; 
  graph_type graph(dc, clopts);  
  graph.load(input_dir, graph_loader); 
  dc.cout() << "Loading graph. Finished in " 
            << timer.current_time() << std::endl;
  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
            << timer.current_time() << std::endl;


  dc.cout() 
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
 
  dc.cout() << "Creating engine" << std::endl;
  engine_type engine(dc, graph, clopts, "synchronous");

  // Add error reporting to the engine
  const bool success = engine.add_edge_aggregator<error_aggregator>
    ("error", error_aggregator::map, error_aggregator::finalize) &&
    engine.aggregate_periodic("error", interval);
  ASSERT_TRUE(success);
  

  // Signal all vertices on the vertices on the left (liberals) 
  engine.map_reduce_vertices<graphlab::empty>(als_vertex_program::signal_left);
 

  // Run the PageRank ---------------------------------------------------------
  dc.cout() << "Running ALS" << std::endl;
  timer.start();
  engine.start();  

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime 
            << std::endl
            << "Updates executed: " << engine.num_updates() << std::endl
            << "Update Rate (updates/second): " 
            << engine.num_updates() / runtime << std::endl;

  // Compute the final training error -----------------------------------------
  dc.cout() << "Final error: " << std::endl;
  engine.aggregate_now("error");

  // Make predictions ---------------------------------------------------------
  if(!predictions.empty()) {
    std::cout << "Saving predictions" << std::endl;
    const bool gzip_output = false;
    const bool save_vertices = false;
    const bool save_edges = true;
    const size_t threads_per_machine = 2;
    graph.save(predictions, prediction_saver(),
               gzip_output, save_vertices, 
               save_edges, threads_per_machine);
    
  }
             


  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



