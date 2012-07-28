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
 *      http://graphlab.org
 *
 */


/**
 * \file
 * 
 * \brief The main file for the BIAS-SGD matrix factorization algorithm.
 *
 * This file contains the main body of the BIAS-SGD matrix factorization
 * algorithm. 
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab.hpp>
#include "eigen_serialization.hpp"
#include "biassgd_vertex_program.hpp"

#include <graphlab/macros_def.hpp>



size_t vertex_data::NLATENT = 20;
double biassgd_vertex_program::TOLERANCE = 1e-3;
double biassgd_vertex_program::LAMBDA = 0.001;
double biassgd_vertex_program::GAMMA = 0.001;
size_t biassgd_vertex_program::MAX_UPDATES = -1;
double biassgd_vertex_program::MAXVAL = 1e+100;
double biassgd_vertex_program::MINVAL = 1e-100;
double biassgd_vertex_program::STEP_DEC = 0.9;
bool biassgd_vertex_program::debug = false;
uint biassgd_vertex_program::USERS = 0;
double biassgd_vertex_program::GLOBAL_MEAN = 0;
size_t biassgd_vertex_program::NUM_TRAINING_EDGES = 0;

/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef graphlab::omni_engine<biassgd_vertex_program> engine_type;

double calc_global_mean(const graph_type::edge_type & edge){
  if (edge.data().role == edge_data::TRAIN)
     return edge.data().obs;
  else return 0;
}

size_t count_edges(const graph_type::edge_type & edge){
  if (edge.data().role == edge_data::TRAIN)
     return 1;
  else return 0;
}


int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir, output_dir;
  std::string predictions;
  size_t interval = 0;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D", vertex_data::NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("engine", exec_type, 
                       "The engine type synchronous or asynchronous");
  clopts.attach_option("max_iter", biassgd_vertex_program::MAX_UPDATES,
                       "The maxumum number of udpates allowed for a vertex");
  clopts.attach_option("lambda", biassgd_vertex_program::LAMBDA, 
                       "SGD regularization weight"); 
  clopts.attach_option("gamma", biassgd_vertex_program::GAMMA, 
                       "SGD step size"); 
  clopts.attach_option("debug", biassgd_vertex_program::debug, 
                       "debug - additional verbose info"); 
  clopts.attach_option("tol", biassgd_vertex_program::TOLERANCE,
                       "residual termination threshold");
  clopts.attach_option("maxval", biassgd_vertex_program::MAXVAL, "max allowed value");
  clopts.attach_option("minval", biassgd_vertex_program::MINVAL, "min allowed value");
  clopts.attach_option("step_dec", biassgd_vertex_program::STEP_DEC, "multiplicative step decrement");
  clopts.attach_option("interval", interval, 
                       "The time in seconds between error reports");
  clopts.attach_option("predictions", predictions,
                       "The prefix (folder and filename) to save predictions.");
  clopts.attach_option("output", output_dir,
                       "Output results");
  clopts.attach_option("users", biassgd_vertex_program::USERS, "number of users");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
    if (biassgd_vertex_program::USERS == 0){
    logstream(LOG_FATAL)<<"Please specify the number of users using the --users=XX command line argument"<<std::endl;
  }
 debug = biassgd_vertex_program::debug;
  //  omp_set_num_threads(clopts.get_ncpus());
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
  engine_type engine(dc, graph, exec_type, clopts);

  // Add error reporting to the engine
  const bool success = engine.add_edge_aggregator<error_aggregator>
    ("error", error_aggregator::map, error_aggregator::finalize) &&
    engine.aggregate_periodic("error", interval);
  ASSERT_TRUE(success);
  

  biassgd_vertex_program::GLOBAL_MEAN = graph.map_reduce_edges<double>(calc_global_mean);
  biassgd_vertex_program::NUM_TRAINING_EDGES = graph.map_reduce_edges<size_t>(count_edges);
  biassgd_vertex_program::GLOBAL_MEAN /= biassgd_vertex_program::NUM_TRAINING_EDGES;
  dc.cout() << "Global mean is: " <<biassgd_vertex_program::GLOBAL_MEAN << std::endl;

  // Signal all vertices on the vertices on the left (libersgd) 
  engine.map_reduce_vertices<graphlab::empty>(biassgd_vertex_program::signal_left);
 

  // Run the PageRank ---------------------------------------------------------
  dc.cout() << "Running Bias-SGD" << std::endl;
  dc.cout() << "(C) Code by Danny Bickson, CMU " << std::endl;
  dc.cout() << "Please send bug reports to danny.bickson@gmail.com" << std::endl;
  dc.cout() << "Time   Training    Validation" <<std::endl;
  dc.cout() << "       RMSE        RMSE " <<std::endl;
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



