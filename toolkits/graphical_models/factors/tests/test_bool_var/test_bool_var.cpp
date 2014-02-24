/**  
 *  Software submitted by 
 *  Systems & Technology Research / Vision Systems Inc., 2013
 *
 *  Approved for public release; distribution is unlimited. [DISTAR Case #21428]
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
 * This file contains an example of graphlab belief propagation on
 * a factor graph designed to test a unary factor node.
 *
 *  \author Scott Richardson 
 */

// INCLUDES ===================================================================>

// Including Standard Libraries


#include <cstdlib>
#include <cassert>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>

#include <graphlab.hpp>

#include <factors/factor_graph.hpp>
#include <factors/bp_vertex_program.hpp>


// Include the macro for each operation
#include <graphlab/macros_def.hpp>


const size_t MAX_DIM = 4;
typedef graphlab::dense_table<MAX_DIM>   dense_table_t;
typedef graphlab::discrete_variable      variable_t;


struct clopts_vals { 
  clopts_vals(double bound = 1E-4, double damping = 0.0, std::string exec_t="sync") : 
      BOUND(bound), DAMPING(damping), exec_type(exec_t) { }

  double BOUND;
  double DAMPING;
  std::string exec_type;
};

int setup_cli(graphlab::command_line_options& clopts, clopts_vals& clvals, 
    int argc, char** argv);
template<size_t MAX_DIM>
void run_engine(graphlab::distributed_control& dc, 
    typename belief_prop::graph_type<MAX_DIM>::type& graph, 
    const std::string& exec_type, const graphlab::command_line_options& clopts);

// MAIN
// ============================================================================>
int main(int argc, char** argv) {
  std::cout << "This program solves the sum task."
            << std::endl;

  graphlab::mpi_tools::init(argc, argv);
  ///! Create a distributed control object (must come after mpi_tools::init())
  graphlab::distributed_control dc; 

  // Parse command line arguments --------------------------------------------->
  graphlab::command_line_options clopts("Run Loopy BP on a Network");
  clopts_vals clvals;
  if( setup_cli(clopts, clvals, argc, argv) != EXIT_SUCCESS ) return EXIT_FAILURE;

  ///! Create a distributed graph object 
  belief_prop::graph_type<MAX_DIM>::type graph(dc, clopts);



  // Create the factor graph ------------------------------------------>
  std::cout << "Loading Factor Graph" << std::endl;
  belief_prop::factor_graph<MAX_DIM> fgraph;

  // Create the variable
  size_t nlabels = 2;
  variable_t bool_var_b = fgraph.add_variable(nlabels, "bool_var_b");
  dense_table_t& prior = fgraph.prior_for_variable(bool_var_b);
  prior.zero();


  // Set the weights
  std::vector<double> logf(2); logf[0] = std::log(0.6); logf[1] = std::log(0.4);
  // Build a unary factor
  dense_table_t bool_obs( bool_var_b, logf );
  // Save the factor to the factor graph
  fgraph.add_factor(bool_obs, "bool_obs");


  const size_t num_variables = fgraph.num_variables();
  const size_t num_factors = fgraph.num_factors();
  std::cout << "num_variables: " << num_variables << " num_factors: " << num_factors << std::endl;
  std::cout << "Finished!" << std::endl;


  // Build the BP graph from the factor graph---------------------------------->
  std::cout << "Building BP graph from the factor graph" << std::endl;
  fgraph.make_bp_graph( graph, clvals.BOUND, clvals.DAMPING ); 
  run_engine<MAX_DIM>(dc, graph, clvals.exec_type, clopts);
  fgraph.pull_beliefs_for_variables( graph );


  // Saving the output -------------------------------------------------------->
  std::cout << fgraph.belief_for_variable(bool_var_b) << std::endl;
  double bobs = fgraph.belief_for_variable(bool_var_b).logP(1);
  double err = abs(bobs - .405465);
  ASSERT_LT(err, .01);
  std::cout << "All tests passed" << std::endl;
} // end of main


// UTILS
// ============================================================================>
int setup_cli(graphlab::command_line_options& clopts, clopts_vals& clvals,
    int argc, char** argv) {

  clopts.attach_option("bound", clvals.BOUND,
                       "Residual termination bound");
  clopts.attach_option("damping", clvals.DAMPING,
                       "The amount of message damping (higher = more damping)");
//  clopts.attach_option("beliefs", &beliefs_filename,
//                       "The file to save the belief predictions"); 
  clopts.attach_option("engine", clvals.exec_type,
                       "The type of engine to use {async, sync}.");
  clopts.set_scheduler_type("fifo");

  bool success = clopts.parse(argc, argv);
  if(!success) {    
    std::cout << "Error parsing command line arguments!"
              << std::endl;
    graphlab::mpi_tools::finalize();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template<size_t MAX_DIM>
void run_engine(graphlab::distributed_control& dc, 
    typename belief_prop::graph_type<MAX_DIM>::type& graph, 
    const std::string& exec_type, 
    const graphlab::command_line_options& clopts) 
{
  size_t num_vertices = graph.num_vertices();
  size_t num_edges = graph.num_edges();
  std::cout << "Loaded: " << num_vertices << " vertices "
            << "and " << num_edges << " edges." << std::endl;
  std::cout << "Finished!" << std::endl;

  // Create the engine -------------------------------------------------------->
  std::cout << "Creating the engine. " << std::endl;
  typedef graphlab::omni_engine<belief_prop::bp_vertex_program<MAX_DIM> > engine_type;
  engine_type engine(dc, graph, exec_type, clopts);

  std::cout << "Scheduling all vertices" << std::endl;
  engine.signal_all();
  std::cout << "Starting the engine" << std::endl;
  engine.start();
  const float runtime = engine.elapsed_seconds();
  size_t update_count = engine.num_updates();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;
}

