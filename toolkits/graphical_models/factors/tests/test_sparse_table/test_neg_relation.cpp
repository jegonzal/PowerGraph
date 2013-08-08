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
 * a factor node to constrain a negation relation a = -b
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
typedef graphlab::dense_table<MAX_DIM>          dense_table_t;
typedef graphlab::sparse_table<MAX_DIM>         sparse_table_t;
typedef graphlab::discrete_assignment<MAX_DIM>  assignment_t;
typedef graphlab::discrete_domain<MAX_DIM>      domain_t;
typedef graphlab::discrete_variable             variable_t;


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
std::vector<double> compute_labels(size_t n_labels, 
    double max_range, double min_range);
void compute_normal_dist(double mean, double std_dev, 
    const std::vector<double>& labels, dense_table_t& prior);


// MAIN
// ============================================================================>
int main(int argc, char** argv) {
  std::cout << "This program solves the sum task."
            << std::endl;

  global_logger().set_log_level(LOG_DEBUG);

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

  // Create the variables
  size_t n_labels = 10; 
  double cost_scale = 100;

  double min_range = -3.0; double max_range = 3.0; 
  std::vector<double> labels = compute_labels(n_labels, max_range, min_range);

  variable_t var_a = fgraph.add_variable(n_labels, "var_a");
  {
  double mean = 2.3; double std_dev = 0.1; 
  dense_table_t& prior = fgraph.prior_for_variable(var_a);
  compute_normal_dist(mean, std_dev, labels, prior);
  std::cout << "var_a_prior=" << prior << std::endl;
  }

  variable_t var_b = fgraph.add_variable(n_labels, "var_b");
  {
  dense_table_t& prior = fgraph.prior_for_variable(var_b);
  prior.zero();
  std::cout << "var_b_prior=" << prior << std::endl;
  }


  // Create a factor
  std::vector<variable_t> args;
  // connect vertical neighbors
  args.push_back(var_a);
  args.push_back(var_b);
  // Build the factor
  domain_t dom(args);
  sparse_table_t neg(dom);
  // Set the weights
  assignment_t asg(dom);
  for(size_t i=0; i<n_labels; ++i) { 
    asg.set_asg(var_a, i);
    asg.set_asg(var_b, n_labels-(i+1));
    float err = 0;
    neg.set_logP(asg, -1*(cost_scale*err*err));
  }
  // Save the factor to the factor graph
  fgraph.add_factor(neg, "neg");
  std::cout << neg << std::endl;


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
  fgraph.print_variable(var_a, labels);
  fgraph.print_variable(var_b, labels);

  double a = labels[fgraph.belief_for_variable(var_a).max_index()];
  double b = labels[fgraph.belief_for_variable(var_b).max_index()];
  std::cout << "var_a: " << a << std::endl;
  std::cout << "var_b: " << b << std::endl;

  double b_prime = -1*a;
  double err = std::abs(b - b_prime);
  std::cout << "b: " << b << " b_prime: " << b_prime << " err: " << err << std::endl;
  ASSERT_LT(err, 1E-4);
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

std::vector<double> compute_labels(size_t n_labels, 
    double max_range, double min_range) 
{
  std::vector<double> labels(n_labels, 0.0);

  double step = (max_range - min_range)/(n_labels-1);
  for(unsigned i = 0; i < n_labels; ++i) {
    labels[i] = min_range + i*step;
  }
  return labels;
}

void compute_normal_dist(double mean, double std_dev, 
    const std::vector<double>& labels, dense_table_t& prior) 
{
  domain_t::const_iterator asg = prior.domain().begin();
  domain_t::const_iterator end = prior.domain().end();
  std::vector<double>::const_iterator label = labels.begin();
  for( ; asg != end; ++asg, ++label) {
    double nv = (*label-mean)/std_dev;
    prior.set_logP( *asg, -1*nv*nv );
  }
}
