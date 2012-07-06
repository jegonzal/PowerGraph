#include <cassert>

#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <graphlab.hpp>

#include "factor_graph.hpp"


#include <graphlab/macros_def.hpp>



// Types from the factor_graph.hpp
typedef factor_graph::factor_type factor_type;
typedef factor_graph::domain_type domain_type;
typedef factor_type::assignment_type  assignment_type;
typedef factor_graph::variable_type   variable_type;

// Structs
// ===========================================================================>

struct vertex_data {
  factor_type potential;
  factor_type belief;
  void load(graphlab::iarchive& arc) {
    arc >> potential;
    arc >> belief;
  }
  void save(graphlab::oarchive& arc) const {
    arc << potential;
    arc << belief;
  }
}; // End of vertex data

struct edge_data {
  factor_type message;
  factor_type old_message;
  void load(graphlab::iarchive& arc) {
    arc >> message;
    arc >> old_message;
  }
  void save(graphlab::oarchive& arc) const {
    arc << message;
    arc << old_message;
  }
}; // End of edge data


// Define the graph types
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;

double bound = 1e-5;
double damping = 0.3;


// Update Functions
// ===========================================================================>
void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler) {
   // Grab the state from the scope
  // ------------------------------------------------------->
  // Get the vertex data
  vertex_data& vdata = scope.vertex_data();
  
  // Get the in and out edges by reference
  const gl_types::edge_list in_edges = scope.in_edge_ids();
  const gl_types::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size());

  // Compute the belief
  // ------------------------------------------------------->
  // Initialize the belief as the value of the factor
  vdata.belief = vdata.potential;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the in and out edge data
    edge_data& edata = scope.edge_data(ineid);
    edata.old_message = edata.message;
    vdata.belief *= edata.old_message;
  }
  vdata.belief.normalize();


  // Compute outbound messages
  // ------------------------------------------------>
  // Send outbound messages
  factor_type cavity, tmp_msg;
  for(size_t i = 0; i < in_edges.size(); ++i) {
    // Get the edge ids
    const graphlab::edge_id_t outeid = out_edges[i];
    const graphlab::edge_id_t ineid = in_edges[i];
    // CLEVER HACK: Here we are expoiting the sorting of the edge ids
    // to do fast O(1) time edge reversal
    assert(scope.target(outeid) == scope.source(ineid));
    // Get the in and out edge data
    const edge_data& in_edge = scope.edge_data(ineid);
    edge_data& out_edge = scope.edge_data(outeid);
    
    // Create stack allocated cavity factor
    cavity = vdata.belief;   
    // Compute cavity
    cavity /= in_edge.old_message; // Make the cavity a cavity
    cavity.normalize();

    // Create stack allocated temporary message
    tmp_msg.set_args(out_edge.message.args());
    assert(tmp_msg.num_vars() == 1);
    tmp_msg.marginalize(cavity);
    tmp_msg.normalize();

    // Compute message residual
    double residual = tmp_msg.l1_diff(out_edge.old_message);
     // Damp the message
    tmp_msg.damp(out_edge.message, damping);   
    // Assign the out message
    out_edge.message = tmp_msg;
    
    if(residual > bound) {
      gl_types::update_task task(scope.target(outeid), bp_update);      
      scheduler.add_task(task, residual);
    }    
  }
} // end of BP_update


/**
 * Construct a belief propagation graph from a factor graph
 */
void make_bp_graph(const factor_graph& fgraph,
                   gl_types::graph& graph) {
  assert(!fgraph.variables().empty());
  assert(!fgraph.factors().empty());
  assert(fgraph.variables().rbegin()->id()
         == fgraph.variables().size()-1);
  vertex_data vdata;
  // Add all the variables 
  foreach(factor_graph::variable_type variable, fgraph.variables()) {
    factor_graph::domain_type domain(variable);
    vdata.potential.set_args(domain);
    vdata.potential.uniform();
    vdata.belief = vdata.potential;  
    graphlab::vertex_id_t vid = graph.add_vertex(vdata);
    assert(vid == variable.id());
  }
  assert(graph.num_vertices() == fgraph.variables().size());
  // Add all the factors and all the edges
  size_t factor_index = graph.num_vertices();
  edge_data edata;
  foreach(const factor_graph::factor_type& factor, fgraph.factors()) {
    // Setup the vertex data for a factor
    vdata.potential = factor;
    vdata.belief.set_args(factor.args());
    vdata.belief.uniform();
    // Add the factor to the graph
    graphlab::vertex_id_t vid = graph.add_vertex(vdata);
    assert(vid == factor_index);    
    // Attach all the edges
    for(size_t i = 0; i < factor.num_vars(); ++i) {
      factor_graph::variable_type variable = factor.args().var(i);
      factor_graph::domain_type domain(variable);
      edata.message.set_args(domain);
      edata.message.uniform();
      edata.old_message = edata.message;     
      graph.add_edge(factor_index, variable.id(), edata);                     
      graph.add_edge(variable.id(), factor_index, edata);
    }
    ++factor_index;
  }
  graph.finalize();
} // end of make_bp_graph


/**
 * Save the final belief estimates to a text file.
 */
void save_beliefs(const std::string& filename,
                  factor_graph& fgraph,
                  gl_types::graph& graph,
                  size_t num_variables) {
  // Open the file to store the final belief vectors
  std::ofstream fout(filename.c_str());
  assert(fout.good());
  // Save all the beliefs
  for(size_t i = 0; i < num_variables; ++i) {
    vertex_data& vdata = graph.vertex_data(i);
    vdata.belief.normalize();
    fout << fgraph.var_name(i) << " // ";
    // Save the normalized belief estimates to the file
    for(size_t j = 0; j < vdata.belief.size(); ++j) {
      fout << std::exp(vdata.belief.logP(j));
      if((j + 1) < vdata.belief.size()) fout << ", ";
    }
    fout << "\n";
  }
  fout.close();   
} // end of save beliefs




// MAIN
// ============================================================================>
int main(int argc, char** argv) {
  std::cout << "This program solves the sidechain prediction task."
            << std::endl;

  // Parse command line arguments --------------------------------------------->
  bound = 1E-5; // <-- Defined globally
  damping = 0.3; // <-- Defined globally
  std::string network_filename;
  std::string beliefs_filename = "beliefs.txt";
 
  
  graphlab::command_line_options clopts("Run Loopy BP on an Alchemy Network");
  clopts.attach_option("graph",
                       &network_filename,
                       "The Alchemy factor graph file.");
  clopts.add_positional("graph");
  clopts.attach_option("bound",
                       &bound, bound,
                       "Residual termination bound");
  clopts.attach_option("damping",
                       &damping, damping,
                       "The amount of message damping (higher = more damping)");
  clopts.attach_option("beliefs",
                       &beliefs_filename, beliefs_filename,
                       "The file to save the belief predictions"); 
  clopts.set_scheduler_type("splash(splash_size=100)");
  clopts.set_scope_type("edge");

  // set the global logger
  // global_logger().set_log_level(LOG_WARNING);
  // global_logger().set_log_to_console(true);

  bool success = clopts.parse(argc, argv);
  if(!success && !clopts.is_set("graph")) {    
    return EXIT_FAILURE;
  }
 
  gl_types::core core;
  core.set_engine_options(clopts);
  core.sched_options().add_option("update_function", bp_update);

  
  // Load the factor graph from file ------------------------------------------>
  std::cout << "Loading Factor Graph in Alchemy Format" << std::endl;
  factor_graph fgraph;
  fgraph.load_alchemy(network_filename);
  const size_t num_variables = fgraph.variables().size();
  const size_t num_factors = fgraph.factors().size();
  std::cout << "Finished!" << std::endl;

  
  // Build the BP graph from the factor graph---------------------------------->
  std::cout << "Building BP graph from the factor graph" << std::endl;
  make_bp_graph( fgraph, core.graph() ); 
  size_t num_vertices = core.graph().num_vertices();
  assert(num_vertices == num_variables + num_factors);
  size_t num_edges = core.graph().num_edges();
  std::cout << "Loaded: " << num_vertices << " vertices "
            << "and " << num_edges << " edges." << std::endl;
  std::cout << "Finished!" << std::endl;

  
  // Tell the scheduler that the bp_update function should be applied
  // to all vertices with priority:
  double initial_priority = 100.0;
  core.add_task_to_all(bp_update, initial_priority);
  

  // Running the engine ------------------------------------------------------->
  std::cout << "Running the engine. " << std::endl;
  // Run the engine (this blocks until their are no tasks left).
  // Convergence is determined when there are no update tasks left.
  const double runtime = core.start();
  std::cout << "Done!" << std::endl;
  
  // Print some fun facts------------------------------------------------------>
  size_t update_count = core.last_update_count();
  std::cout << "Finished Running engine in " << runtime 
            << " seconds." << std::endl
            << "Total updates: " << update_count << std::endl
            << "Efficiency: " << (double(update_count) / runtime)
            << " updates per second "
            << std::endl;

  
  
  // Save the beliefs --------------------------------------------------------->
  std::cout << "Saving the beliefs. " << std::endl;
  save_beliefs(beliefs_filename, fgraph, core.graph(), num_variables);
  std::cout << "Finished saving beliefs." << std::endl;


  // // Save statistics about the run -------------------------------------------->
  // std::cout << "Saving the statistics file. " << std::endl;
  // std::ofstream stats(opts.stats_filename.c_str(), std::ios::app);
  // assert(stats.good());
  // stats << opts.scheduler << ", "
  //       << opts.network_filename << ", "
  //       << opts.scope << ", "
  //       << opts.ncpus << ", "
  //       << runtime << ", "
  //       << update_count << ", "
  //       << num_vertices << ", "
  //       << num_edges << ", "
  //       << num_variables << ", "
  //       << num_factors << ", "
  //       << opts.splash_size << std::endl;
  // stats.close();
  // std::cout << "Finished saving statistics file." << std::endl;

 
  return EXIT_SUCCESS;
} // end of main

