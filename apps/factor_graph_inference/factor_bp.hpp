#ifndef FACTOR_BP_HPP
#define FACTOR_BP_HPP

#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>


#include <graphlab.hpp>

#include "factor_graph.hpp"


#include <graphlab/macros_def.hpp>

enum constants {BOUND_ID, DAMPING_ID};

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



// Update Functions
// ===========================================================================>


void bp_update(gl_types::iscope& scope, 
               gl_types::icallback& scheduler,
               gl_types::ishared_data* shared_data) {
  assert(shared_data != NULL);
  // Get the shared data
  double bound =   shared_data->get_constant(BOUND_ID).as<double>();
  double damping = shared_data->get_constant(DAMPING_ID).as<double>();

  // Grab the state from the scope
  // ------------------------------------------------------->
  // Get the vertex data
  vertex_data& v_data = scope.vertex_data();
  
  // Get the in and out edges by reference
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();
  assert(in_edges.size() == out_edges.size());

  // Flip the old and new messages to improve safety when using the
  // unsynch scope
  foreach(graphlab::edge_id_t ineid, in_edges) {   
    // Get the in and out edge data
    edge_data& in_edge = scope.edge_data(ineid);
    // copy the message onto the out edge
    in_edge.old_message = in_edge.message;
  }

  // Compute the belief
  // ------------------------------------------------------->
  // Initialize the belief as the value of the factor
  v_data.belief = v_data.potential;
  foreach(graphlab::edge_id_t ineid, in_edges) {
    // Get the message
    const edge_data& e_data = scope.edge_data(ineid);
    v_data.belief *= e_data.old_message;
  }
  v_data.belief.normalize();


  // Compute outbound messages
  // ------------------------------------------------>
  // Send outbound messages
  factor_type cavity, tmp_msg;
  for(size_t i = 0; i < in_edges.size(); ++i) {
    // Get the edge ids
    graphlab::edge_id_t outeid = out_edges[i];
    graphlab::edge_id_t ineid = in_edges[i];
    // CLEVER HACK: Here we are expoiting the sorting of the edge ids
    // to do fast O(1) time edge reversal
    assert(scope.target(outeid) == scope.source(ineid));
    // Get the in and out edge data
    const edge_data& in_edge = scope.edge_data(ineid);
    edge_data& out_edge = scope.edge_data(outeid);
    
    // Create stack allocated cavity factor
    cavity = v_data.belief;   
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



void make_bp_graph(const factor_graph& fgraph,
                   gl_types::graph& graph) {
  assert(!fgraph.variables().empty());
  assert(!fgraph.factors().empty());
  assert(fgraph.variables().rbegin()->id
         == fgraph.variables().size()-1);
  vertex_data vdata;
  // Add all the variables 
  foreach(factor_graph::variable_type variable, fgraph.variables()) {
    factor_graph::domain_type domain(variable);
    vdata.potential.set_args(domain);
    vdata.potential.uniform();
    vdata.belief = vdata.potential;  
    graphlab::vertex_id_t vid = graph.add_vertex(vdata);
    assert(vid == variable.id);
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
      graph.add_edge(factor_index, variable.id, edata);                     
      graph.add_edge(variable.id, factor_index, edata);
    }
    ++factor_index;
  }
}


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



#include <graphlab/macros_undef.hpp>
#endif
