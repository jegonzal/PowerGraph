#ifndef SEQUENTIAL_GIBBS_HPP
#define SEQUENTIAL_GIBBS_HPP

#include "data_structures.hpp"

// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>

void update_vertex(gl::graph& graph,
                   const binary_factor& edge_factor, 
                   vertex_id_t vertex) {
  vertex_data& vdata = graph.vertex_data(vertex);
  unary_factor conditional = vdata.potential;
  foreach(edge_id_t eid, graph.in_edge_ids(vertex)) {
    vertex_id_t source = graph.source(eid);
    const vertex_data& neighbor = 
      graph.vertex_data(source);
    conditional.condition(edge_factor, neighbor.asg);
  }
  conditional.normalize();
  size_t sample = conditional.sample();
  vdata.asg = sample;
  vdata.updates++;
  vdata.belief.plus(conditional);
} // end of update vertex


void forward_sweep(gl::graph& graph,
                   const binary_factor& edge_factor) {
  for(size_t i = 0; i < graph.num_vertices(); ++i) {
    update_vertex(graph, edge_factor, i);
  }
}


void forward_sweep(gl::graph& graph,
                   const binary_factor& edge_factor,
                   size_t nsamples) {
  for(size_t i = 0; i < nsamples; ++i) {
    update_vertex(graph, edge_factor,
                  i % graph.num_vertices() );
  }
}


void random_sweep(gl::graph& graph,
                   const binary_factor& edge_factor,
                   size_t nsamples) {
  for(size_t i = 0; i < nsamples; ++i) {
    vertex_id_t vid = graphlab::random::rand_int(graph.num_vertices() - 1);

    update_vertex(graph, edge_factor, vid);
                  
  }
}



// Include the macro for the foreach operation
#include <graphlab/macros_undef.hpp>
#endif 
