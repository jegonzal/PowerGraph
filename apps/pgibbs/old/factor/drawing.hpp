#ifndef DRAWING_HPP
#define DRAWING_HPP


#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>




#include "data_structures.hpp"



// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>




void draw_beliefs(const graph_type& graph,
                  const std::string& filename,
                  size_t rows, size_t cols,
                  bool expectation = true) {
  image img(rows, cols);
  factor_t marginal;
  std::vector<double> exp_values;
  for(size_t v = 0; v < graph.num_vertices(); ++v) {
    marginal = graph.vertex_data(v).belief;
    marginal.normalize();
    double value = 0;
    if(expectation) {
      exp_values.clear();
      marginal.expectation(exp_values);
      assert(exp_values.size() == 1);
      value = exp_values[0];
    } else {
      assignment_t asg = marginal.max_asg();
      assert(asg.num_vars() == 1);
      value = asg.asg(v);
    }
    img.pixel(v) = value;
  } 
  img.save(filename.c_str());
} // End of draw_beliefs


void draw_asg(const graph_type& graph,
              const std::string& filename,
              size_t rows, size_t cols) {
  image img(rows, cols);
  for(size_t v = 0; v < graph.num_vertices(); ++v)  {
    assignment_t asg = graph.vertex_data(v).asg;
    assert(asg.num_vars() == 1);
    img.pixel(v) = asg.asg(v); 

  }
  img.save(filename.c_str());
} // End of draw_beliefs




void draw_edge_weights(gl::graph& graph,
                       const std::string& filename,
                       size_t rows, size_t cols) {
  image img(rows, cols);

  for(size_t i = 0; i < rows*cols; ++i) 
    img.pixel(i) = 0;

  for(graphlab::edge_id_t eid = 0; eid < graph.num_edges(); ++eid) {
    graphlab::vertex_id_t source = graph.source(eid);
    graphlab::vertex_id_t target = graph.target(eid);
    img.pixel(source) =
      std::max(img.pixel(source), graph.edge_data(eid).weight);
    img.pixel(target) =
      std::max(img.pixel(target), graph.edge_data(eid).weight);
  }
  img.save(filename.c_str());
}




void draw_tree(graph_type& graph,
               const std::string& name,
               size_t rows, size_t cols) { 
  image img(rows, cols);
  for(size_t i = 0; i < img.pixels(); ++i) {
    vertex_data& vdata = graph.vertex_data(i);
    if(vdata.parent != NULL_VID) { // in tree
      //      img.pixel(i) = std::log( vdata.height + 2) ;
      img.pixel(i) = vdata.height + 100;
      if(vdata.parent == i) img.pixel(i) *= 2;
    } else  // not in tree
      img.pixel(i) = 0;
  }
  img.save(name.c_str());
}


void draw_state(graph_type& graph,
                const std::string& name,
                size_t rows, size_t cols) { 
  image img(rows, cols);
  for(size_t i = 0; i < img.pixels(); ++i) {
    vertex_data& vdata = graph.vertex_data(i);
    img.pixel(i) = vdata.state;
  }
  img.save(name.c_str());
}













#include <graphlab/macros_undef.hpp>
#endif
