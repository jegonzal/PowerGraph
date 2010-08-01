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

#include "sequential_tree_gibbs.hpp"

// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>




void draw_beliefs(const graph_type& graph,
                  const std::string& filename,
                  size_t rows, size_t cols,
                  bool expectation = true) {
  image img(rows, cols);
  graphlab::unary_factor marginal;
  for(size_t v = 0; v < graph.num_vertices(); ++v) {
    marginal = graph.vertex_data(v).belief;
    marginal.normalize();
    double value = 0;
    if(expectation) {
      value = marginal.expectation();
    } else {
      value = marginal.max_asg();
    }
    img.pixel(v) = value;
  } 
  img.save(filename.c_str());
} // End of draw_beliefs


void draw_asg(const graph_type& graph,
              const std::string& filename,
              size_t rows, size_t cols) {
  image img(rows, cols);
  for(size_t v = 0; v < graph.num_vertices(); ++v) 
    img.pixel(v)  = graph.vertex_data(v).asg;
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



void draw_tree(std::vector<vertex_id_t>& tree_order,
               const std::string& name,
               size_t rows, size_t cols) { 
  image img(rows, cols);
  for(size_t i = 0; i < img.pixels(); ++i) 
    img.pixel(i) = 0;
  for(size_t i = 0; i < tree_order.size(); ++i) 
    img.pixel(tree_order[i]) = (tree_order.size() - i) + 1;
  img.save(name.c_str());
}


void draw_tree(graph_type& graph,
               const std::string& name,
               size_t rows, size_t cols) { 
  image img(rows, cols);
  for(size_t i = 0; i < img.pixels(); ++i) {
    vertex_data& vdata = graph.vertex_data(i);
    if(vdata.parent != NULL_VID) { // in tree
      //      img.pixel(i) = std::log( vdata.height + 2) ;
      img.pixel(i) = vdata.height + 10;
      if(vdata.parent == i) img.pixel(i) += 10;
    } else  // not in tree
      img.pixel(i) = 0;
  }
  img.save(name.c_str());
}



void draw_trees(gl::graph& graph,
                size_t max_tree_size,
                size_t count,
                size_t rows, size_t cols) { 

  for(size_t i = 0; i < count; ++i) {
    vertex_id_t root = graphlab::random::rand_int(graph.num_vertices() - 1);
    std::cout << "Drawing tree " << i
              << " out of " << count
              << " at vertex " << root
              << ". " << std::endl;

    std::vector<vertex_id_t> tree_order;
    // grow_bfs_tree(graph, root, max_tree_size, tree_order);
    grow_dfs_tree(graph, root, max_tree_size, tree_order);    
    std::stringstream strm;
    strm << "tree_struct_" << i << ".pgm";
    draw_tree(tree_order, strm.str(), rows, cols);
  }
}












#include <graphlab/macros_undef.hpp>
#endif
