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


/*
 *  PageRank tutorial application.
 *  pagerank.cpp
 *
 *  For description of the PageRank algorithm, see Wikipedia article
 *  http://en.wikipedia.org/wiki/Pagerank
 */



#ifndef GRAPHLAB_PAGERANK_HPP
#define GRAPHLAB_PAGERANK_HPP

#include <string>
#include <graphlab.hpp>


/// Declare Constants --------------------------------------------------------->

// Constants for the algorithm. Better way would be to
// pass them in the shared_data to the update function, but for
// the simplicity of this example, we simply define them here.

// The termination bound
extern double termination_bound; // Defined in pagerank.cpp
// PageRank random reset probability
extern double random_reset_prob; // Defined in pagerank.cpp

/// Define Core Data Types---------------------------------------------------->
/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  float weight;
  float old_source_value;
  edge_data(float weight = 1) :
    weight(weight), old_source_value(0) { } 
}; // End of edge data


#ifdef DIFFABLE
/**
 * Stores the value and the self weight
 */
struct vertex_data : public graphlab::idiffable<vertex_data> {
  float value;
  float self_weight; // GraphLab does not support edges from vertex to itself, so
  // we save weight of vertex's self-edge in the vertex data
  vertex_data(float value = 1) : value(value), self_weight(0) { }
  void apply_diff(const vertex_data& changed, const vertex_data& old) {
    value += (changed.value - old.value);
  }
}; // End of vertex data
#else
/**
 * Stores the value and the self weight
 */
struct vertex_data {
  float value;
  float self_weight; // GraphLab does not support edges from vertex to itself, so
  // we save weight of vertex's self-edge in the vertex data
  vertex_data(float value = 1) : value(value), self_weight(0) { }
}; // End of vertex data
#endif

std::ostream& operator<<(std::ostream& out, const edge_data& edata);
std::ostream& operator<<(std::ostream& out, const vertex_data& vdata);



//! The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> pagerank_graph;

//! Predeclar the pagerank update functor
class pagerank_update;

/**
 * The collection of graphlab types restricted to the graph type used
 * in this program.
 */
typedef graphlab::types<pagerank_graph, pagerank_update> gl;


//! Save the graph to tsv file
void save_edges_as_tsv(const std::string& fname, 
                       const pagerank_graph& graph);

//! save the pagerank as a tsv file
void save_pagerank(const std::string& fname,
                   const pagerank_graph& graph);

/**
 * Load a graph in metis format
 * line  0: 5 7
 * line  1: 1 2
 * line  2: 3
 *         ...
 * line  5: 5
 *
 */
bool load_graph_from_metis_file(const std::string& filename,
                                pagerank_graph& graph);

bool load_graph_from_jure_file(const std::string& filename,
                               pagerank_graph& graph);


/**
 * Predecleration of the graph file loading function.  Defined at
 * bottom of file for clarity.
 *
 * Load a graph file specified in the format:
 *
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph_from_tsv_file(const std::string& filename,
                              pagerank_graph& graph);

/**
 * Makes a small to graph.
 */
void make_toy_graph(pagerank_graph& graph);



#endif
