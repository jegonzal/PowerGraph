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


/// Declare Constants
/// --------------------------------------------------------->

// Constants for the algorithm. Better way would be to
// pass them in the shared_data to the update function, but for
// the simplicity of this example, we simply define them here.

// The termination bound
extern double termination_bound; // Defined in pagerank.cpp

// PageRank random reset probability
extern double random_reset_prob; // Defined in pagerank.cpp

/// Define Core Data
/// Types---------------------------------------------------->
/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  float weight;
  edge_data(float weight = 1) : weight(weight) { } 
}; // End of edge data

/**
 * Stores the value and the self weight
 */
struct vertex_data {
  float value, old_value, self_weight; 
  vertex_data(float value = 1) : 
    value(value), old_value(0), self_weight(0) { }
}; // End of vertex data


std::ostream& operator<<(std::ostream& out, const edge_data& edata);
std::ostream& operator<<(std::ostream& out, const vertex_data& vdata);

//! The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> graph_type;
//! Save the graph to tsv file
void save_edges_as_tsv(const std::string& fname, 
                       const graph_type& graph);
//! save the pagerank as a tsv file
void save_pagerank(const std::string& fname,
                   const graph_type& graph);

bool load_graph_from_metis_file(const std::string& filename,
                                graph_type& graph);
bool load_graph_from_jure_file(const std::string& filename,
                               graph_type& graph);
bool load_graph_from_tsv_file(const std::string& filename,
                              graph_type& graph);
void make_toy_graph(graph_type& graph);



#endif
