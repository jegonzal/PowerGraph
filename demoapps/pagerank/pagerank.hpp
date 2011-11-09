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


/// Types------------------------------------------------------------------>

/**
 * Stores the value and the self weight
 */
struct vertex_data {
  float value, old_value, self_weight; 
  vertex_data(float value = 1) : 
    value(value), old_value(0), self_weight(0) { }
}; // End of vertex data
//! Print the vertex data
std::ostream& operator<<(std::ostream& out, const vertex_data& vdata);


/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  float weight;
  edge_data(float weight = 1) : weight(weight) { } 
}; // End of edge data
//! Print the edge data
std::ostream& operator<<(std::ostream& out, const edge_data& edata);


//! The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> graph_type;



/// Utility routines defined in utility.cpp ------------------------------->

//! Save the graph to tsv file
void save_graph_as_edge_list(const std::string& fname, 
                             const graph_type& graph);

//! save the pagerank vector as a tsv file
void save_pagerank(const std::string& fname,
                   const graph_type& graph);

//! Return the ids of the top k pages
void get_top_pages(const graph_type& graph, size_t num_pages,
                   std::vector<graph_type::vertex_id_type>& ret);



//! Load the graph from a file with a given format
bool load_graph(const std::string& filename,
                const std::string& format,
                graph_type& graph);


//! Load the graph from a metis (adjacency format file)
bool load_graph_from_metis_file(const std::string& filename,
                                graph_type& graph);
//! Load the graph in Jure Leskovec's file format
bool load_graph_from_jure_file(const std::string& filename,
                               graph_type& graph);
//! Load the graph from a tab separated file
bool load_graph_from_tsv_file(const std::string& filename,
                              graph_type& graph);
//! Make a toy graph to quickly test pagerank
void make_toy_graph(graph_type& graph);



#endif
