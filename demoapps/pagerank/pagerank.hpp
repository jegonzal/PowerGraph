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

#include <vector>
#include <string>

#include <graphlab.hpp>


/// Types------------------------------------------------------------------>


/**
 * Stores the value and the self weight
 */
struct vertex_data {
  uint32_t nupdates;
  double value, old_value;
  vertex_data(double value = 1) : 
    nupdates(0), value(value), old_value(0) { }
}; // End of vertex data
SERIALIZABLE_POD(vertex_data);

//! Print the vertex data
std::ostream& operator<<(std::ostream& out, const vertex_data& vdata);

/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  double weight;
  edge_data(double weight = 1) : weight(weight) { } 
}; // End of edge data
SERIALIZABLE_POD(edge_data);

//! Print the edge data
std::ostream& operator<<(std::ostream& out, const edge_data& edata);

//! The type of graph used in this program
typedef graphlab::graph2<vertex_data, edge_data> graph_type;
//typedef graphlab::graph<vertex_data, edge_data> graph_type;







/// Utility routines defined in utility.cpp ------------------------------->
std::ostream& operator<<(std::ostream& out, const edge_data& edata) {
  return out << "E(w: " << edata.weight << ")";
}

std::ostream& operator<<(std::ostream& out, const vertex_data& vdata) {
  return out << "Rank=" << vdata.value;
}



void save_pagerank(const std::string& fname,
                   const graph_type& graph) {
  std::ofstream fout;
  fout.open(fname.c_str());
  fout << std::setprecision(10);
  for(graph_type::vertex_id_type vid = 0; 
      vid < graph.num_vertices(); ++vid) {
    fout << graph.vertex_data(vid).value << "\n";
  }
  fout.close();
} // end of save_pagerank


void get_top_pages(const graph_type& graph, size_t num_pages,
                   std::vector<graph_type::vertex_id_type>& ret) {
  typedef std::pair<float, graph_type::vertex_id_type> pair_type;
  std::priority_queue<pair_type> top;
  for(graph_type::vertex_id_type vid = 0; vid < graph.num_vertices(); ++vid) {
    const graph_type::vertex_data_type& vdata = graph.vertex_data(vid);
    top.push(std::make_pair(-vdata.value, vid));
    if(top.size() > num_pages) top.pop();
  }
  if(top.empty()) return;
  ret.resize(top.size());
  for(size_t i = top.size()-1; i < top.size(); --i) {
    ret[i] = top.top().second;
    top.pop();
  }
} // end of top pages


void normalize_graph(graph_type& graph) {
  logstream(LOG_INFO)
    << "Optimizing graph layout in memory." << std::endl;
  graph.finalize();
  logstream(LOG_INFO)
    << "Renormalizing transition probabilities." << std::endl;
  typedef graph_type::vertex_id_type vertex_id_type;
  for(vertex_id_type vid = 0; vid < graph.num_vertices(); ++vid) {  
    double sum = 0;
    const graph_type::edge_list_type out_edges = graph.out_edges(vid);
    // Sum up weight on out edges
    for(size_t i = 0; i < out_edges.size(); ++i) 
      sum += graph.edge_data(out_edges[i]).weight;
    for(size_t i = 0; i < out_edges.size(); ++i) 
      graph.edge_data(out_edges[i]).weight /= sum;
  }
  logstream(LOG_INFO)
    << "Finished normalizing transition probabilities." << std::endl;
} // end of normalize_graph



#endif
