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

#include "dist_pagerank.hpp"
#include <graphlab/macros_def.hpp>



void normalize_graph(graph_type& graph) {
  logstream(LOG_INFO)
    << "Optimizing graph layout in memory." << std::endl;
  graph.finalize();
  logstream(LOG_INFO)
    << "Coloring the graph." << std::endl;
  const size_t ncolors = graphlab::graph_ops<graph_type>::color(graph);
  std::cout << "Ncolors: " << ncolors << std::endl;
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





int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts
    ("Make atoms file for distribute pagerank", true);
  std::string graph_file;
  std::string format = "snap";
  std::string part_file;
  std::string idx_fname = "pagerank";
  size_t nparts = 64;

  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.");
  clopts.add_positional("graph");
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv}");
  clopts.attach_option("part", &part_file, part_file,
                       "The partition file");
  clopts.attach_option("output", &idx_fname, idx_fname,
                       "The output filename");
  clopts.attach_option("nparts", &nparts, nparts,
                       "The number of parts.");
  if(!clopts.parse(argc, argv) || !clopts.is_set("graph") || 
     !clopts.is_set("format") ) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  clopts.print();

  graph_type graph;
  if(format == "bin") {
    std::cout << "Loading binary graph." << std::endl;
    graph.load(graph_file);
  } else {
    std::cout << "Loading graph from structure file." << std::endl;
    const bool success = graphlab::graph_ops<graph_type>::    
      load_structure(graph_file, format, graph);
    if(!success) {
      std::cout << "Error in reading file: " << graph_file << std::endl;
    }
    normalize_graph(graph);
  }

  
  // call init_graph to create the graph
  std::vector<graphlab::graph_partitioner::part_id_type> parts(graph.num_vertices());
  if(!part_file.empty()) {
    graphlab::graph_partitioner::random_partition(graph, nparts, parts);
  } else {
    std::ifstream fin(part_file.c_str());
    for(size_t i = 0; i < parts.size(); ++i) { 
      fin >> parts[i];
      nparts = std::max(nparts, size_t(parts[i]) + 1);
    }
    fin.close();
  }
  graphlab::disk_graph<vertex_data, edge_data> 
    dg(idx_fname, nparts,
       graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM);
  dg.create_from_graph(graph, parts);




  return EXIT_SUCCESS;
} // End of main
