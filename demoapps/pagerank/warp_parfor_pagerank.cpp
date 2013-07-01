/*
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

#include <vector>
#include <string>
#include <fstream>

#include <graphlab.hpp>
#include <graphlab/engine/warp_graph_mapreduce.hpp>
#include <graphlab/engine/warp_parfor_all_vertices.hpp>
#include <graphlab/engine/warp_graph_transform.hpp>
using namespace graphlab;

// The graph type is determined by the vertex and edge data types
typedef distributed_graph<float , float> graph_type;

/*
 * A simple function used by graph.transform_vertices(init_vertex);
 * to initialize the vertes data.
 */
void init_vertex(graph_type::vertex_type& vertex) { vertex.data() = 1; }


float pagerank_map(graph_type::edge_type edge, graph_type::vertex_type other) {
  return other.data() / other.num_out_edges();
}

void pagerank(graph_type::vertex_type vertex) {
  vertex.data() = 0.15 + 0.85 * warp::map_reduce_neighborhood(vertex,
                                                              IN_EDGES,
                                                              pagerank_map);
}

/*
 * We want to save the final graph so we define a write which will be
 * used in graph.save("path/prefix", pagerank_writer()) to save the graph.
 */
struct pagerank_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data() << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of pagerank writer


int main(int argc, char** argv) {
  // Initialize control plain using mpi
  mpi_tools::init(argc, argv);
  distributed_control dc;

  // Parse command line options -----------------------------------------------
  command_line_options clopts("PageRank algorithm.");
  std::string graph_dir;
  std::string format = "tsv";
  clopts.attach_option("graph", graph_dir, "The graph file. ");
  clopts.attach_option("format", format, "The graph file format");
  size_t iterations = 10;
  clopts.attach_option("iterations", iterations,
                       "Number of asynchronous iterations to run");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "Prefix to save the output pagerank in");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  graph.load_format(graph_dir, format);
  // must call finalize before querying the graph
  graph.finalize();

  // Initialize the vertex data
  graph.transform_vertices(init_vertex);

  timer ti;
  for (size_t i = 0;i < iterations; ++i) {
    warp::parfor_all_vertices(graph, pagerank);
    std::cout << "Iteration " << i << " complete\n";
  }

  dc.cout() << "Finished Running in " << ti.current_time()
            << " seconds." << std::endl;


  // Save the final graph -----------------------------------------------------
  if (saveprefix != "") {
    graph.save(saveprefix, pagerank_writer(),
               false,    // do not gzip
               true,     // save vertices
               false);   // do not save edges
  }

  mpi_tools::finalize();
} // End of main




