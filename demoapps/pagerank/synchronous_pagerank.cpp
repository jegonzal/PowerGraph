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



#include <vector>
#include <string>
#include <fstream>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

//! Global random reset probability
float RESET_PROB = 0.15;

/**
 * The type of data associated with each vertex
 */
typedef float vertex_data_type;

/**
 * The type of data associated with each edge.  Since there is no data
 * we set the edge type to graphlab empty which lets graphlab know not
 * to allocate memory for edge data.
 */
typedef graphlab::empty edge_data_type;


/**
 * The graph type naturally depends on the vertex and edge data type.
 */
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;


/**
 * The factorized page rank update function
 * extends ivertex_program specifying the:
 *   1) graph_type: the type of graph 
 *   2) gather_type: float (returned by the gather function). Note
 *      that the gather type is not strictly needed here since it is
 *      assumed to be the same as the vertex_data_type unless
 *      otherwise specified
 *
 * In addition ivertex program also takes a message type which is
 * assumed to be empty. Since we do not need messages no message type
 * is provided.
 *
 * pagerank also extends graphlab::IS_POD_TYP (is plain old data type)
 * which tells graphlab that the pagerank program can be serialized
 * (converted to a byte stream) by directly reading its in memory
 * representation.  If a vertex program does not exted
 * graphlab::IS_POD_TYPE it must implement load and save functions
 * (\todo see ref).
 *
 */
class pagerank :
  public graphlab::ivertex_program<graph_type, 
                                   float /* gather_type */ >,
  public graphlab::IS_POD_TYPE {
public:  
  /** Initialize the vertex program and vertex data */
  void init(icontext_type& context, vertex_type& vertex) { 
    vertex.data() = 1.0; 
  }

  /** Gather the weighted rank of the adjacent page   */
  float gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
    return ((1.0 - RESET_PROB) / edge.source().num_out_edges()) * 
      edge.source().data();
  }
  
  /** Use the total rank of adjacent pages to update this page */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& total) {
    vertex.data() = total + RESET_PROB;
    // Schedule this vertex to run again in the future
    context.signal(vertex);
  }
  
  /** Skip the scatter phase */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
}; // end of factorized_pagerank update functor


/**
 * Simple function used at the end of pagerank to extract the rank of
 * each page.  See: graph.map_reduce_vertices(float_identity);
 */
float float_identity(graph_type::vertex_type f) { return f.data(); }

struct pagerank_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data() << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) {
    return "";
  }
};

int main(int argc, char** argv) {
  //global_logger().set_log_level(LOG_DEBUG);
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
 
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_dir; 
  std::string format = "adj";
  clopts.attach_option("graph", &graph_dir, graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("format", &format, format,
                       "The graph file format: {metis, snap, tsv, adj, bin}");
  size_t powerlaw = 0;
  clopts.attach_option("powerlaw", &powerlaw, powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  std::string saveprefix;
  clopts.attach_option("saveprefix", &saveprefix, saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << dc.procid() << ": Starting." << std::endl;
  graph_type graph(dc, clopts);
  if(powerlaw > 0) {
    graphlab::graph_ops::load_synthetic_powerlaw(graph, powerlaw);
  } else {
    graphlab::graph_ops::load(graph, graph_dir, format);
  }
  graph.finalize();
  std::cout << "#vertices: " << graph.num_vertices() << " #edges:" << graph.num_edges() << std::endl;
  
  std::cout << dc.procid() << ": Creating engine" << std::endl;
  graphlab::synchronous_engine<pagerank> engine(dc, graph, clopts);

  std::cout << dc.procid() << ": Scheduling all" << std::endl;
  engine.signal_all();
  engine.start();
  
  
  float sum_of_graph = graph.map_reduce_vertices<float>(float_identity);
  std::cout << "Sum of graph: " << sum_of_graph << std::endl;

  if (saveprefix != "") {
    graphlab::graph_ops::save(graph,
                              saveprefix,
                              pagerank_writer(),
                              false,    // do not gzip
                              true,     // save vertices
                              false);   // do not save edges
  }
  
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


