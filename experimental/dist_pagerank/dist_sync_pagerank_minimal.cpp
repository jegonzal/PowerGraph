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
 * The factorized page rank update function
 */
class pagerank :
  public graphlab::ivertex_program<float, char>,
  public graphlab::IS_POD_TYPE {
public:
  
  void init(icontext_type& context,
            vertex_type& vertex) { vertex.data() = 1.0; }

  float gather(icontext_type& context, const vertex_type& vertex,
         edge_type& edge) const {
    return (edge.source().data() / edge.source().num_out_edges()) * (1.0 - RESET_PROB);
  }
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& total) {
    vertex.data() = total + RESET_PROB;
    if (vertex.id() == 0) std::cout << "v0: " << total << " " << vertex.data() << std::endl;
    context.signal(vertex);
  }
  
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
}; // end of factorized_pagerank update functor

typedef graphlab::distributed_graph<float, char> graph_type;

float float_identity(graph_type::vertex_type f) { return f.data(); }


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
  engine.initialize();

  std::cout << dc.procid() << ": Scheduling all" << std::endl;
  engine.signal_all();
  engine.start();
  
  
  float sum_of_graph = graph.map_reduce_vertices<float>(float_identity);
  std::cout << "Sum of graph: " << sum_of_graph << std::endl;
  
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


