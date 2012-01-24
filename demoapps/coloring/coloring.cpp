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

/**
 * This application demonstrates a simple graph coloring algorithm
 * that uses a simple greedy coloring heuristic with a first fit.
 */


#include <iostream>
#include <string>
#include <set>

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>

typedef size_t color_type;
const color_type UNCOLORED(-1);


/**
 * Vertex representation. Each vertex has a color.
 */
struct vertex_data {
  size_t color, saturation;
  vertex_data () : color(UNCOLORED), saturation(0) {};
};

/**
 * No edge data
 */
struct edge_data { };

typedef graphlab::graph<vertex_data, edge_data> graph_type;



/**
 * Update function for each vertex.
 */
struct coloring_update :
  public graphlab::iupdate_functor<graph_type, coloring_update> {
  /**
   * Collect all neighbors' colors and determine a color for
   * this vertex.
   */
  void operator()(icontext_type& context){
    vertex_data& vdata = context.vertex_data();
    std::set<color_type> neighbor_colors;
    // collect neighbor colors
    foreach (const edge_type& edge, context.in_edges()) {
      const color_type& color = context.const_vertex_data(edge.source()).color;
      if (UNCOLORED == color) continue;
      neighbor_colors.insert(color);
    }
    foreach (const edge_type& edge, context.out_edges()) {
      const color_type& color = context.const_vertex_data(edge.target()).color;
      if (UNCOLORED == color) continue;
      neighbor_colors.insert(color);
    }
    // find a unique color 
    for (vdata.color = 0; neighbor_colors.count(vdata.color); vdata.color++);
    vdata.saturation = neighbor_colors.size();
  }
};// end of update functor




/**
 * This aggregator tracks the color distribution and validates the
 * coloring.
 */       
class aggregator :
  public graphlab::iaggregator<graph_type, coloring_update, aggregator> {
private:
  color_type max_color, total_saturation;
  std::map<color_type, size_t> color_distribution;
public:
  aggregator() : max_color(0), total_saturation(0) { }
  void operator()(icontext_type& context) {
    const color_type color = context.const_vertex_data().color;
    max_color = std::max(max_color, color);
    total_saturation += context.const_vertex_data().saturation;
    color_distribution[color]++;
    // Check the coloring
    foreach (const edge_type& edge, context.in_edges()) 
      ASSERT_NE(context.const_vertex_data(edge.source()).color, color);
    foreach (const edge_type& edge, context.out_edges()) 
      ASSERT_NE(context.const_vertex_data(edge.target()).color, color);
  } // end of operator()
  void operator+=(const aggregator& other) { 
    max_color = std::max(max_color, other.max_color);
    total_saturation += other.total_saturation;
    typedef std::pair<color_type, size_t> pair_type;
    foreach(const pair_type& pair, other.color_distribution)
      color_distribution[pair.first] += pair.second;
  }
  void finalize(iglobal_context_type& context) {
    std::cout << "Total colors:       " << (max_color + 1) << std::endl 
              << "Average Saturation: " 
              << (total_saturation / context.num_vertices()) << std::endl;
    std::cout << "(color, #vertices): ";
    typedef std::pair<color_type, size_t> pair_type;
    foreach(const pair_type& pair, color_distribution)
      std::cout << '(' << pair.first << ", " << pair.second << ")  ";
    std::cout << std::endl;
  }
}; // end of aggregator




int main (int argc, char *argv[]){
  // set logging levels
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // parse command line options
  graphlab::command_line_options
    opts("Naive Greedy Graph Coloring");
  std::string graph_file;
  std::string format = "snap";
  // add graph and format
  opts.attach_option("graph", &graph_file,
                     "The graph file (required).");
  opts.add_positional("graph");
  opts.attach_option("format", &format, format,
                     "Format of the graph file.");                     
  if(!opts.parse (argc, argv)) return EXIT_FAILURE;
  if(!opts.is_set("graph")) {
    std::cout << "Graph file is required." << std::endl;
    return EXIT_FAILURE;
  }  

  // print command line options
  opts.print();
  std::cout << "Coloring Options -------------------" << std::endl;
  std::cout << "Input file:\t" << graph_file << std::endl;
 
  
  // set up graphlab execution core
  graphlab::core<graph_type, coloring_update> core;
  core.set_options(opts);
  if(!graphlab::graph_ops<graph_type>::
     load_structure (graph_file, format, core.graph())){
    std::cout << "Error loading graph from " << graph_file << "." << std::endl;
    return EXIT_FAILURE;
  }
  
  // execute graph updates
  core.schedule_all(coloring_update());
  std::cout << "Running graph coloring..." << std::endl;
  const double runtime = core.start();
  std::cout << "Done. Took " << runtime << " seconds." << std::endl;
  std::cout << "Checking coloring." << std::endl;
  // Compute the max colors
  core.add_aggregator("aggregator", aggregator(), 0);
  core.aggregate_now("aggregator");

  return EXIT_SUCCESS;
}


