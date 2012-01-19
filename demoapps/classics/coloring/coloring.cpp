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
 * This application demonstrates a simple graph coloring algorithm that uses
 * a simple greedy coloring heuristic with a first fit.
 */
#include "coloring.hpp"
#include <graphlab/macros_def.hpp>

using namespace std;

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
      set<color_type> neighbor_colors;
      color_type color;
      
      // collect neighbor colors
      foreach (edge_type edge, context.in_edges()){
        color = context.const_vertex_data(edge.source()).color;
        if (UNCOLORED == color) continue;
        neighbor_colors.insert(color);
      }
      
      foreach (edge_type edge, context.out_edges()){
        color = context.const_vertex_data(edge.target()).color;
        if (UNCOLORED == color) continue;
        neighbor_colors.insert(color);
      }
      
      // find a unique color
      for (color=0; ; color++){
        if (!neighbor_colors.count(color)) break;
      }
      
      vdata.color = color;
      vdata.saturation = neighbor_colors.size();
      
    }
  
};// end of update functor

/**
 * Initializes program configuration and parses the command line arguments
  for a graph file and a file format.
 * @param opts        command line options struct
 * @param graph_file  reference to variable to store path to graph file
 * @param format      reference to variable to store file format;
 *                    value will be used as default value
 */
static int init_options  (int argc, char *argv[],
                          graphlab::command_line_options &opts,
                          string &graph_file,
                          string &format);
static int print_results (const graph_type &graph);

ostream& operator<< (ostream& out, const vertex_data& vdata) {
  return out << "C=" << vdata.color << ", S=" << vdata.saturation;
}

int main (int argc, char *argv[]){

  // set logging levels
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // parse command line options
  graphlab::command_line_options
    opts("Naive Greedy Graph Coloring");
  string graph_file = NONE;
  string format = DEFAULT_FORMAT;
  if (0 > init_options (argc, argv, opts, graph_file, format)){
    return EXIT_FAILURE;
  }
  
  // set up graphlab execution core
  graphlab::core<graph_type, coloring_update> core;
  core.set_options(opts);
  if (!graphlab::
        graph_ops<graph_type>::
        load_structure (graph_file, format, core.graph())){
    cout << "Error loading graph from " << graph_file << "." << endl;
    return EXIT_FAILURE;
  }
  
  // execute graph updates
  core.schedule_all(coloring_update());
  cout << "Running graph coloring..." << endl;
  const double runtime = core.start();
  cout << "Done. Took " << runtime << " seconds." << endl;
  
  // print results
  print_results (core.graph());
  
  return SUCCESS;

}

static int init_options  (int argc, char *argv[],
                          graphlab::command_line_options &opts,
                          string &graph_file,
                          string &format){
  
  // add graph and format
  opts.attach_option(OPT_GRAPH_FILE, &graph_file,
                     "The graph file (required).");
  opts.add_positional("graph");
  opts.attach_option(OPT_FORMAT, &format, format,
                     "Format of the graph file.");
                     
  // fix scope
  opts.set_scope_type("edge");

  if (!opts.parse (argc, argv)) return ERR_INPUT;
  if (!opts.is_set(OPT_GRAPH_FILE)){
    cout << "Graph file is required." << endl;
    return ERR_INPUT;
  }
  
  // print command line options
  opts.print();
  cout << "Coloring Options -------------------" << endl;
  cout << "Input file:\t" << graph_file << endl;

  return SUCCESS;

}

static int print_results (const graph_type &graph){
  for (graph_type::vertex_id_type vid = 0;
      vid < graph.num_vertices();
      vid++){
    const graph_type::vertex_data_type& vdata = graph.vertex_data(vid);
    cout << "[" << int(vid) << "]: "  << vdata << endl;
  }
  return SUCCESS;
}

#include <graphlab/macros_undef.hpp>