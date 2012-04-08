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

/**
 * Vertex representation. Each vertex has an integer.
 */
struct vertex_data {
  int number;
  vertex_data () : number(0) {};
  vertex_data (int n) : number(n) {};
};

/**
 * No edge data
 */
struct edge_data { };

typedef graphlab::graph<vertex_data, edge_data> graph_type;
bool load_structure(std::string, graph_type&);

class sum_updater :
  public graphlab::iupdate_functor<graph_type, sum_updater> {
};

/**
 * This aggregator finds the sum of all integer values
 */       
class sum_aggregator :
  public graphlab::iaggregator<graph_type, sum_updater, sum_aggregator> {
private:
  int sum;
public:
  sum_aggregator() : sum(0) { }
  void operator()(icontext_type& context) {
    sum += context.const_vertex_data().number;
  } // end of operator()
  void operator+=(const sum_aggregator& other) {
    sum += other.sum;
  }
  void finalize(iglobal_context_type& context) {
    std::cout << "Total:\t\t" << sum << std::endl;
  }
}; // end of aggregator

int main (int argc, char *argv[]){

  // set logging levels
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // parse command line options
  graphlab::command_line_options opts("Simple Summation");
  std::string graph_file;

  // add graph and format
  opts.attach_option("graph", &graph_file,
                     "The graph file (required).");
  opts.add_positional("graph");                  
  if(!opts.parse (argc, argv)) return EXIT_FAILURE;
  if(!opts.is_set("graph")) {
    std::cout << "Graph file is required." << std::endl;
    return EXIT_FAILURE;
  } 

  // print command line options
  opts.print();
  std::cout << "Sum Options -------------------" << std::endl;
  std::cout << "Input file:\t" << graph_file << std::endl;
  
  // set up graphlab execution core
  graphlab::core<graph_type, sum_updater> core;
  core.set_options(opts);
  if(!load_structure (graph_file, core.graph())){
    std::cout << "Error loading graph from " << graph_file << "." << std::endl;
    return EXIT_FAILURE;
  }
  
  // execute graph aggregation
  core.add_aggregator("aggregator", sum_aggregator(), 0);
  core.aggregate_now("aggregator");

  return EXIT_SUCCESS;

}

bool load_structure(const std::string fname, graph_type& graph){
  
  logstream(LOG_INFO) << "Loading values file" << std::endl;
  
  namespace bios = boost::iostreams;
  std::ifstream in_file(fname.c_str(), 
    std::ios_base::in | std::ios_base::binary);
  bios::filtering_stream<bios::input> fin;
  fin.push(in_file);
  fin.set_auto_close(true);
  assert(fin.good());
   
  if(!fin.good()) {
    logstream(LOG_WARNING) << "file open failed" << std::endl;
    return false;
  }

  logstream(LOG_INFO) << "file open successful" << std::endl;
  
  // Loop over the contents
  size_t i = 0;
  while(fin.good()) {
    // Load a vertex
    int value = 0;
    try { fin >> value; } catch ( ... ) { 
      logstream(LOG_WARNING) << "Error reading value." << std::endl;
      return false;
    }
    graph.add_vertex(i++, vertex_data(value));
  }
  
  return true;
}