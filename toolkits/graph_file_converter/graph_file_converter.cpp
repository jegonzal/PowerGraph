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

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

typedef graphlab::graph<char, char> graph_type;

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Graph File Converter.", true);
  std::string input_graph_file;
  std::string input_format = "metis";
  std::string output_graph_file;
  std::string output_format = "metis";

  clopts.attach_option("input",
                       &input_graph_file, input_graph_file,
                       "Input Graph File");
  clopts.add_positional("input");
  clopts.attach_option("inputformat",
                       &input_format, input_format,
                       "The input graph file format: {metis, snap, tsv}");
  clopts.add_positional("inputformat");
  clopts.attach_option("output",
                       &output_graph_file, output_graph_file,
                       "Output Graph File");
  clopts.add_positional("output");
  clopts.attach_option("outputformat",
                       &output_format, output_format,
                       "The output graph file format: {metis, patoh, zoltan}");
  clopts.add_positional("outputformat");

  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  if (input_graph_file.length() == 0) {
    std::cout << "Input Graph not specified.\n";
    return EXIT_FAILURE;
  }
  if (output_graph_file.length() == 0) {
    std::cout << "Output Graph not specified.\n";
    return EXIT_FAILURE;
  }


  if (!(output_format == "metis" || 
        output_format == "patoh" || 
        output_format == "zoltan")) {
    std::cout << "Invalid output format \"" << output_format << "\"\n";
  }
  
  if (!(input_format == "metis" || 
        input_format == "snap"  || 
        input_format == "tsv")) {
    std::cout << "Invalid input format \"" << input_format << "\"\n";
  }

  graph_type graph;
  const bool success = graphlab::graph_ops<graph_type>::
                    load_structure(input_graph_file, input_format, graph);
  if(!success) {
    std::cout << "Error in reading file: " << input_graph_file << std::endl;
  }
  graph.finalize();

  if (output_format == "metis") {
    graphlab::graph_ops<graph_type>::save_metis_structure(output_graph_file, 
                                                          graph);
  }
  else if (output_format == "patoh") {
    graphlab::
        graph_ops<graph_type>::
            save_patoh_hypergraph_structure(output_graph_file, graph);
  }
  else if (output_format == "zoltan") {
    graphlab::
        graph_ops<graph_type>::
            save_zoltan_hypergraph_structure(output_graph_file, graph);
  }
  return EXIT_SUCCESS;
} // End of main
