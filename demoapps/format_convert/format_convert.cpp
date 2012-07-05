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
#include <graphlab.hpp>

int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  std::string ingraph, informat;
  std::string outgraph, outformat;
  bool gzip = true;
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Graph Format Conversion.", true);

  clopts.attach_option("ingraph", ingraph,
                       "The input graph file. Required ");
  clopts.attach_option("informat", informat,
                       "The input graph file format");
  clopts.attach_option("outgraph", outgraph,
                       "The output graph file. Required ");
  clopts.attach_option("outformat", outformat,
                       "The output graph file format");
  clopts.attach_option("outgzip", gzip,
                       "If output is to be gzip compressed"); 

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  
  typedef graphlab::distributed_graph<graphlab::empty, graphlab::empty> graph_type;
  graph_type graph(dc, clopts);

  dc.cout() << "Loading graph in format: "<< ingraph << std::endl;
  graph.load_format(ingraph, informat);
  graph.finalize();

  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  graph.save_format(outgraph, outformat, gzip);

  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main



