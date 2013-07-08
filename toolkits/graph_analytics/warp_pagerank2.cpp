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
#include <graphlab/engine/gl3engine.hpp>
// #include <graphlab/macros_def.hpp>
using namespace graphlab;

#define PAGERANK_MAP_REDUCE 0

// Global random reset probability
float RESET_PROB = 0.15;

float TOLERANCE = 1E-2;

// The vertex data is just the pagerank value (a float)
typedef float vertex_data_type;

// There is no edge data in the pagerank application
typedef empty edge_data_type;

// The graph type is determined by the vertex and edge data types
typedef distributed_graph<vertex_data_type, edge_data_type> graph_type;

typedef gl3engine<graph_type> engine_type;
/*
 * A simple function used by graph.transform_vertices(init_vertex);
 * to initialize the vertes data.
 */
void init_vertex(graph_type::vertex_type& vertex) { vertex.data() = 1; }



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


float pagerank_map(const graph_type::vertex_type& v) {
  return v.data() / v.num_out_edges();
}

void pagerank_combine(float& v1, const float& v2) {
  v1 += v2;
}

void update_function(engine_type::context_type& context,
                     graph_type::vertex_type& vertex) {
  vertex.data() = 0.15 + 0.85 *
      context.map_reduce<float>(PAGERANK_MAP_REDUCE, IN_EDGES);
}

float pagerank_sum(graph_type::vertex_type v) {
  return v.data();
}


int main(int argc, char** argv) {
  // Initialize control plain using mpi
  mpi_tools::init(argc, argv);
  distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  command_line_options clopts("PageRank algorithm.");
  std::string graph_dir;
  std::string format = "adj";
  clopts.set_scheduler_type("fifo");
  clopts.attach_option("graph", graph_dir,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.add_positional("graph");
  clopts.attach_option("tol", TOLERANCE,
                       "The permissible change at convergence.");
  clopts.attach_option("format", format,
                       "The graph file format");
  size_t powerlaw = 0;
  clopts.attach_option("powerlaw", powerlaw,
                       "Generate a synthetic powerlaw out-degree graph. ");
  size_t iterations = 10;
  clopts.attach_option("iterations", iterations,
                       "Number of asynchronous iterations to run");
  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  if(powerlaw > 0) { // make a synthetic graph
    dc.cout() << "Loading synthetic Powerlaw graph." << std::endl;
    graph.load_synthetic_powerlaw(powerlaw, false, 2.1, 100000000);
  }
  else if (graph_dir.length() > 0) { // Load the graph from a file
    dc.cout() << "Loading graph in format: "<< format << std::endl;
    graph.load_format(graph_dir, format);
  }
  else {
    dc.cout() << "graph or powerlaw option must be specified" << std::endl;
    clopts.print_description();
    return 0;
  }
  // must call finalize before querying the graph
  graph.finalize();
  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  // Initialize the vertex data
  graph.transform_vertices(init_vertex);

  // Running The Engine -------------------------------------------------------
  engine_type engine(dc, graph, clopts);
  engine.register_map_reduce(PAGERANK_MAP_REDUCE,
                             pagerank_map,
                             pagerank_combine);

  timer ti;
  for (size_t i = 0;i < iterations; ++i) {
    engine.parfor_all_local_vertices(update_function);
    std::cout << "Iteration " << i << " complete\n";
    engine.wait();
  }

  dc.cout() << "Finished Running engine in " << ti.current_time()
            << " seconds." << std::endl;
  dc.cout() << engine.num_updates()
            << " updates." << std::endl;


  // Save the final graph -----------------------------------------------------
  if (saveprefix != "") {
    graph.save(saveprefix, pagerank_writer(),
               false,    // do not gzip
               true,     // save vertices
               false);   // do not save edges
  }

  mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main


// We render this entire program in the documentation


