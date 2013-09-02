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

// The vertex data is just the id value
typedef graphlab::vertex_id_type vertex_data_type;
typedef graphlab::vertex_id_type edge_data_type;

typedef double gather_type;

// The graph type is determined by the vertex and edge data types
typedef graphlab::distributed_graph<vertex_data_type, edge_data_type> graph_type;

/*
 * A simple function used by graph.transform_vertices(init_vertex);
 * to initialize the vertex data to the vertex id
 */
void init_vertex(graph_type::vertex_type& vertex) { vertex.data() = -1; }

struct min_combiner : public graphlab::IS_POD_TYPE {
  vertex_data_type value;

  min_combiner() {
    value = -1;
  }
    
  min_combiner& operator+=(const min_combiner& other) { 
    if (other.value < value) {
        value = other.value;
    }
    return *this; 
  }
};

class concomp :
  public graphlab::ivertex_program<graph_type, gather_type, min_combiner>,
  public graphlab::IS_POD_TYPE {

  // set changed to determine which edges to scatter on
  bool changed;

  // local copy of the message value
  vertex_data_type message_value; 
public:
  // Receive inbound message (minimum data of adjacent vertices)
  void init(icontext_type& context, const vertex_type& vertex, const message_type& message) {
    // message.value == 4294967295 on first run, so init message_value to vertex data.
    if (message.value == 4294967295) {
      message_value = vertex.id();
    } else {
      // else, set the local copy to the message parameter.
      message_value = message.value;
    }

  }

  edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }

  // Change the vertex data if any of its neighbors have a lower data value.
  void apply(icontext_type& context, vertex_type& vertex, const gather_type& total) {
    // mark if values differ to determine which edges to scatter on.
    if (message_value < vertex.data()) {
        changed = true;
        vertex.data() = message_value;
    } else {
        changed = false;
    }
  }

  edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
    // If the vertex data changed, scatter along all edges. Otherwise stop.
    if (changed) {
        return graphlab::ALL_EDGES;
    } else {
        return graphlab::NO_EDGES;
    }
  }

  // Scatter to scatter_edges edges with the new message value.
  void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    bool isEdgeSource = (vertex.id() == edge.source().id());
    bool hasSameData = isEdgeSource ? (vertex.data() == edge.target().data()) : (vertex.data() == edge.source().data()) ;
    if (!hasSameData) {
      min_combiner combiner;
      combiner.value = message_value;

      context.signal(isEdgeSource ? edge.target() : edge.source(), combiner);
    }
  }
};

/* We want to save the final graph so we define a write which will be
* used in graph.save("path/prefix", concomp_writer()) to save the graph.
*/
struct concomp_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t" << v.data() << "\n";
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
}; // end of concomp writer
          

int main(int argc, char** argv) {
  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);
  
  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("Connected Components algorithm.");
  std::string graph_dir;
  std::string format = "snap";
  std::string execution_type = "synchronous";
  clopts.attach_option("graph", graph_dir, "The graph file. Required ");
  clopts.add_positional("graph");
  clopts.attach_option("format", format, "The graph file format");
  clopts.attach_option("execution", execution_type, "Execution type (synchronous or asynchronous)");

  std::string saveprefix;
  clopts.attach_option("saveprefix", saveprefix,
                       "If set, will save the resultant pagerank to a "
                       "sequence of files with prefix saveprefix");

  if(!clopts.parse(argc, argv)) {
    dc.cout() << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }
  if (graph_dir == "") {
    dc.cout() << "Graph not specified. Cannot continue";
    return EXIT_FAILURE;
  }
 
  // Build the graph ----------------------------------------------------------
  graph_type graph(dc, clopts);
  dc.cout() << "Loading graph in format: "<< format << std::endl;
  graph.load_format(graph_dir, format);
  // must call finalize before querying the graph
  graph.finalize();

  graph.transform_vertices(init_vertex);

  dc.cout() << "#vertices: " << graph.num_vertices() << " #edges:" << graph.num_edges() << std::endl;

  graphlab::omni_engine<concomp> engine(dc, graph, execution_type, clopts);

  min_combiner initial_message;
  initial_message.value = -1;

  engine.signal_all(initial_message);

  engine.start();

  const float runtime = engine.elapsed_seconds();
  dc.cout() << "Finished Running engine in " << runtime << " seconds." << std::endl;

  if (saveprefix != "") {
    graph.save(saveprefix, concomp_writer(),
                false,    // do not gzip
                true,     // save vertices
                false);   // do not save edges
  }

  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}
