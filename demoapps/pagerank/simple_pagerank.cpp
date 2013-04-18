#include <vector>
#include <string>
#include <fstream>

#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>

#include <dlfcn.h>

using namespace graphlab;

typedef distributed_graph<float, graphlab::empty> graph_type;
typedef gl3engine<graph_type> engine_type;

// // define a neighborhood map reduce operation
// #define PAGERANK_MAP_REDUCE 0

// // Global random reset probability
// float RESET_PROB = 0.15;

// float TOLERANCE = 1.0E-2;

//TODO: turn these into function pointers
float pagerank_map(const graph_type::vertex_type& v);
void pagerank_combine(float& v1, const float& v2);
void pagerank_program(engine_type::context_type& context,
		      graph_type::vertex_type& vertex,
		      const engine_type::message_type& unused);

void init_vertex(graph_type::vertex_type& vertex) { vertex.data() = 1; }

int main(int argc, char** argv) {
  

  // Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  global_logger().set_log_level(LOG_INFO);

  // Parse command line options -----------------------------------------------
  graphlab::command_line_options clopts("PageRank algorithm.");
  std::string graph_dir;
  std::string format = "adj";
  std::string exec_type = "synchronous";
  clopts.attach_option("graph", graph_dir,
                       "The graph file. Required ");
  clopts.add_positional("graph");
  clopts.attach_option("format", format,
                       "The graph file format");
  clopts.attach_option("engine", exec_type, 
                       "The engine type synchronous or asynchronous");
  clopts.attach_option("tol", TOLERANCE,
                       "The permissible change at convergence.");
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
  dc.cout() << "#vertices: " << graph.num_vertices()
            << " #edges:" << graph.num_edges() << std::endl;

  // Initialize the vertex data
  graph.transform_vertices(init_vertex);

    // Running The Engine -------------------------------------------------------
  engine_type engine(dc, graph, clopts);

  // register the map reduce operation before usage
  // Each task registration must have a distinct ID ranging fro 0 to 223
  engine.register_map_reduce(PAGERANK_MAP_REDUCE,
			     pagerank_map,
			     pagerank_combine);

  engine.set_vertex_program(pagerank_program);
  engine.signal_all();
  engine.wait();
  
  // Save the final graph -----------------------------------------------------
  // if (saveprefix != "") {
  //   graph.save(saveprefix, pagerank_writer(),
  //              false,    // do not gzip
  //              true,     // save vertices
  //              false);   // do not save edges
  // }

  // Tear-down communication layer and quit -----------------------------------
  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // End of main
