/*
 *  PageRank tutorial application.
 *  pagerank.cpp
 *
 *  For description of the PageRank algorithm, see Wikipedia article
 *  http://en.wikipedia.org/wiki/Pagerank
 */

#include <string>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>


// Constants for the algorithm. Better way would be to
// pass them in the shared_data to the update function, but for
// the simplicity of this example, we simply define them here.

#define termination_bound 1e-5
#define damping_factor 0.85   // PageRank damping factor

/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  float weight;
  float old_source_value;
  edge_data(float weight) :
    weight(weight), old_source_value(0) { } 
    
}; // End of edge data


/**
 * Stores the value and the self weight
 */
struct vertex_data {
  float value;
  float self_weight; // GraphLab does not support edges from vertex to itself, so
  // we save weight of vertex's self-edge in the vertex data
  vertex_data(float value = 1) : value(value), self_weight(0) { }
}; // End of vertex data



//! The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> pagerank_graph;


/**
 * The collection of graphlab types restricted to the graph type used
 * in this program.
 */
typedef graphlab::types<pagerank_graph> gl_types;




/**
 * The Page rank update function
 */
void pagerank_update(gl_types::iscope &scope,
                     gl_types::icallback &scheduler,
                     gl_types::ishared_data* shared_data) {
  
  
                     
  // Get the data associated with the vertex
  vertex_data& vdata = scope.vertex_data();
  
  // Sum the incoming weights; start by adding the 
  // contribution from a self-link.
  float sum = vdata.value * vdata.self_weight;
  
  foreach(graphlab::edge_id_t eid, scope.in_edge_ids()) {
    // Get the neighobr vertex value
    const vertex_data& neighbor_vdata =
      scope.const_neighbor_vertex_data(scope.source(eid));
    double neighbor_value = neighbor_vdata.value;
    
    // Get the edge data for the neighbor
    edge_data& edata = scope.edge_data(eid);
    // Compute the contribution of the neighbor
    double contribution = edata.weight * neighbor_value;
    
    // Add the contribution to the sum
    sum += contribution;
    
    // Remember this value as last read from the neighbor
    edata.old_source_value = neighbor_value;
  }

  // compute the jumpweight
  sum = (1-damping_factor)/scope.num_vertices() + damping_factor*sum;
  vdata.value = sum;
  
  // Schedule the neighbors as needed
  foreach(graphlab::edge_id_t eid, scope.out_edge_ids()) {
    edge_data& outedgedata = scope.edge_data(eid);
    
    // Compute edge-specific residual by comparing the new value of this
    // vertex to the previous value seen by the neighbor vertex.
    double residual =
      outedgedata.weight *
      std::fabs(outedgedata.old_source_value - vdata.value);
    // If the neighbor changed sufficiently add to scheduler.
    if(residual > termination_bound) {
      gl_types::update_task task(scope.target(eid), pagerank_update);
      scheduler.add_task(task, residual);
    }
  }
} // end of pagerank update function





// Creates simple 5x5 graph
void create_graph(pagerank_graph& graph) {
  // Create 5 vertices
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());
  graph.add_vertex(vertex_data());

	
  // Page 0 links to page 3 only, so weight is 1
  graph.add_edge(0, 3, edge_data(1));
	
  // Page 1 links to 0 and 2
  graph.add_edge(1, 0, edge_data(0.5));
  graph.add_edge(1, 2, edge_data(0.5));
	
  // ... and so on
  graph.add_edge(2, 0, edge_data(1.0/3));
  graph.add_edge(2, 1, edge_data(1.0/3));
  graph.add_edge(2, 3, edge_data(1.0/3));
  graph.add_edge(3, 0, edge_data(0.25));
  graph.add_edge(3, 1, edge_data(0.25));
  graph.add_edge(3, 2, edge_data(0.25));
  graph.add_edge(3, 4, edge_data(0.25));
  graph.add_edge(4, 0, edge_data(0.2));
  graph.add_edge(4, 1, edge_data(0.2));
  graph.add_edge(4, 2, edge_data(0.2));
  graph.add_edge(4, 3, edge_data(0.2));
  // and self edge which must be handled specially from 4 to 4
  graph.vertex_data(4).self_weight = 0.2;

}



int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");
   
  // Setup the parser
  graphlab::command_line_options
    clopts("Run the PageRank algorithm.");

  // Create a graphlab core
  gl_types::core core;

  // Set the engine options
  core.set_engine_options(clopts);
  
  // Create a synthetic graph
  create_graph(core.graph());

  // Schedule all vertices to run pagerank update on the
  // first round.
  core.add_task_to_all(pagerank_update, 100.0);
  
  // Run the engine
  double runtime = core.start();

  // We are done, now output results.
  std::cout << "Graphlab finished, runtime: " << runtime << " seconds." << std::endl;
  
  
  // First we need to compute a normalizer. This could be
  // done with the sync facility, but for simplicity, we do
  // it by hand.
  double norm = 0.0;
  for(graphlab::vertex_id_t vid=0; vid<core.graph().num_vertices(); vid++) {
    norm += core.graph().vertex_data(vid).value;
  }
  
  // And output 5 first vertices pagerank after dividing their value
  // with the norm.
  for(graphlab::vertex_id_t vid=0; vid<5 && vid<core.graph().num_vertices(); vid++) {
    std::cout << "Page " << vid << " pagerank = " <<
      core.graph().vertex_data(vid).value / norm << std::endl;
  }
  
	  
  return EXIT_SUCCESS;
} // End of main


  // Configuration information
// std::string filename;

//   clopts.attach_option("infile", &filename,
//                        "PageRank input file. In src, dest, weight format.");
// 
//   // Parse the command line input
//   if(!clopts.parse(argc, argv)) {
//     std::cout << "Error in parsing input." << std::endl;
//     return EXIT_FAILURE;
//   }
//   if(!clopts.is_set("infile")) {
//     std::cout << "Input file no provided!" << std::endl;
//     clopts.print_description();
//     return EXIT_FAILURE;
//   }

// Load the graph
//   if(!load_graph(filename, core.graph())) {
//     std::cout << "Error in parsing graph!" << std::endl;
//     return EXIT_FAILURE;
//   }



/**
 * Load a graph file specified in the format:
 *
 *   source_id, target_id, weight
 *   source_id, target_id, weight
 *   source_id, target_id, weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph(const std::string& filename,
                pagerank_graph& graph) {
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good()) {
    size_t source = 0;
    size_t target = 0;
    float weight = -1;
    fin >> source;
    fin.ignore(1); // skip comma
    fin >> target;
    fin.ignore(1); // skip comma
    fin >> weight;
    //    fin.ignore(1); // skip comma
    // Ensure that the number of vertices is correct
    if(source >= graph.num_vertices() ||
       target >= graph.num_vertices())
      graph.resize(std::max(source, target) + 1);
    if(source != target) {
      // Add the edge
      edge_data edata(weight);
      graph.add_edge(source, target, weight);
    } else {
      // add the self edge by updating the vertex weight
      graph.vertex_data(source).self_weight = weight;
    }       
  }
  graph.finalize();
  return true;
} // end of load graph
