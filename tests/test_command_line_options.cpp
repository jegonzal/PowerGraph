// // #include <iostream>
// // #include <vector>


// // #include <graphlab/util/command_line_options.hpp>

// #include <string>

// #include <graphlab.hpp>

// #include <graphlab/macros_def.hpp>

// int main(int argc, char** argv) {
//   std::string filename;
//   size_t dimensions = 20;
//   double bound = 1E-5;
//   bool use_x = false;
//   std::vector<size_t> nsamples(1,10000);

//   // Parse command line options
//   graphlab::command_line_options clopts("Welcome to a the HelloWorld");
//   clopts.attach_option("file", &filename, 
//                      "The input filename (required)");
//   clopts.add_positional("file");
//   clopts.attach_option("dim",
//                        &dimensions, dimensions,
//                        "the dimension of the grid");
//   clopts.attach_option("bound",
//                        &bound, bound,
//                        "The termination bound");
//   clopts.attach_option("usex",
//                        &use_x, use_x,
//                        "Use algorithm x");
//   clopts.attach_option("nsamples",
//                        &nsamples, nsamples,
//                        "A vector of the number of samples"); 
//   clopts.set_scheduler_type("fifo");
//   clopts.set_scope_type("edge");
//   if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
//   if(!clopts.is_set("file")) {
//     std::cout << "Input file not provided" << std::endl;
//     clopts.print_description();
//     return EXIT_FAILURE;
//   }
// }

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
// #include <graphlab/metrics/metrics.hpp>


// Constants for the algorithm. Better way would be to
// pass them in the shared_data to the update function, but for
// the simplicity of this example, we simply define them here.

double termination_bound = 1e-5;
double random_reset_prob = 0.15;   // PageRank random reset probability

/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  float weight;
  float old_source_value;
  edge_data(float weight = 1) :
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
 * Predecleration of the graph file loading function.  Defined at
 * bottom of file for clarity.
 *
 * Load a graph file specified in the format:
 *
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph_from_file(const std::string& filename,
                          pagerank_graph& graph);

/**
 * Makes a small to graph.
 */
void make_toy_graph(pagerank_graph& graph);


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
  sum = random_reset_prob/scope.num_vertices() + 
    (1-random_reset_prob)*sum;
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






int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");

  // Metrics  
  graphlab::metrics app_metrics("app::pagerank");

  // Setup the parser
  graphlab::command_line_options
    clopts("Run the PageRank algorithm.");

  // Add some command line options
  std::string graph_file;
  clopts.attach_option("graph",
                       &graph_file, graph_file,
                       "The graph file.  If none is provided "
                       "then a toy graph will be created");
  clopts.attach_option("bound",
                       &termination_bound, termination_bound,
                       "residual termination threshold");
  clopts.attach_option("resetprob",
                       &random_reset_prob, random_reset_prob,
                       "Random reset probability");

  // Parse the command line input
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }

  // Create a graphlab core
  gl_types::core core;

  // Set the engine options
  core.set_engine_options(clopts);
  
  // Create or load graph depending on if the file was set
  app_metrics.start_time("load");
  if(graph_file.empty()) {
    // Create a synthetic graph
    make_toy_graph(core.graph());
  } else {
    // load the graph from the file
    bool success = load_graph_from_file(graph_file, core.graph());
    if(!success) {
      std::cout << "Error in reading file: " << graph_file
                << std::endl;
      return EXIT_FAILURE;
    }
  }
  app_metrics.stop_time("load");

  //Normalize the vertices
  app_metrics.start_time("normalize");
  double sum = 0;
  for(size_t i = 0; i < core.graph().num_vertices(); ++i) {
    core.graph().vertex_data(i).value = 
      graphlab::random::rand01() + 1;
    sum += core.graph().vertex_data(i).value;
  }
  for(size_t i = 0; i < core.graph().num_vertices(); ++i) {
    core.graph().vertex_data(i).value = 
      core.graph().vertex_data(i).value / sum;
  }
  app_metrics.stop_time("normalize");
  
  // Schedule all vertices to run pagerank update on the
  // first round.
  core.add_task_to_all(pagerank_update, 100.0);
  
  // Run the engine
  double runtime = core.start();
  
  // We are done, now output results.
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  
  // First we need to compute a normalizer. This could be done with
  // the sync facility, but for simplicity, we do it by hand.
  double norm = 0.0;
  for(graphlab::vertex_id_t vid = 0; 
      vid < core.graph().num_vertices(); vid++) {
            norm += core.graph().vertex_data(vid).value;
  }
  
  // And output 5 first vertices pagerank after dividing their value
  // with the norm.
  for(graphlab::vertex_id_t vid = 0; 
      vid < 5 && vid < core.graph().num_vertices(); vid++) {
    std::cout << "Page " << vid << " pagerank = " <<
      core.graph().vertex_data(vid).value / norm << '\n';
  }
  app_metrics.report(core.get_reporter());
	  
  return EXIT_SUCCESS;
} // End of main



// Creates simple 5x5 graph
void make_toy_graph(pagerank_graph& graph) {
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

} // end of make_toy_graph



/**
 * Load a graph file specified in the format:
 *
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *   source_id <tab> target_id <tab> weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph_from_file(const std::string& filename,
                          pagerank_graph& graph) {
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good() && !fin.eof()) {
    size_t source = 0;
    size_t target = 0;
    float weight = -1;
    fin >> source;
    if(!fin.good()) break;
    //  fin.ignore(1); // skip comma
    fin >> target;
    assert(fin.good());
    //  fin.ignore(1); // skip comma
    fin >> weight;
    assert(fin.good());
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
  std::cout 
    << "Finished loading graph with: " << std::endl
    << "\t Vertices: " << graph.num_vertices() << std::endl
    << "\t Edges: " << graph.num_edges() << std::endl;

  std::cout << "Normalizing out edge weights." << std::endl;
  // This could be done in graphlab but the focus of this app is
  // demonstrating pagerank
  for(graphlab::vertex_id_t vid = 0; 
      vid < graph.num_vertices(); ++vid) {
    vertex_data& vdata = graph.vertex_data(vid);
    // Initialze with self out edge weight
    double sum = vdata.self_weight;
    const graphlab::edge_list& out_eids = 
      graph.out_edge_ids(vid);
    // Sum up weight on out edges
    for(size_t i = 0; i < out_eids.size(); ++i) {
      graphlab::edge_id_t out_eid = out_eids[i];
      sum += graph.edge_data(out_eid).weight;
      
    }
    if (sum == 0) {
        vdata.self_weight = 1.0;
        sum = 1.0; //Dangling page
    }
    assert(sum > 0);
    // divide everything by sum
    vdata.self_weight /= sum;
    for(size_t i = 0; i < out_eids.size(); ++i) {
      graphlab::edge_id_t out_eid = out_eids[i];
      graph.edge_data(out_eid).weight /= sum;
    } 
  }
  std::cout << "Finished normalizing edes." << std::endl;


  std::cout 
    << "Finalizing graph." << std::endl
    << "\t This is required for the locking protocol to function correctly"
    << std::endl;
  graph.finalize();
  std::cout << "Finished finalization!" << std::endl;
  return true;
} // end of load graph
