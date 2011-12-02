

#include <iostream>

#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>

typedef graphlab::graph<char, char> graph_type;
typedef graphlab::graph_ops<graph_type> graph_ops;

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
 
  std::string in_fname;
  std::string out_fname;
  size_t threshold = 0;
  graphlab::command_line_options opts("Generate Graph files.");
  opts.attach_option("in_fname", &in_fname, in_fname,
                     "Input file");
  opts.attach_option("out_fname", &out_fname, out_fname,
                     "Output file");
  opts.attach_option("threshold", &threshold, threshold,
                     "Max neighbors");
  if(!opts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  
  graph_type graph;
  bool success = graph_ops::load_structure(in_fname, graph);
  assert(success);
  std::cout << "Finished loading!" << std::endl;

  std::cout << "Finalizing" << std::endl;
  graph.finalize();
  std::cout << "Finished finalizing" << std::endl;

  if(threshold > 0) {
    std::cout << "Allocating secondary graph datastructures" << std::endl;
    graph_type graph2;
    std::vector<bool> drop(graph.num_vertices(), false);
    graph2.resize(graph.num_vertices());
    std::cout << "Counting degree" << std::endl;
    size_t drop_count = 0;
    for(size_t i = 0; i < graph.num_vertices(); ++i) { 
      drop[i] = graph_ops::num_neighbors(graph, i) >  threshold;
      drop_count += drop[i];
    }
    std::cout << "Dropping: " << drop_count << std::endl;
    std::cout << "Building second graph" << std::endl;
    for(size_t i = 0; i < graph.num_vertices(); ++i) {
      const graph_type::edge_list_type out_edges = graph.out_edges(i);
      foreach(graph_type::edge_type edge, out_edges)
        if(!drop[edge.target()]) graph2.add_edge(i, edge.target());
    }
    std::cout << "Finalizing" << std::endl;
    graph2.finalize();
    std::cout << "Finished finalizing" << std::endl;
    std::cout << "Saving as metis" << std::endl;
    success = graph_ops::
      save_metis_structure(out_fname, graph2);
    assert(success);
    std::cout << "Finished saving" << std::endl;
  } else {
    std::cout << "Saving as metis" << std::endl;
    success = graph_ops::
      save_metis_structure(out_fname, graph);
    assert(success);
    std::cout << "Finished saving" << std::endl;
  }

  return EXIT_SUCCESS;
}

