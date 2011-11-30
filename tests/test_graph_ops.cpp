

#include <iostream>

#include <graphlab.hpp>

typedef graphlab::graph<char, char> graph_type;

int main(int argc, char** argv) {
  const std::string in_fname = argv[1];
  const std::string format = argv[2];
  const std::string out_fname = argv[3];
  
  std::cout << "Loading graph: " << std::endl 
            << "\t In:     " << in_fname << std::endl
            << "\t Format: " << format << std::endl
            << "\t Out:    " << out_fname << std::endl;
  graph_type graph;
  bool success = 
    graphlab::graph_ops<graph_type>::load_structure(in_fname, format, graph);
  assert(success);
  std::cout << "Finished loading!" << std::endl;

  std::cout << "Finalizing" << std::endl;
  graph.finalize();
  std::cout << "Finished finalizing" << std::endl;

  std::cout << "Saving as metis" << std::endl;
  success = 
    graphlab::graph_ops<graph_type>::save_metis_structure(out_fname, graph);
  assert(success);
  std::cout << "Finished saving" << std::endl;


  return EXIT_SUCCESS;
}

