#include "extensions.hpp"
using namespace graphlab::extension;

int main(int argc, char** argv) {
  extension_graph graph;
  if (argc < 3) {
    std::cout << argv[0] << " [input prefix] [output prefix]\n"; 
    return 0;
  }
  graph.load_structure(argv[1], "snap");
  pagerank(graph, "pr", 0.01);
  graph.save_vertices(argv[2], "pr");
}
