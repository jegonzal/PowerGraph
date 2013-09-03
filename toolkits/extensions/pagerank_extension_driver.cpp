#include "extensions.hpp"
#include <graphlab/util/timer.hpp>
using namespace graphlab::extension;

int main(int argc, char** argv) {
  extension_graph graph;
  if (argc < 2) {
    std::cout << argv[0] << " [input prefix] optional:[output prefix]\n"; 
    return 0;
  }
  graph.load_structure(argv[1], "snap");
  pagerank(graph, "pr", 0.01);
  if (argc > 2) {
    graph.save_vertices(argv[2], "pr");
  }
}
