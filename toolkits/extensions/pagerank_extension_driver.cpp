#include "extensions.hpp"
#include <graphlab/util/timer.hpp>
using namespace graphlab::extension;

int main(int argc, char** argv) {
  extension_graph graph;
  if (argc < 3) {
    std::cout << argv[0] << " [input prefix] [output prefix]\n"; 
    return 0;
  }
  graph.load_structure(argv[1], "snap");
  graphlab::timer ti;
  pagerank(graph, "pr", 0.01);
  std::cout << "PageRank in " << ti.current_time() << "s\n";
  graph.save_vertices(argv[2], "pr");
}
