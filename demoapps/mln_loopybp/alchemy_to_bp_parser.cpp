
#include <cassert>

#include <iostream>


#include "factor_graph.hpp"
#include "factor_bp.hpp"

using namespace std;

int main(int argc, char** argv) {
  if(argc != 2) {
    cout << "<alchemy graph file>" << endl;
    return EXIT_FAILURE;
  }

  cout << "Making bp factor graph from alchemy graph." << endl;
  gl_types::graph graph;
  { // scoping to save some memory
    cout << "\t (1/5) Loading file: " << argv[1] << ". " << endl;
    factor_graph model;
    model.load_alchemy(argv[1]);
    cout << "\t (2/5) Rendering BP graphlab graph: " << endl ;
    make_bp_graph(model, graph);
  }
  cout << "\t (3/5) Finalizing graph." << endl;
  graph.finalize();
  cout << "\t (4/5) Saving Adjacency file." << endl;
  graph.save_adjacency(std::string("adjacency.csv"));
  cout << "\t (5/5) Saving graph." << endl;
  graph.save("graphlab_bp.bin");
  return EXIT_SUCCESS;
}
