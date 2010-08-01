#include <string>
#include <boost/program_options.hpp>
#include <graph/graph.hpp>

using namespace graphlab;

int main(int argc,  char *argv[]) {
  std::string filename;
  namespace boost_po = boost::program_options;
  // Create a description for this program
  boost_po::options_description
  desc("PageRank input file. Can be a Jure text file, or a graphlab serialized binary file");
  // Set the program options
  desc.add_options()
  ("infile",  boost_po::value<std::string>(&(filename))->default_value(""),
   "Input filename ");

  boost_po::variables_map vm;
  boost_po::store(boost_po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);
  boost_po::notify(vm);

  if (filename == "") {
    std::cout << "Input Graph needed\n";
    return 0;
  }

  graph g1, g2;
  g1.load(filename);
  g1.optimize_blob_layout();
  g2.load(filename);
  for (size_t i = 0;i < g1.num_vertices(); ++i) {
    if (g1.vertex_blob(i).size != g2.vertex_blob(i).size) {
      std::cout << "Vertex " << i << " differs in size\n";
    }
    if (memcmp(g1.vertex_blob(i).data, g2.vertex_blob(i).data, g2.vertex_blob(i).size) != 0)  {
      std::cout << "Vertex " << i << " differs in content\n";
    }
  }

  for (size_t i = 0;i < g1.num_edges(); ++i) {
    if (g1.edge_blob(i).size != g2.edge_blob(i).size) {
      std::cout << "Edge " << i << " differs in size\n";
    }
    if (memcmp(g1.edge_blob(i).data, g2.edge_blob(i).data, g2.edge_blob(i).size) != 0)  {
      std::cout << "Edge " << i << " differs in content\n";
    }
  }
  std::cout << "No differences!\n";
}
