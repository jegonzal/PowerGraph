
#include <graphlab.hpp>


#include <stdio.h>
#include <unistd.h>


#include <graphlab/macros_def.hpp>



struct vertex_data {
  float dist;
  vertex_data() : dist(std::numeric_limits<float>::max()) { }
}; // end of vertex_data

struct edge_data { 
  float dist; 
  edge_data() : dist(1) { } 
}; // end of edge data

typedef graphlab::graph<vertex_data, edge_data> graph_type;


class shortest_path_update : public
graphlab::iupdate_functor<graph_type, shortest_path_update> {
public:
  void operator()(icontext_type& context) {
    vertex_data& vdata = context.vertex_data();    
    foreach(edge_id_type eid, context.in_edge_ids()) {
      const vertex_id_type nbr_id = context.source(eid);
      const vertex_data& nbr = context.const_vertex_data(nbr_id);
      const edge_data& edata = context.const_edge_data(eid);
      vdata.dist = std::min(vdata.dist, nbr.dist + edata.dist);
    }
    // Reschedule any affected neighbors
    foreach(edge_id_type eid, context.out_edge_ids()) {
      const vertex_id_type nbr_id = context.target(eid);
      const vertex_data& nbr = context.const_vertex_data(nbr_id);
      const edge_data& edata = context.const_edge_data(eid);
      if(nbr.dist > (vdata.dist + edata.dist)) 
        context.schedule(nbr_id, *this);    
    }
  }
}; // end of shortest path update functor


bool load_graph_from_file(const std::string& fname,
                          graph_type& graph);

int main(int argc, char** argv) {
  // Parse input
  graphlab::command_line_options clopts("Run Shortest Path Algorithm.");
  std::string graph_file;
  clopts.attach_option("graph", &graph_file, graph_file,
                       "The graph file.");
  clopts.add_positional("graph");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }

  // Create a graphlab core
  graphlab::core<graph_type, shortest_path_update> core;
  core.set_options(clopts);
  
  std::cout << "Loading graph from file" << std::endl;
  const bool success = load_graph_from_file(graph_file, core.graph());
  if(!success) {
    std::cout << "Error in reading file: " << graph_file
              << std::endl;
    return EXIT_FAILURE;
  }

  // Set the root
  const graph_type::vertex_id_type root = 0;
  core.graph().vertex_data(root).dist = 0;

  core.schedule(root, shortest_path_update());
  std::cout << "Running!" << std::endl;
  const double runtime = core.start();
  std::cout << "Runtime:      " << runtime << std::endl;
  std::cout << "Update Count: " << core.last_update_count() << std::endl;
  std::cout << "Update Freq:  " << (core.last_update_count() / runtime) 
            << std::endl;

}



// file format from BFS test
#define MAGIC_WORD  0x10102048
bool load_graph_from_file(const std::string& fname, 
                          graph_type& graph) {
  graphlab::binary_input_stream bin(fname.c_str());
  if(!bin.good()) return false;

  // write it 4B wise?
  uint32_t key = 0; bin.read(key);
  if (key != MAGIC_WORD) {
    std::cout << "Invalid file format!" << std::endl;
    return false;
  }

  bin.read(key); const bool need_back = key;
  assert(!need_back);
  
  uint32_t nverts = 0;
  bin.read(nverts);
  uint32_t nedges = 0;
  bin.read(nedges);

  std::cout << "Nverts:  " << nverts << std::endl;
  std::cout << "Nedges:  " << nedges << std::endl;  

  std::vector<uint32_t> offsets(nverts + 1, -1);
  bin.read_vector(offsets);

  graph.resize(nverts);
  for(size_t i = 0; i < nverts; ++i) {
    const size_t from = offsets[i];
    const size_t to = offsets[i+1]; 
    for(size_t j = from; j < to; ++j) {
      const graph_type::vertex_id_type from_node = i;
      uint32_t to_node = -1; bin.read(to_node);
      if ((from_node == 2187551) && (to_node == 2868359)) {
        std::cout << "invalid node id!" << std::endl;
      }
      graph.add_edge(from_node, to_node);
    }
  }
  graph.finalize();
  return true;
} // end of load graph from file





