
#include <cassert>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <graphlab/graph/graph.hpp>

#include <graphlab/util/binary_parser.hpp>

struct edge_data {
  float weight;
  float lastusedval;
  float srcvalue;
  float residual;
};

struct vertex_data {
  float value;
  float selfweight;
};


struct tripple{
  graphlab::vertex_id_t source;
  graphlab::vertex_id_t target;
  uint32_t count;
};


void make_graph(const std::vector<tripple>& data,
                graphlab::vertex_id_t max_vertex_id,
                const char* dest) {

  std::cout << "Generating graphlab graph." << std::endl;
  const size_t num_vertices = max_vertex_id + 1;
  const float unif_const = 1.0 / num_vertices;
  graphlab::blob_graph graph;
  std::cout << "\t (1/4) -> Adding vertices" << std::endl;
  for(size_t i = 0; i < num_vertices; ++i) {
    vertex_data vdata;
    vdata.value = unif_const;
    vdata.selfweight = 0;
    graph.add_vertex(graphlab::blob( sizeof(vertex_data), &vdata));
  }
  std::cout << "\t (2/4) -> Adding edges" << std::endl;
  for(size_t i = 0; i < data.size(); ++i) {
    const tripple& trip = data[i];
    if(trip.source != trip.target) {
      edge_data edata;
      edata.weight = trip.count;
      edata.lastusedval = unif_const;
      edata.srcvalue = unif_const;
      edata.residual = 0;
      graph.add_edge(trip.source, trip.target,
                     graphlab::blob( sizeof(edge_data), &edata));
    } else {
      vertex_data* vdata =
        graph.vertex_data(trip.source).as_ptr<vertex_data>();
      assert(vdata != NULL);
      vdata->selfweight = trip.count;
      
    }
  }
  std::cout << "\t (3/4) -> Normalizing the graph" << std::endl;
  for(graphlab::vertex_id_t v = 0; v < graph.num_vertices(); ++v) {
    vertex_data* vdata =
      graph.vertex_data(v).as_ptr<vertex_data>();
    assert(vdata != NULL);
    const std::vector<graphlab::edge_id_t>& out_edges =
      graph.out_edge_ids(v);
    float total = vdata->selfweight;
    for(size_t i = 0; i < out_edges.size(); ++i) {
      edge_data* edata =
        graph.edge_data(out_edges[i]).as_ptr<edge_data>();
      assert(edata != NULL);
      total += edata->weight;
    }
    assert(total > 0.0);
    // Divide through by total;
    vdata->selfweight /= total;
    for(size_t i = 0; i < out_edges.size(); ++i) {
      edge_data* edata =
        graph.edge_data(out_edges[i]).as_ptr<edge_data>();
      assert(edata != NULL);
      edata->weight /= total;
    }
  }
  std::cout << "\t (4/4) -> Finalizng the graph." << std::endl;
  graph.finalize();
  
  std::cout << "Saving an archive of the graph." << std::endl;
  graph.save(dest);
  std::cout << "Finished!" << std::endl;


}

void read_from_tsv(const char* src, const char* dest) {
  std::cout << "Reading file: " << src << std::endl;
  std::ifstream fin(src);
  assert(fin.good());
  std::vector<tripple> data;
  graphlab::vertex_id_t max_vertex_id(0);
  while(fin.good()) {
    tripple trip;
    fin >> trip.source >> trip.target >> trip.count;
    // Decrement to start at zero
    --trip.source;
    --trip.target;
    max_vertex_id = std::max(max_vertex_id,
                             std::max(trip.source, trip.target));
    data.push_back(trip);
  }
  fin.close();
  std::cout << "Read " << data.size() << " edges and "
            << max_vertex_id + 1 << " vertices.";

  make_graph(data, max_vertex_id, dest);
}

void read_from_binary(const char* src, const char* dest) {
  std::cout << "Reading file: " << src << std::endl;
  graphlab::binary_input_stream bin(src);
  assert(bin.good());
  uint32_t num_vertices = bin.read<uint32_t>();  
  uint32_t num_edges = bin.read<uint32_t>();
  std::cout << "Going to read " << num_vertices << " vertices and "
            << num_edges << " edges." << std::endl;
  std::vector<tripple> data(num_edges);
  graphlab::vertex_id_t max_vertex_id(0);
  while(bin.good()) {
    tripple trip;
    trip.source = bin.read<uint32_t>();
    trip.target = bin.read<uint32_t>();
    trip.count = bin.read<uint32_t>();
    // Decrement to start at zero
    --trip.source;
    --trip.target;
    max_vertex_id = std::max(max_vertex_id,
                             std::max(trip.source, trip.target));
    data.push_back(trip);
  }
  bin.close();
  std::cout << "Read " << data.size() << " edges and "
            << max_vertex_id + 1 << " vertices.";

  make_graph(data, max_vertex_id, dest);
}


int main(int argc, char** argv) {
  if(argc != 3) {
    std::cout << "Usage: " << std::endl
              << argv[0] << " input_graph.tsv output_graph.bin"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string srcfile(argv[1]);
  if(srcfile.substr(srcfile.length() -3, 3) == "tsv") {
    read_from_tsv(argv[1], argv[2]);
  } else {
    read_from_binary(argv[1], argv[2]);
  }
 
  return EXIT_SUCCESS;
}
