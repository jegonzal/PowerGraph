#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include <graphlab.hpp>

size_t num_data = 0;
size_t current_cluster = 0;
size_t num_clusters = 0;

//helper function to normalize a vector;
void normalize_eivec(std::vector<float>& vec) {
  float sum = 0.0;
  for (size_t i = 0; i < vec.size(); ++i) {
    sum += vec[i] * vec[i];
  }
  sum = sqrt(sum);
  for (size_t i = 0; i < vec.size(); ++i) {
    vec[i] /= sum;
  }
}

struct evec_vertex_data {
  std::vector<float> vec;
  evec_vertex_data() :
      vec(num_clusters, 0.0) {
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << vec;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> vec;
  }
};

struct evec_edge_data {
  float val;
  evec_edge_data() :
      val(0.0) {
  }
  explicit evec_edge_data(float in_val) :
      val(in_val) {
  }
  void save(graphlab::oarchive& oarc) const {
    oarc << val;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> val;
  }
};

typedef graphlab::distributed_graph<evec_vertex_data, evec_edge_data> evec_graph_type;

// Read a line from a file and add values to vertices
bool evec_line_parser(evec_graph_type& graph, const std::string& filename,
    const std::string& textline) {
  std::stringstream strm(textline);
  size_t vid = 0;
  strm >> vid;
  float val = 0.0;
  strm >> val;
  graph.add_edge(vid, num_data + current_cluster + 1, evec_edge_data(val));

  return true;
}

//gather values belonging to the same data
struct id_and_val {
  std::vector<size_t> ids;
  std::vector<float> vals;

  id_and_val() :
      ids(), vals() {
  }
  id_and_val(size_t id, float val) :
      ids(), vals() {
    ids.push_back(id);
    vals.push_back(val);
  }

  id_and_val& operator+=(const id_and_val& other) {
    for (size_t i = 0; i < other.ids.size(); ++i) {
      ids.push_back(other.ids[i]);
      vals.push_back(other.vals[i]);
    }
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    size_t num = ids.size();
    oarc << num;
    for (size_t a = 0; a < num; ++a)
      oarc << ids[a];
    for (size_t a = 0; a < num; ++a)
      oarc << vals[a];
  }
  void load(graphlab::iarchive& iarc) {
    ids.clear();
    vals.clear();
    size_t num = 0;
    iarc >> num;
    for (size_t a = 0; a < num; ++a) {
      size_t id = 0;
      iarc >> id;
      ids.push_back(id);
    }
    for (size_t a = 0; a < num; ++a) {
      float val = 0;
      iarc >> val;
      vals.push_back(val);
    }
  }
};

//gather values belonging to the same data
class aggregate_values: public graphlab::ivertex_program<evec_graph_type,
    id_and_val>, public graphlab::IS_POD_TYPE {
public:
  //gather on out edges
  edge_dir_type gather_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::OUT_EDGES;
  }
  id_and_val gather(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    return id_and_val(edge.target().id() - num_data - 1, edge.data().val);
  }

  //get values and make a vector
  void apply(icontext_type& context, vertex_type& vertex,
      const gather_type& total) {
    const std::vector<size_t>& ids = total.ids;
    const std::vector<float>& vals = total.vals;
    for (size_t i = 0; i < ids.size(); ++i) {
      vertex.data().vec[ids[i]] = vals[i];
    }
  }

  edge_dir_type scatter_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::NO_EDGES;
  }
  void scatter(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
  }
};

void normalize_eigen_vector(evec_graph_type::vertex_type& v) {
  normalize_eivec(v.data().vec);
}

class evec_graph_writer {
public:
  std::string save_vertex(evec_graph_type::vertex_type v) {
    if(v.id()>num_data) return"";

    std::stringstream strm;
    strm << v.id() << " ";
    const std::vector<float>& vec = v.data().vec;
    for (size_t i = 0; i < vec.size(); ++i) {
      if (i != 0)
        strm << " ";
      strm << vec[i];
    }
    strm << "\n";

    return strm.str();
  }

  std::string save_edge(evec_graph_type::edge_type e) {
    return "";
  }
};

//read and normalize eigen vectors
int main(int argc, char** argv) {
  std::string datafile;
  size_t rank = 0;

  //parse command line
  graphlab::command_line_options clopts(
      "Normalize eigen vectors for graph partitioning");
  clopts.attach_option("data", datafile, "Input file prefix");
  clopts.attach_option("rank", rank, "Rank of Lanczos method");
  clopts.attach_option("clusters", num_clusters, "Number of clusters");
  clopts.attach_option("data-num", num_data, "Number of data points");
  if (!clopts.parse(argc, argv))
    return EXIT_FAILURE;
  if (datafile == "") {
    std::cout << "--data is not optional\n";
    return EXIT_FAILURE;
  }
  if (rank == 0) {
    std::cout << "--rank is not optional\n";
    return EXIT_FAILURE;
  }
  if (num_clusters == 0) {
    std::cout << "--clusters is not optional\n";
    return EXIT_FAILURE;
  }
  if (num_data == 0) {
    std::cout << "--data-num is not optional\n";
    return EXIT_FAILURE;
  }

  //load and normalize eigen vectors
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  evec_graph_type graph(dc, clopts);
  //special vertices corresponding to clusters
  for (size_t i = 0; i < num_clusters; ++i) {
    graph.add_vertex(num_data + i + 1, evec_vertex_data());
  }
  for (size_t i = 0; i < num_clusters; ++i) {
    current_cluster = i;
    std::stringstream vec_filename;
    vec_filename << datafile;
    vec_filename << ".U.";
    vec_filename << i;
    vec_filename << "_";
    graph.load(vec_filename.str(), evec_line_parser);
  }
  graph.finalize();

  graphlab::omni_engine<aggregate_values> engine(dc, graph, "sync", clopts);
  engine.signal_all();
  engine.start();

  graph.transform_vertices(normalize_eigen_vector);

  graph.save(datafile + ".compressed", evec_graph_writer(), false, true, false,
      1);

  graphlab::mpi_tools::finalize();

  return EXIT_SUCCESS;
}

