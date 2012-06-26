/**
 * This implements the classical "k-means" clustering algorithm.
 *
 * It takes as input file a series of lines where each line is a comma separated
 * or space separated list of values representing a vector. For instance:
 *
 * \verbatim
 * 1.1, 1.5, 0.9
 * 0.3, 0.4, -1.1
 * ...
 * \endverbatim
 *
 * It constructs a graph with a single vertex for each data point and simply
 * uses the "Map-Reduce" scheme to perform a k-means clustering of all
 * the datapoints.
 */


#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <vector>
#include <iostream>

#include <graphlab.hpp>


size_t NUM_CLUSTERS = 0;

struct cluster {
  std::vector<double> center;
  bool changed;
};
std::vector<cluster> CLUSTERS;

// the current cluster to initialize
size_t KMEANS_INITIALIZATION;

struct vertex_data{
  std::vector<double> point;
  size_t cluster;
  void save(graphlab::oarchive& oarc) const {
    oarc << point << cluster;
  }
  void load(graphlab::iarchive& iarc) {
    oarc >> point >> cluster;
  }
};




// helper function to compute distance between points
double sqr_distance(const std::vector<double>& a, 
                    const std::vector<double>& b) {
  ASSERT_EQ(a.size(), b.size());
  double total = 0;
  for (size_t i = 0;i < a.size(); ++i) {
    double d = a[i] - b[i];
    total += d * d;
  }
  return total;
}


// helper function to add two vectors 
std::vector<double>& operator+=(std::vector<double>& a, 
                                std::vector<double>& b) {
  ASSERT_EQ(a.size(), b.size());
  for (size_t i = 0;i < a.size(); ++i) {
    a[i] += b[i];
  }
  return a;
}

// helper function to scale a vector vectors 
std::vector<double>& scale_vector(std::vector<double>& a, double d) {
  for (size_t i = 0;i < a.size(); ++i) {
    a[i] *= d;
  }
  return a;
}


typedef graphlab::distributed_graph<vertex_data, graphlab::empty> graph_type;


// Read a line from a file and creates a vertex
bool vertex_loader(graph_type& graph, const std::string& fname, 
                   const std::string& line) {
  if (line.empty()) return;
  namespace phoenix = boost::phoenix;
  vertex_data vtx;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),       
     //  Begin grammar
     (
      (qi::double_[phoenix::push_back(phoenix::ref(vtx.point), qi::_1)] % -qi::char_(",") )
      )
     ,
     //  End grammar
     ascii::space); 

  vtx.cluster = 0;
  graph.add_vertex(vtx);
}







// A set of Map Reduces to compute the maximum and minimum vector sizes
// to ensure that all vectors have the same length
struct max_point_size_reducer {
  size_t max_point_size;

  static max_point_size_reducer get_max_point_size(const graph_type::vertex_type& v) {
    max_point_size_reducer r;
    r.max_point_size = v.data().point.size();
  }

  max_point_size_reducer& operator+=(const max_point_size_reducer& other) {
    max_point_size = std::max(max_point_size, other.max_point_size);
    return *this;
  }
};

struct min_point_size_reducer {
  size_t min_point_size;

  static min_point_size_reducer get_min_point_size(const graph_type::vertex_type& v) {
    min_point_size_reducer r;
    r.min_point_size = v.data().point.size();
  }

  min_point_size_reducer& operator+=(const min_point_size_reducer& other) {
    min_point_size = std::min(min_point_size, other.min_point_size);
    return *this;
  }
};



struct random_sample_reducer {
  std::vector<double> vtx;
  double weight;
 
  random_sample_reducer():weight(0) { }
  random_sample_reducer(graph_type::vertex_id_type vtx, 
                        double weight):vtx(vtx),weight(weight) { }

  static random_sample_reducer get_weight(const graph_type::vertex_type& v) {
    if (KMEANS_INITIALIZATION == -1) {
      return random_sample_reducer(v.data().point, 1);
    }
    else {
      return random_sample_reducer(v.data().point,
                                   sqr_distance(v.data().point, 
                                        cluster[KMEANS_INITIALIZATION].center));
    }
  }

  random_sample_reducer& operator+=(const random_sample_reducer& other) {
    double myp = weight / (weight + other.weight);
    if (graphlab::random::bernoulli(myp)) {
      weight += other.weight;
      return *this;
    }
    else {
      vtx = other.vtx;
      weight += other.weight;
      return *this;
    }
  }
};













int main(int argc, char** argv) {
  std::cout << "Computes a K-means clustering of data.\n\n";

  graphlab::command_line_options clopts
    ("K-means clustering. The input data file is provided by the "
     "--data argument which is non-optional. The format of the data file is a "
     "collection of lines, where each line contains a comma or white-space " 
     "separated lost of numeric values representing a vector. Every line "
     "must have the same number of values. The non-optional --cluster=N "
     "argument denotes the number of clusters to generate. To store the output "
     "see the --output-cluster and --output-data arguments");

  std::string datafile;
  std::string outcluster_file;
  std::string outdata_file;
  clopts.attach_option("data",
                       &datafile, datafile,
                       "Input file. Each line hold a white-space or comma separated numeric vector");
  clopts.attach_option("cluster",
                       &NUM_CLUSTERS, NUM_CLUSTERS,
                       "The number of clusters to create.");
  clopts.attach_option("output-clusters",
                       &outcluster_file, outcluster_file,
                       "If set, will write a file containing cluster centers "
                       "to this filename. This must be on the local filesystem "
                       "and must be accessible to the root node.");
  clopts.attach_option("output-data",
                       &outdata_file, outdata_file,
                       "If set, will output a copy of the input data with an additional "
                       "last column denoting the assigned cluster centers. The output "
                       "will be written to a sequence of filenames where each file is "
                       "prefixed by this value. This may be on HDFS.");

  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (datafile == "") {
    std::cout << "--data is not optional\n";
    return EXIT_FAILURE;
  }
  if (NUM_CLUSTERS == 0) {
    std::cout << "--cluster is not optional\n";
    return EXIT_FAILURE;
  }

  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  // load graph
  graph_type graph(dc, clopts);
  graph.load(datafile, vertex_loader);
  graph.finalize();
  dc.cout() << "Number of datapoints: " << graph.num_vertices() << std::endl
  dc.cout() << "Validating data..."; 
  // make sure all have the same array length
 
  size_t max_p_size = graph.map_reduce_vertices<max_point_size_reducer>
                                (max_point_size_reducer::get_max_point_size).max_point_size;

  size_t min_p_size = graph.map_reduce_vertices<min_point_size_reducer>
                                (min_point_size_reducer::get_min_point_size).min_point_size;

  if (max_p_size != min_p_size) {
    dc.cout() << "Point sizes range from " << min_p_size << " to " << max_p_size 
              << "! K-means cannot proceed!" << std::endl;
    return EXIT_FAILURE;
  }

  // allocate clusters
  CLUSTERS.resize(NUM_CLUSTERS);
  for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
    CLUSTERS[i].center.resize(max_p_size);
  }

  // ok. perform kmeans++ initialization
  for (KMEANS_INITIALIZATION = -1; 
       KMEANS_INITIALIZATION < NUM_CLUSTERS;
       ++KMEANS_INITIALIZATION) {
    random_sample_reducer rs = graph.map_reduce_vertices<random_sample_reducer>
                                      (random_sample_reducer::get_weight);
    CLUSTERS[i].center = rs.vtx;
  } 

 // perform Kmeans iteration 
 
}


