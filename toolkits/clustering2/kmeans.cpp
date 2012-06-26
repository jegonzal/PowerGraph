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

#include <limits>
#include <vector>
#include <iostream>

#include <graphlab.hpp>


size_t NUM_CLUSTERS = 0;

struct cluster {
  cluster(): count(0), changed(false) { }
  std::vector<double> center;
  size_t count;
  bool changed;

  void save(graphlab::oarchive& oarc) const {
    oarc << center << count << changed;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> center >> count >> changed;
  }
};

std::vector<cluster> CLUSTERS;

// the current cluster to initialize
size_t KMEANS_INITIALIZATION;

struct vertex_data{
  std::vector<double> point;
  size_t best_cluster;
  double best_distance;
  bool changed;

  void save(graphlab::oarchive& oarc) const {
    oarc << point << best_cluster << best_distance << changed;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> point >> best_cluster >> best_distance >> changed;
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
std::vector<double>& plus_equal_vector(std::vector<double>& a, 
                                       const std::vector<double>& b) {
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

graphlab::atomic<graphlab::vertex_id_type> NEXT_VID;

// Read a line from a file and creates a vertex
bool vertex_loader(graph_type& graph, const std::string& fname, 
                   const std::string& line) {
  if (line.empty()) return true;
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
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

  if (!success) return false;
  vtx.best_cluster = (size_t)(-1);
  vtx.best_distance = std::numeric_limits<double>::infinity();
  vtx.changed = false;
  graph.add_vertex(NEXT_VID.inc_ret_last(graph.numprocs()), vtx);
  return true;
}







// A set of Map Reduces to compute the maximum and minimum vector sizes
// to ensure that all vectors have the same length
struct max_point_size_reducer: public graphlab::IS_POD_TYPE {
  size_t max_point_size;

  static max_point_size_reducer get_max_point_size(const graph_type::vertex_type& v) {
    max_point_size_reducer r;
    r.max_point_size = v.data().point.size();
    return r;
  }

  max_point_size_reducer& operator+=(const max_point_size_reducer& other) {
    max_point_size = std::max(max_point_size, other.max_point_size);
    return *this;
  }
};

struct min_point_size_reducer: public graphlab::IS_POD_TYPE {
  size_t min_point_size;

  static min_point_size_reducer get_min_point_size(const graph_type::vertex_type& v) {
    min_point_size_reducer r;
    r.min_point_size = v.data().point.size();
    return r;
  }

  min_point_size_reducer& operator+=(const min_point_size_reducer& other) {
    min_point_size = std::min(min_point_size, other.min_point_size);
    return *this;
  }
};


/*
 * This transform vertices call is only used during
 * the initialization phase. IT computes distance to 
 * cluster[KMEANS_INITIALIZATION] and assigns itself
 * to the new cluster KMEANS_INITIALIZATION if the new distance
 * is smaller that its previous cluster asssignment
 */
void kmeans_pp_initialization(graph_type::vertex_type& v) {
  double d = sqr_distance(v.data().point, 
                          CLUSTERS[KMEANS_INITIALIZATION].center);
  if (v.data().best_distance > d) {
    v.data().best_distance = d;
    v.data().best_cluster = KMEANS_INITIALIZATION;
  }
}

/*
 * Draws a random sample from the data points that is 
 * proportionate to the "best distance" stored in the vertex.
 */
struct random_sample_reducer {
  std::vector<double> vtx;
  double weight;
 
  random_sample_reducer():weight(0) { }
  random_sample_reducer(const std::vector<double>& vtx, 
                        double weight):vtx(vtx),weight(weight) { }

  static random_sample_reducer get_weight(const graph_type::vertex_type& v) {
    if (v.data().best_cluster == (size_t)(-1)) {
      return random_sample_reducer(v.data().point, 1);
    }
    else {
      return random_sample_reducer(v.data().point,
                                   v.data().best_distance);
    }
  }

  random_sample_reducer& operator+=(const random_sample_reducer& other) {
    double totalweight = weight + other.weight;
    // if any weight is too small, just quit
    if (totalweight <= 0) return *this;

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

  void save(graphlab::oarchive &oarc) const {
    oarc << vtx << weight;
  }
  
  void load(graphlab::iarchive& iarc) {
    iarc >> vtx >> weight;
  }
};





/*
 * This transform vertices call is used during the 
 * actual k-means iteration. It computes distance to 
 * all "changed" clusters and reassigns itself if necessary
 */
void kmeans_iteration(graph_type::vertex_type& v) {
  // if current vertex's cluster was modified, we invalidate the distance.
  // and we need to recompute to all existing clusters
  // otherwise, we just need to recompute to changed cluster centers.
  size_t prev_asg = v.data().best_cluster;
  if (CLUSTERS[v.data().best_cluster].changed) {
    // invalidate. recompute to all
    v.data().best_cluster = (size_t)(-1);
    v.data().best_distance = std::numeric_limits<double>::infinity();
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      if (CLUSTERS[i].center.size() > 0) {
        double d = sqr_distance(v.data().point, 
                                CLUSTERS[i].center);
        if (d < v.data().best_distance) {
          v.data().best_distance = d;
          v.data().best_cluster = i;
        }
      }
    }
  }
  else {
    // just compute distance to what has changed
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      if (CLUSTERS[i].changed && CLUSTERS[i].center.size() > 0) {
        double d = sqr_distance(v.data().point, 
                                CLUSTERS[i].center);
        if (d < v.data().best_distance) {
          v.data().best_distance = d;
          v.data().best_cluster = i;
        }
      }
    }
  }
  v.data().changed = (prev_asg != v.data().best_cluster);
}







/*
 * computes new cluster centers
 * Also accumulates a counter counting the number of vertices which
 * assignments changed.
 */
struct cluster_center_reducer {
  std::vector<cluster> new_clusters;
  size_t num_changed;

  cluster_center_reducer():new_clusters(NUM_CLUSTERS), num_changed(0) { }

  static cluster_center_reducer get_center(const graph_type::vertex_type& v) {
    cluster_center_reducer cc;
    ASSERT_NE(v.data().best_cluster, (size_t)(-1));
    
    cc.new_clusters[v.data().best_cluster].center = v.data().point;
    cc.new_clusters[v.data().best_cluster].count = 1;
    cc.num_changed = v.data().changed;
    return cc;
  }

  cluster_center_reducer& operator+=(const cluster_center_reducer& other) {
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      if (new_clusters[i].count == 0) new_clusters[i] = other.new_clusters[i];
      else if (other.new_clusters[i].count > 0) {
        plus_equal_vector(new_clusters[i].center, other.new_clusters[i].center);
        new_clusters[i].count += other.new_clusters[i].count;
      }
    }
    num_changed += other.num_changed;
    return *this;
  }

  void save(graphlab::oarchive& oarc) const { 
    oarc << new_clusters << num_changed;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> new_clusters >> num_changed;
  }
};

struct vertex_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    for (size_t i = 0;i < v.data().point.size(); ++i) {
      strm << v.data().point[i] << "\t";
    }
    strm << v.data().best_cluster << "\n";
    strm.flush();
    return strm.str();
  }
 
  std::string save_edge(graph_type::edge_type e) { return ""; }
};


int main(int argc, char** argv) {
  std::cout << "Computes a K-means clustering of data.\n\n";

  graphlab::command_line_options clopts
    ("K-means clustering. The input data file is provided by the "
     "--data argument which is non-optional. The format of the data file is a "
     "collection of lines, where each line contains a comma or white-space " 
     "separated lost of numeric values representing a vector. Every line "
     "must have the same number of values. The required --clusters=N "
     "argument denotes the number of clusters to generate. To store the output "
     "see the --output-cluster and --output-data arguments");

  std::string datafile;
  std::string outcluster_file;
  std::string outdata_file;
  clopts.attach_option("data",
                       &datafile, datafile,
                       "Input file. Each line hold a white-space or comma separated numeric vector");
  clopts.attach_option("clusters",
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
  dc.cout() << "Number of datapoints: " << graph.num_vertices() << std::endl;

  if (graph.num_vertices() < NUM_CLUSTERS) {
    dc.cout() << "More clusters than datapoints! Cannot proceed" << std::endl;
    return EXIT_FAILURE;
  }

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


  dc.cout() << "Initializing using Kmeans++\n";
  // allocate clusters
  CLUSTERS.resize(NUM_CLUSTERS);
  for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
    CLUSTERS[i].center.resize(max_p_size);
  }

  // ok. perform kmeans++ initialization
  for (KMEANS_INITIALIZATION = 0; 
       KMEANS_INITIALIZATION < NUM_CLUSTERS;
       ++KMEANS_INITIALIZATION) {
    random_sample_reducer rs = graph.map_reduce_vertices<random_sample_reducer>
                                      (random_sample_reducer::get_weight);
    CLUSTERS[KMEANS_INITIALIZATION].center = rs.vtx;
    graph.transform_vertices(kmeans_pp_initialization);
  } 

  // "reset" all clusters
  for (size_t i = 0; i < NUM_CLUSTERS; ++i) CLUSTERS[i].changed = true;
  // perform Kmeans iteration 
  
  dc.cout() << "Running Kmeans...\n";
  bool clusters_changed = true;
  size_t iteration_count = 0;
  while(clusters_changed) {

    graph.transform_vertices(kmeans_iteration);
    cluster_center_reducer cc = graph.map_reduce_vertices<cluster_center_reducer>
                                    (cluster_center_reducer::get_center);  

    ++iteration_count;
    dc.cout() << "Kmeans iteration " << iteration_count << ": " <<
                 "# points with changed assignments = " << cc.num_changed << std::endl;

    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      double d = cc.new_clusters[i].count;
      if (d > 0) scale_vector(cc.new_clusters[i].center, 1.0 / d);
      if (cc.new_clusters[i].count == 0 && CLUSTERS[i].count > 0) {
        dc.cout() << "Cluster " << i << " lost" << std::endl;
        CLUSTERS[i].center.clear();
        CLUSTERS[i].count = 0;
        CLUSTERS[i].changed = false;
      }
      else {
        CLUSTERS[i] = cc.new_clusters[i];
        CLUSTERS[i].changed = true;
      }
    }
    clusters_changed = cc.num_changed > 0;

  }


  if (!outcluster_file.empty() && dc.procid() == 0) {
    dc.cout() << "Writing Cluster Centers..." << std::endl;
    std::ofstream fout(outcluster_file.c_str());
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      for (size_t j = 0; j < CLUSTERS[i].center.size(); ++j) {
        fout << CLUSTERS[i].center[j] << "\t";
      }
      fout << "\n";
    }
  }

  if (!outdata_file.empty()) {
    dc.cout() << "Writing Data with clister assignments...\n" << std::endl;
    graph.save(outdata_file, vertex_writer(), false, true, false, 1);
  }

  graphlab::mpi_tools::finalize();
}


