/*  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


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
#include <boost/tokenizer.hpp>

#include <limits>
#include <vector>
#include <map>
#include <iostream>
#include <stdlib.h>

#include <graphlab.hpp>


size_t NUM_CLUSTERS = 0;
bool IS_SPARSE = false;

struct cluster {
  cluster(): count(0), changed(false) { }
  std::vector<double> center;
  std::map<size_t, double> center_sparse;
  size_t count;
  bool changed;

  void save(graphlab::oarchive& oarc) const {
    oarc << center << count << changed << center_sparse;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> center >> count >> changed >> center_sparse;
  }
};

std::vector<cluster> CLUSTERS;

// the current cluster to initialize
size_t KMEANS_INITIALIZATION;

struct vertex_data{
  std::vector<double> point;
  std::map<size_t, double> point_sparse;
  size_t best_cluster;
  double best_distance;
  bool changed;

  void save(graphlab::oarchive& oarc) const {
    oarc << point << best_cluster << best_distance << changed << point_sparse;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> point >> best_cluster >> best_distance >> changed >> point_sparse;
  }
};

//use edges when edge weight file is given
struct edge_data {
  double weight;

  edge_data() :
      weight(0.0) {
  }
  explicit edge_data(double w) :
      weight(w) {
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << weight;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> weight;
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

double sqr_distance(const std::map<size_t, double>& a,
                    const std::map<size_t, double>& b) {
  double total = 0.0;
  for(std::map<size_t, double>::const_iterator iter = a.begin();
      iter != a.end(); ++iter){
    size_t id = (*iter).first;
    double val = (*iter).second;
    if(b.find(id) != b.end()){
      double d = val - b.at(id);
      total += d*d;
    }else{
      total += val * val;
    }
  }
  for(std::map<size_t, double>::const_iterator iter = b.begin();
      iter != b.end(); ++iter){
    double val = (*iter).second;
    if(a.find((*iter).first) == a.end()){
      total += val * val;
    }
  }

  return total;

////   cosine distance is better for sparse datapoints?
//    double ip = 0.0;
//    double lenA = 0.0;
//    double lenB = 0.0;
//    for(std::map<size_t, double>::const_iterator iter = a.begin();
//        iter != a.end(); ++iter){
//      size_t id = (*iter).first;
//      double val = (*iter).second;
//      if(b.find(id) != b.end()){
//        ip += val * b.at(id);
//      }
//      lenA += val*val;
//    }
//
//    if(ip == 0.0 || lenA == 0.0)
//      return 1.0;
//
//    for(std::map<size_t, double>::const_iterator iter = b.begin();
//        iter != b.end(); ++iter){
//      double val = (*iter).second;
//      lenB += val * val;
//    }
//
//    if(lenB == 1.0)
//      return 1.0;
//
//    return 1.0 - ip/(sqrt(lenA)*sqrt(lenB));

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

// helper function to add two vectors
std::map<size_t, double>& plus_equal_vector(std::map<size_t, double>& a,
                                       const std::map<size_t, double>& b) {
  for(std::map<size_t, double>::const_iterator iter = b.begin();
    iter != b.end(); ++iter){
    size_t id = (*iter).first;
    double val = (*iter).second;
    if(a.find(id) != a.end()){
      a[id] += b.at(id);
    }else{
      a.insert(std::make_pair(id, val));
    }
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

// helper function to scale a vector vectors
std::map<size_t, double>& scale_vector(std::map<size_t, double>& a, double d) {
  for(std::map<size_t, double>::iterator iter = a.begin();
    iter != a.end(); ++iter){
  size_t id = (*iter).first;
  double val = (*iter).second;
  a[id] = val*d;
//    (*iter).second *= d;
  }
  return a;
}


typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;

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
  graph.add_vertex(NEXT_VID.inc_ret_last(1), vtx);
  return true;
}

// Read a line from a file and creates a vertex
bool vertex_loader_sparse(graph_type& graph, const std::string& fname,
                   const std::string& line) {
  if (line.empty()) return true;

  vertex_data vtx;
  boost::char_separator<char> sep(" ");
  boost::tokenizer< boost::char_separator<char> > tokens(line, sep);
  BOOST_FOREACH (const std::string& t, tokens) {
    std::string::size_type pos = t.find(":");
    if(pos > 0){
      size_t id = (size_t)std::atoi(t.substr(0, pos).c_str());
      double val = std::atof(t.substr(pos+1, t.length() - pos -1).c_str());
      vtx.point_sparse.insert(std::make_pair(id, val));
    }
  }
  vtx.best_cluster = (size_t)(-1);
  vtx.best_distance = std::numeric_limits<double>::infinity();
  vtx.changed = false;
  graph.add_vertex(NEXT_VID.inc_ret_last(1), vtx);
  return true;
}

// Read a line from a file and creates a vertex
bool vertex_loader_with_id(graph_type& graph, const std::string& fname,
                   const std::string& line) {
  if (line.empty()) return true;
  size_t id = 0;
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  vertex_data vtx;
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(id) = qi::_1] >> -qi::char_(",") >>
      (qi::double_[phoenix::push_back(phoenix::ref(vtx.point), qi::_1)] % -qi::char_(",") )
      )
     ,
     //  End grammar
     ascii::space);

  if (!success) return false;
  vtx.best_cluster = (size_t)(-1);
  vtx.best_distance = std::numeric_limits<double>::infinity();
  vtx.changed = false;
  graph.add_vertex(id, vtx);
  return true;
}

// Read a line from a file and creates a vertex
bool vertex_loader_with_id_sparse(graph_type& graph, const std::string& fname,
                   const std::string& line) {
  if (line.empty()) return true;

  vertex_data vtx;
  size_t id = 0;
  boost::char_separator<char> sep(" ");
  boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
  bool first = true;
  BOOST_FOREACH (const std::string& t, tokens) {
    if(first){
      id = (size_t)std::atoi(t.c_str());
      first = false;
    }else{
      std::string::size_type pos = t.find(":");
      if(pos > 0){
        size_t id = (size_t)std::atoi(t.substr(0, pos).c_str());
        double val = std::atof(t.substr(pos+1, t.length() - pos -1).c_str());
        vtx.point_sparse.insert(std::make_pair(id, val));
      }
    }
  }
  vtx.best_cluster = (size_t)(-1);
  vtx.best_distance = std::numeric_limits<double>::infinity();
  vtx.changed = false;
  graph.add_vertex(id, vtx);
  return true;
}



//call this when edge weight file is given.
//each line should be [source id] [target id] [weight].
//directions of edges are ignored.
bool edge_loader(graph_type& graph, const std::string& filename,
    const std::string& textline) {
  if (textline.empty())
    return true;
  std::stringstream strm(textline);
  size_t source_vid = 0;
  size_t target_vid = 0;
  double weight = 0.0;
  strm >> source_vid;
  strm.ignore(1);
  strm >> target_vid;
  strm.ignore(1);
  strm >> weight;
  if(source_vid != target_vid)
    graph.add_edge(source_vid, target_vid, edge_data(weight));
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
 * the initialization phase. It computes distance to
 * cluster[KMEANS_INITIALIZATION] and assigns itself
 * to the new cluster KMEANS_INITIALIZATION if the new distance
 * is smaller that its previous cluster assignment
 */
void kmeans_pp_initialization(graph_type::vertex_type& v) {
  double d = sqr_distance(v.data().point,
                          CLUSTERS[KMEANS_INITIALIZATION].center);
  if (v.data().best_distance > d) {
    v.data().best_distance = d;
    v.data().best_cluster = KMEANS_INITIALIZATION;
  }
}

void kmeans_pp_initialization_sparse(graph_type::vertex_type& v) {
  double d = sqr_distance(v.data().point_sparse,
                          CLUSTERS[KMEANS_INITIALIZATION].center_sparse);
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

struct random_sample_reducer_sparse{
  std::map<size_t, double> vtx;
  double weight;

  random_sample_reducer_sparse():weight(0) { }
  random_sample_reducer_sparse(const std::map<size_t, double>& vtx,
                        double weight):vtx(vtx),weight(weight) { }

  static random_sample_reducer_sparse get_weight(const graph_type::vertex_type& v) {
    if (v.data().best_cluster == (size_t)(-1)) {
      return random_sample_reducer_sparse(v.data().point_sparse, 1);
    }
    else {
      return random_sample_reducer_sparse(v.data().point_sparse,
                                   v.data().best_distance);
    }
  }

  random_sample_reducer_sparse& operator+=(const random_sample_reducer_sparse& other) {
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
      if (CLUSTERS[i].center.size() > 0 || CLUSTERS[i].center_sparse.size() > 0) {
        double d = 0.0;
        if(IS_SPARSE == true)
          d = sqr_distance(v.data().point_sparse, CLUSTERS[i].center_sparse);
        else
          d = sqr_distance(v.data().point, CLUSTERS[i].center);
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
      if (CLUSTERS[i].changed &&
          (CLUSTERS[i].center.size() > 0 || CLUSTERS[i].center_sparse.size() > 0)) {
        double d = 0.0;
        if(IS_SPARSE == true)
          d = sqr_distance(v.data().point_sparse, CLUSTERS[i].center_sparse);
        else
          d= sqr_distance(v.data().point, CLUSTERS[i].center);
        if (d < v.data().best_distance) {
          v.data().best_distance = d;
          v.data().best_cluster = i;
        }
      }
    }
  }
  v.data().changed = (prev_asg != v.data().best_cluster);
}

//gathered information
//used when edge weight file is given
struct neighbor_info {
  std::map<size_t, double> cw_map;

  neighbor_info() :
      cw_map() {
  }
  neighbor_info(size_t clst, double weight) :
      cw_map() {
    cw_map.insert(std::make_pair(clst, weight));
  }

  neighbor_info& operator+=(const neighbor_info& other) {
    for (std::map<size_t, double>::const_iterator iter = other.cw_map.begin();
        iter != other.cw_map.end(); iter++) {
      size_t clst = iter->first;
      if (cw_map.find(clst) == cw_map.end()) {
        cw_map.insert(std::make_pair(clst, iter->second));
      } else {
        cw_map[clst] += iter->second;
      }
    }
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << cw_map;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> cw_map;
  }
};

//used when edge weight file is given
class cluster_assignment: public graphlab::ivertex_program<graph_type,
    neighbor_info>, public graphlab::IS_POD_TYPE {
public:
  //gather on all the edges
  edge_dir_type gather_edges(icontext_type& context,
      const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }

  //for each edge gather the weights and the assigned clusters of the neighbors
  neighbor_info gather(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
    if (edge.source().id() == vertex.id()) { //out edge
      return neighbor_info(edge.target().data().best_cluster,
          edge.data().weight);
    } else { //in edge
      return neighbor_info(edge.source().data().best_cluster,
          edge.data().weight);
    }
  }

  //assign a cluster, considering the clusters of neighbors
  void apply(icontext_type& context, vertex_type& vertex,
      const gather_type& total) {
    size_t past_clst = vertex.data().best_cluster;
    vertex.data().best_cluster = (size_t) (-1);
    vertex.data().best_distance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < NUM_CLUSTERS; ++i) {
      if (CLUSTERS[i].center.size() > 0 || CLUSTERS[i].center_sparse.size() > 0) {
        double d = 0.0;
        if(IS_SPARSE == true)
          d = sqr_distance(vertex.data().point_sparse, CLUSTERS[i].center_sparse);
        else
          d = sqr_distance(vertex.data().point, CLUSTERS[i].center);
        //consider neighbors
        const std::map<size_t, double>& cw_map = total.cw_map;
        for (std::map<size_t, double>::const_iterator iter = cw_map.begin();
            iter != cw_map.end(); iter++) {
          size_t neighbor_cluster = iter->first;
          double total_wieght = iter->second;
          if (i == neighbor_cluster)
            d -= total_wieght;
        }
        if (d < vertex.data().best_distance) {
          vertex.data().best_distance = d;
          vertex.data().best_cluster = i;
        }
      }
    }
    vertex.data().changed = (past_clst != vertex.data().best_cluster);
  }

  //send signals to the neighbors when the cluster assignment has changed
  edge_dir_type scatter_edges(icontext_type& context,
      const vertex_type& vertex) const {
    if (vertex.data().changed)
      return graphlab::ALL_EDGES;
    else
      return graphlab::NO_EDGES;
  }

  void scatter(icontext_type& context, const vertex_type& vertex,
      edge_type& edge) const {
  }
};



/*
 * computes new cluster centers
 * Also accumulates a counter counting the number of vertices which
 * assignments changed.
 */
struct cluster_center_reducer {
  std::vector<cluster> new_clusters;
  size_t num_changed;
  double cost;

  cluster_center_reducer():new_clusters(NUM_CLUSTERS), num_changed(0), cost(0) { }

  static cluster_center_reducer get_center(const graph_type::vertex_type& v) {
    cluster_center_reducer cc;
    ASSERT_NE(v.data().best_cluster, (size_t)(-1));

    if(IS_SPARSE == true)
      cc.new_clusters[v.data().best_cluster].center_sparse = v.data().point_sparse;
    else
      cc.new_clusters[v.data().best_cluster].center = v.data().point;
    cc.new_clusters[v.data().best_cluster].count = 1;
    cc.num_changed = v.data().changed;
    cc.cost = v.data().best_distance;
    return cc;
  }

  cluster_center_reducer& operator+=(const cluster_center_reducer& other) {
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      if (new_clusters[i].count == 0) new_clusters[i] = other.new_clusters[i];
      else if (other.new_clusters[i].count > 0) {
        if(IS_SPARSE == true)
          plus_equal_vector(new_clusters[i].center_sparse, other.new_clusters[i].center_sparse);
        else
          plus_equal_vector(new_clusters[i].center, other.new_clusters[i].center);
        new_clusters[i].count += other.new_clusters[i].count;
      }
    }
    num_changed += other.num_changed;
    cost += other.cost;
    return *this;
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << new_clusters << num_changed <<cost;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> new_clusters >> num_changed >> cost;
  }
};

struct vertex_writer {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    for (size_t i = 0;i < v.data().point.size(); ++i) {
      strm << v.data().point[i] << "\t";
    }
    strm << v.data().best_distance << "\t";
    strm << v.data().best_cluster << "\n";
    strm.flush();
    return strm.str();
  }

  std::string save_edge(graph_type::edge_type e) { return ""; }
};

struct vertex_writer_sparse {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    for(std::map<size_t, double>::iterator iter = v.data().point_sparse.begin();
        iter != v.data().point_sparse.end();++iter){
      strm << (*iter).first << ":" << (*iter).second << " ";
    }
    strm << v.data().best_cluster << "\n";
    strm.flush();
    return strm.str();
  }

  std::string save_edge(graph_type::edge_type e) { return ""; }
};

struct vertex_writer_with_id {
  std::string save_vertex(graph_type::vertex_type v) {
    std::stringstream strm;
    strm << v.id() << "\t";
    strm << v.data().best_cluster+1 << "\n";
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
  std::string edgedata_file;
  size_t MAX_ITERATION = 0;
  bool use_id = false;
  clopts.attach_option("data", datafile,
                       "Input file. Each line holds a white-space or comma separated numeric vector");
  clopts.attach_option("clusters", NUM_CLUSTERS,
                       "The number of clusters to create.");
  clopts.attach_option("output-clusters", outcluster_file,
                       "If set, will write a file containing cluster centers "
                       "to this filename. This must be on the local filesystem "
                       "and must be accessible to the root node.");
  clopts.attach_option("output-data", outdata_file,
                       "If set, will output a copy of the input data with an additional "
                       "two columns. The first added column is the distance to assigned "
		       "center and the last is the assigned cluster centers. The output "
                       "will be written to a sequence of filenames where each file is "
                       "prefixed by this value. This may be on HDFS.");
  clopts.attach_option("sparse", IS_SPARSE,
                       "If set to true, will use a sparse vector representation."
                       "The file format is [feature id]:[value] [feature id]:[value] ..."
                       ", where [feature id] must be positive integer or zero.");
  clopts.attach_option("id", use_id,
                       "If set to true, will use ids for data points. The id of a data point "
                       "must be written at the head of each line of the input data. "
                       "The output data will consist of two columns: the first one "
                       "denotes the ids; the second one denotes the assigned clusters.");
  clopts.attach_option("pairwise-reward", edgedata_file,
                       "If set, will consider pairwise rewards when clustering. "
                       "Each line of the file beginning with the argument holds [id1] [id2] "
                       "[reward]. This mode must be used with --id option.");
  clopts.attach_option("max-iteration", MAX_ITERATION,
                       "The max number of iterations");

  if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
  if (datafile == "") {
    std::cout << "--data is not optional\n";
    return EXIT_FAILURE;
  }
  if (NUM_CLUSTERS == 0) {
    std::cout << "--clusters is not optional\n";
    return EXIT_FAILURE;
  }
  if(edgedata_file.size() > 0){
    if(use_id == false){
      std::cout << "--id is not optional when you use edge data\n";
      return EXIT_FAILURE;
    }
  }

  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;
  // load graph
  graph_type graph(dc, clopts);
  NEXT_VID = (((graphlab::vertex_id_type)1 << 31) / dc.numprocs()) * dc.procid();
  if(IS_SPARSE == true){
    if(use_id){
      graph.load(datafile, vertex_loader_with_id_sparse);
    }else{
      graph.load(datafile, vertex_loader_sparse);
    }
  }else{
    if(use_id){
      graph.load(datafile, vertex_loader_with_id);
    }else{
      graph.load(datafile, vertex_loader);
    }
  }
  if(edgedata_file.size() > 0){
    graph.load(edgedata_file, edge_loader);
  }
  graph.finalize();
  dc.cout() << "Number of datapoints: " << graph.num_vertices() << std::endl;

  if (graph.num_vertices() < NUM_CLUSTERS) {
    dc.cout() << "More clusters than datapoints! Cannot proceed" << std::endl;
    return EXIT_FAILURE;
  }

  dc.cout() << "Validating data...";


  CLUSTERS.resize(NUM_CLUSTERS);
  // make sure all have the same array length
  if(IS_SPARSE == false){
    size_t max_p_size = graph.map_reduce_vertices<max_point_size_reducer>
                                  (max_point_size_reducer::get_max_point_size).max_point_size;
    size_t min_p_size = graph.map_reduce_vertices<min_point_size_reducer>
                                  (min_point_size_reducer::get_min_point_size).min_point_size;
    if (max_p_size != min_p_size) {
      dc.cout() << "Data has dimensionality ranging from " << min_p_size << " to " << max_p_size
                << "! K-means cannot proceed!" << std::endl;
      return EXIT_FAILURE;
    }
    // allocate clusters
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      CLUSTERS[i].center.resize(max_p_size);
    }
  }

  dc.cout() << "Initializing using Kmeans++\n";
  // ok. perform kmeans++ initialization
  for (KMEANS_INITIALIZATION = 0;
       KMEANS_INITIALIZATION < NUM_CLUSTERS;
       ++KMEANS_INITIALIZATION) {

    if(IS_SPARSE == true){
      random_sample_reducer_sparse rs = graph.map_reduce_vertices<random_sample_reducer_sparse>
                                        (random_sample_reducer_sparse::get_weight);
      CLUSTERS[KMEANS_INITIALIZATION].center_sparse = rs.vtx;
      graph.transform_vertices(kmeans_pp_initialization_sparse);
    }else{
      random_sample_reducer rs = graph.map_reduce_vertices<random_sample_reducer>
                                        (random_sample_reducer::get_weight);
      CLUSTERS[KMEANS_INITIALIZATION].center = rs.vtx;
      graph.transform_vertices(kmeans_pp_initialization);
    }
  }

  // "reset" all clusters
  for (size_t i = 0; i < NUM_CLUSTERS; ++i) CLUSTERS[i].changed = true;
  // perform Kmeans iteration

  dc.cout() << "Running Kmeans...\n";
  bool clusters_changed = true;
  size_t iteration_count = 0;
  while(clusters_changed) {
		if(MAX_ITERATION > 0 && iteration_count >= MAX_ITERATION)
			break;

    cluster_center_reducer cc = graph.map_reduce_vertices<cluster_center_reducer>
                                    (cluster_center_reducer::get_center);
    // the first round (iteration_count == 0) is not so meaningful
    // since I am just recomputing the centers from the output of the KMeans++
    // initialization
    if (iteration_count > 0) {
      dc.cout() << "Kmeans iteration " << iteration_count << ": " <<
                 "# points with changed assignments = " << cc.num_changed << 
		 " total cost: " << cc.cost << std::endl;
    }
    for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
      double d = cc.new_clusters[i].count;
      if(IS_SPARSE){
        if (d > 0) scale_vector(cc.new_clusters[i].center_sparse, 1.0 / d);
        if (cc.new_clusters[i].count == 0 && CLUSTERS[i].count > 0) {
          dc.cout() << "Cluster " << i << " lost" << std::endl;
          CLUSTERS[i].center_sparse.clear();
          CLUSTERS[i].count = 0;
          CLUSTERS[i].changed = false;
        }
        else {
          CLUSTERS[i] = cc.new_clusters[i];
          CLUSTERS[i].changed = true;
        }
      }else{
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
    }
    clusters_changed = iteration_count == 0 || cc.num_changed > 0;

    if(edgedata_file.size() > 0){
      clopts.engine_args.set_option("factorized", true);
      graphlab::omni_engine<cluster_assignment> engine(dc, graph, "async", clopts);
      engine.signal_all();
      engine.start();
    }else{
      graph.transform_vertices(kmeans_iteration);
    }

    ++iteration_count;
  }


  if (!outcluster_file.empty() && dc.procid() == 0) {
    dc.cout() << "Writing Cluster Centers..." << std::endl;
    std::ofstream fout(outcluster_file.c_str());
    if(IS_SPARSE){
      for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
        if(use_id)
          fout << i+1 << "\t";
        for (std::map<size_t, double>::iterator iter = CLUSTERS[i].center_sparse.begin();
             iter != CLUSTERS[i].center_sparse.end();++iter) {
          fout << (*iter).first << ":" << (*iter).second << " ";
        }
        fout << "\n";
      }
    }else{
      for (size_t i = 0;i < NUM_CLUSTERS; ++i) {
        if(use_id)
          fout << i+1 << "\t";
        for (size_t j = 0; j < CLUSTERS[i].center.size(); ++j) {
          fout << CLUSTERS[i].center[j] << " ";
        }
        fout << "\n";
      }
    }
  }

  if (!outdata_file.empty()) {
    dc.cout() << "Writing Data with cluster assignments...\n" << std::endl;
    if(use_id){
      graph.save(outdata_file, vertex_writer_with_id(), false, true, false, 1);
    }else{
      if(IS_SPARSE == true)
        graph.save(outdata_file, vertex_writer_sparse(), false, true, false, 1);
      else
        graph.save(outdata_file, vertex_writer(), false, true, false, 1);
    }
  }

  graphlab::mpi_tools::finalize();
}


