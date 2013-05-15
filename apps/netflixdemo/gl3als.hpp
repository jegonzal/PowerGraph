#ifndef GL3_ALS_HPP
#define GL3_ALS_HPP

#include <Eigen/Dense>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <graphlab/util/stl_util.hpp>
#include "stats.hpp"
// This file defines the serialization code for the eigen types.
#include "eigen_serialization.hpp"

using namespace graphlab;
const int SAFE_NEG_OFFSET = 2; //add 2 to negative node id
//to prevent -0 and -1 which arenot allowed

double TOLERANCE = 1e-3;
double LAMBDA = 0.01;
size_t MAX_ITER= -1;
double MAXVAL = 1e+100;
double MINVAL = -1e+100;
int    REGNORMAL = 0;
double TEST_PERCENT = 0.4;
bool INTERACTIVE = false;

/**
 * \brief We use the eigen library's vector type to represent
 * mathematical vectors.
 */
typedef Eigen::VectorXd vec_type;

/**
 * \brief We use the eigen library's matrix type to represent
 * matrices.
 */
typedef Eigen::MatrixXd mat_type;


/**
 * \ingroup toolkit_matrix_factorization
 *
 * \brief the vertex data type which contains the latent factor.
 *
 * Each row and each column in the matrix corresponds to a different
 * vertex in the ALS graph.  Associated with each vertex is a factor
 * (vector) of latent parameters that represent that vertex.  The goal
 * of the ALS algorithm is to find the values for these latent
 * parameters such that the non-zero entries in the matrix can be
 * predicted by taking the dot product of the row and column factors.
 */
struct vertex_data {
  /**
   * \brief A shared "constant" that specifies the number of latent
   * values to use.
   */
  static size_t NLATENT;
  /** \brief The number of times this vertex has been updated. */
  uint32_t nupdates;

  /** \brief The most recent L1 change in the factor value */
  float residual; //! how much the latent value has changed

  /** \brief The mean raiting of the user or movie */
  double bias;

  /** \brief The latent factor for this vertex */
  vec_type factor;

  std::vector<std::pair<double, graphlab::vertex_id_type> > top_rated;
  std::vector<std::pair<double, graphlab::vertex_id_type> > top_pred;

  /** \brief the user specific item-item similarity matrix */
  mat_type Wu;

  // std::vector<std::vector<std::pair<double, graphlab::vertex_id_type> > > top_explain;
  // std::vector<std::vector<std::pair<double, graphlab::vertex_id_type> > > top_explain2;
  // std::vector<std::vector<std::pair<double, graphlab::vertex_id_type> > > top_explain3;

  /**
   * \brief Simple default constructor which randomizes the vertex
   *  data
   */
  vertex_data() : nupdates(0), residual(1), bias(0) { randomize(); }
  /** \brief Randomizes the latent factor */
  void randomize() { factor.resize(NLATENT); factor.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(oarchive& arc) const {
    arc << nupdates << residual << bias << factor << top_rated << top_pred;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(iarchive& arc) {
    arc >> nupdates >> residual >> bias >> factor >> top_rated >> top_pred;
  }
}; // end of vertex data

size_t vertex_data::NLATENT = 20;

size_t hash_value(vertex_data const& v) {
  return v.nupdates;
}

/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data also stores the most recent error estimate.
 */
struct edge_data : public IS_POD_TYPE {
  /**
   * \brief The type of data on the edge;
   *
   * \li *Train:* the observed value is correct and used in training
   * \li *Validate:* the observed value is correct but not used in training
   * \li *Predict:* The observed value is not correct and should not be
   *        used in training.
   */
  enum data_role_type { TRAIN, VALIDATE, PREDICT  };

  /** \brief the observed value for the edge */
  float obs;

  /** \brief The train/validation/test designation of the edge */
  data_role_type role;

  /** \brief basic initialization */
  edge_data(float obs = 0, data_role_type role = TRAIN) :
    obs(obs), role(role) { }

}; // end of edge data


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */
typedef distributed_graph<vertex_data, edge_data> graph_type;

// Movie names
boost::unordered_map<graph_type::vertex_id_type, std::string> mlist;

stats_info count_edges(const graph_type::edge_type & edge){
  stats_info ret;

  if (edge.data().role == edge_data::TRAIN)
     ret.training_edges = 1;
  else if (edge.data().role == edge_data::VALIDATE)
     ret.validation_edges = 1;
  ret.max_user = (size_t)edge.source().id();
  ret.max_item = (-edge.target().id()-SAFE_NEG_OFFSET);
  return ret;
}



/////////////// Serialization /////////////////
#include "implicit.hpp"
/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph,
                         const std::string& filename,
                         const std::string& line) {
  ASSERT_FALSE(line.empty());
  if (boost::starts_with(line, "#") || boost::starts_with(line, "%")) {
    logstream(LOG_INFO) << line << std::endl;
    return true;
  }
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  namespace phoenix = boost::phoenix;
  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
  // Parse the line
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  const bool success = qi::phrase_parse
    (line.begin(), line.end(),
     //  Begin grammar
     (
      qi::ulong_[phoenix::ref(source_id) = qi::_1] >> -qi::char_(',') >>
      qi::ulong_[phoenix::ref(target_id) = qi::_1] >>
      -(-qi::char_(',') >> qi::float_[phoenix::ref(obs) = qi::_1])
      )
     ,
     //  End grammar
     ascii::space);

  if(!success) return false;

  if(role == edge_data::TRAIN || role == edge_data::VALIDATE){
    if (obs < MINVAL || obs > MAXVAL)
      logstream(LOG_FATAL)<<"Rating values should be between " << MINVAL << " and " << MAXVAL << ". Got value: " << obs << " [ user: " << source_id << " to item: " <<target_id << " ] " << std::endl;
  }

  // map target id into a separate number space
  target_id = -(vertex_id_type(target_id + SAFE_NEG_OFFSET));
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(obs, role));
  return true; // successful load
} // end of graph_loader

/**
 * \brief The prediction saver is used by the graph.save routine to
 * output the final predictions back to the filesystem.
 */
struct prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    return ""; //nop
  }
  std::string save_edge(const edge_type& edge) const {
    if(edge.data().role == edge_data::PREDICT) {
      std::stringstream strm;
      const double prediction =
        edge.source().data().factor.dot(edge.target().data().factor);
      strm << edge.source().id() << '\t';
      strm << (-edge.target().id() - SAFE_NEG_OFFSET) << '\t';
      strm << prediction << '\n';
      return strm.str();
    } else return "";
  }
}; // end of prediction_saver


struct linear_model_saver_U {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
  */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() > 0){
      std::string ret = boost::lexical_cast<std::string>(vertex.id()) + " ";
      for (uint i=0; i< vertex_data::NLATENT; i++)
        ret += boost::lexical_cast<std::string>(vertex.data().factor[i]) + " ";
        ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
};

struct linear_model_saver_V {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  /* save the linear model, using the format:
     nodeid) factor1 factor2 ... factorNLATENT \n
  */
  std::string save_vertex(const vertex_type& vertex) const {
    if (vertex.num_out_edges() == 0){
      std::string ret = boost::lexical_cast<std::string>(-vertex.id()-SAFE_NEG_OFFSET) + " ";
      for (uint i=0; i< vertex_data::NLATENT; i++)
        ret += boost::lexical_cast<std::string>(vertex.data().factor[i]) + " ";
        ret += "\n";
      return ret;
    }
    else return "";
  }
  std::string save_edge(const edge_type& edge) const {
    return "";
  }
};

struct error_aggregator : public graphlab::IS_POD_TYPE {
  double train;
  double test;
  size_t ntrain;
  size_t ntest;
  error_aggregator() : train(0), test(0), ntrain(0), ntest(0) {}
  error_aggregator& operator+=(const error_aggregator& other) {
    train += other.train;
    test += other.test;
    ntrain += other.ntrain;
    ntest += other.ntest;
    return *this;
  }
};
/**
 * \brief Given an edge compute the error associated with that edge
 */
error_aggregator extract_l2_error(const graph_type::edge_type & edge) {
  double pred =
    edge.source().data().factor.dot(edge.target().data().factor);
  pred = std::min(MAXVAL, pred);
  pred = std::max(MINVAL, pred);
  double err = (edge.data().obs - pred) * (edge.data().obs - pred);
  error_aggregator ret;
  if (edge.data().role == edge_data::TRAIN) {
    ret.train = err;
    ret.ntrain = 1;
  } else {
    ret.test = err;
    ret.ntest = 1;
  }
  return ret;
} // end of extract_l2_error

double extract_l2_norm(const graph_type::vertex_type& vertex) {
  return vertex.data().factor.norm();
}

template<typename key_t, typename val_t>
struct map_join {
  boost::unordered_map<key_t, val_t> data; // rating than movie
  void save(graphlab::oarchive& oarc) const {
    oarc << data;
  }
  void load(graphlab::iarchive& iarc) {
    iarc >> data;
  }
  map_join& operator+=(const map_join& other) {
    std::copy(other.data.begin(), other.data.end(),
              std::inserter(data, data.end()));
    return *this;
  }

  std::vector<std::pair<val_t, key_t> > get_top_k(size_t n) {
    std::vector<std::pair<val_t, key_t> > ret;
    typename boost::unordered_map<key_t, val_t>::const_iterator iter = data.begin();
    std::vector<std::pair<val_t, key_t> > all_copy;
    while (iter != data.end()) {
      all_copy.push_back(std::make_pair(iter->second, iter->first));
      ++iter;
    }
    std::sort(all_copy.rbegin(), all_copy.rend());
    size_t limit = all_copy.size() < n ? all_copy.size() : n;
    std::copy(all_copy.begin(), all_copy.begin() + limit,
              std::inserter(ret, ret.end()));
    return ret;
  }

  std::vector<std::pair<val_t, key_t> > get_top_k_exclude(
      size_t n, boost::unordered_map<key_t, val_t>& exclude) {
    std::vector<std::pair<val_t, key_t> > ret;
    typename boost::unordered_map<key_t, val_t>::const_iterator iter = data.begin();
    std::vector<std::pair<val_t, key_t> > all_copy;
    while (iter != data.end()) {
      if (exclude.count(iter->first) == 0) {
	      all_copy.push_back(std::make_pair(iter->second, iter->first));
      }
      ++iter;
    }
    std::sort(all_copy.rbegin(), all_copy.rend());
    size_t limit = all_copy.size() < n ? all_copy.size() : n;
    std::copy(all_copy.begin(), all_copy.begin() + limit,
              std::inserter(ret, ret.end()));
    return ret;
  }
};

bool is_user(const graph_type::vertex_type& vertex) {
  return vertex.num_out_edges() > 0;
}

bool is_movie(const graph_type::vertex_type& vertex) {
  return !is_user(vertex);
}

inline int id2movieid(graphlab::vertex_id_type vid) {
  return  -vid - SAFE_NEG_OFFSET;
}
#endif
