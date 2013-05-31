#ifndef FEATURE_MF_HPP
#define FEATURE_MF_HPP

#include <Eigen/Dense>
#include <graphlab/util/stl_util.hpp>
#include "eigen_serialization.hpp"

using namespace graphlab;
//add 2 to negative node id to prevent -0 and -1 which arenot allowed
const int SAFE_NEG_OFFSET = 2;
double TOLERANCE = 1e-3;
double LAMBDA = 0.01;
double LAMBDA2 = 0.01;
size_t MAX_ITER= 10;
double TEST_PERCENT = -1;
bool INTERACTIVE = false;
double STEP = 0.001;
size_t NEDGES = 0;
size_t TOPK = 5;

int MAX_VAL = 5, MIN_VAL = 1;
size_t NTRAIN, NTEST;

bool USE_BIAS = true;
bool USE_BIAS_LATENT = true;
bool TRUNCATE = true;
bool USE_FEATURE_LATENT = false;
bool USE_FEATURE_WEIGHTS = false;

boost::unordered_map<vertex_id_type, std::string> mv_names;
bool is_initialized = false;
std::string saveprefix="result";


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

  /** \brief The latent factor for this vertex */
  vec_type factor, xfactor;

  double bias, fcache;

  /** \brief The feature for the vertex */
  boost::unordered_map<size_t, double> features;

  /** \brief The number of times this vertex has been updated. */
  uint32_t nupdates;

  /** \brief The L1 change in the factor value */
  double residual;

  /**
   * \brief Simple default constructor which randomizes the vertex
   *  data
   */
  vertex_data() : bias(0), fcache(0), nupdates(0), residual(1) { randomize(); }

  /** \brief Randomizes the latent factor */
  void randomize() { 
    factor.resize(NLATENT); factor.setRandom(); 
    xfactor.resize(NLATENT); xfactor.setZero(); 
  }

  /** \brief Save the vertex data to a binary archive */
  void save(oarchive& arc) const {
    arc << bias << fcache << nupdates << residual << factor 
        << xfactor << features;
  }

  /** \brief Load the vertex data from a binary archive */
  void load(iarchive& arc) {
    arc >> bias >> fcache >> nupdates >> residual >> factor
        >> xfactor >> features;
  }
}; // end of vertex data

size_t vertex_data::NLATENT = 20;

size_t hash_value(vertex_data const& v) {
  return boost::hash<float>()(v.residual);
}

/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data also stores the most recent error estimate.
 */
struct edge_data {
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
  float obs, pred;

  boost::unordered_map<size_t, double> features;

  /** \brief The train/validation/test designation of the edge */
  data_role_type role;

  /** \brief basic initialization */
  edge_data(float obs = 0, data_role_type role = TRAIN) : 
      obs(obs), pred(0), role(role) { }


  /** \brief Save the vertex data to a binary archive */
  void save(oarchive& arc) const {
    arc << obs << pred << role << features;
  }

  /** \brief Load the vertex data from a binary archive */
  void load(iarchive& arc) {
    arc >> obs >> pred >> role >> features;
  }
}; // end of edge data

/**
 * \brief The meta data stores everything needs to be displayed for visualization. 
 */
struct meta_data  {
  // boost::unordered_map<graph_type::vertex_id_type, std::string> rest_names;
};

/**
 * \brief Mapping from feature id to the vector of its latent factors. 
 * This map is maintain fully locally, and synchronized across all nodes every iteration.
 */
struct feature_table_type {
  double w0;
  std::vector<std::string> names;
  std::vector<size_t> keys;
  boost::unordered_map<size_t, size_t> key_lookup;

  std::vector<double> weights;
  std::vector<vec_type> latent;

  boost::hash<std::string> hash;

  void add_feature(const std::string& name) {
    size_t key = hash(name);
    if (!key_lookup.count(key)) {
      key_lookup[key] = keys.size();
      names.push_back(name);
      keys.push_back(key);
      weights.push_back(0);
      vec_type v; v.resize(vertex_data::NLATENT); v.setZero();
      latent.push_back(v);
    }
  }

  size_t get_key(std::string& name) {
    size_t key = hash(name);
    if (key_lookup.count(key))
      return key;
    else
      return 0;
  }

  size_t size() {
    return names.size();
  }

  void print_name() {
    for (size_t i = 0; i < names.size(); ++i) {
      std::cout<< names[i];
      if (i < names.size()-1) {
        std::cout << "\t";
      }
    }
    std::cout << std::endl;
  }

  void save(std::string& fname) {
    std::ofstream out;
    out.open(fname.c_str());
    for (size_t i = 0; i < keys.size(); ++i) {
      out << names[i] << "," << weights[i]; 
      for (size_t j = 0; i < latent[i].size(); ++j) {
        out << "," << latent[i][j];
      }
      out << "\n";
    }
    out.close();
  }
};

feature_table_type feature_table;

/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */
typedef distributed_graph<vertex_data, edge_data> graph_type;

/**
 * \brief The engine type used by the ALS matrix factorization
 * algorithm.
 *
 * The ALS matrix factorization algorithm currently uses the
 * synchronous engine.  However we plan to add support for alternative
 * engines in the future.
 */
typedef gl3engine<graph_type> engine_type;


/// Annoying global variable hacks...
distributed_control* dc_ptr;
graph_type* graph_ptr;



////////////// Helper functions /////////////////////
bool is_user(const graph_type::vertex_type& vertex) {
  return vertex.num_out_edges() > 0;
}

inline size_t vid2movieid(graphlab::vertex_id_type vid) {
  return  (size_t)(-vid - SAFE_NEG_OFFSET);
}

inline graphlab::vertex_id_type movieid2vid(size_t movieid) {
  return (graphlab::vertex_id_type)(-(movieid + SAFE_NEG_OFFSET));
}
#endif
