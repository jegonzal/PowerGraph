/**
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
 * \file
 *
 * Matrix factorization with the Alternative Least Squares (ALS)  - parallel coordinate descent algorithm.
 * See the papers:
 * H.-F. Yu, C.-J. Hsieh, S. Si, I. S. Dhillon, Scalable Coordinate Descent Approaches to Parallel Matrix Factorization for Recommender Systems. IEEE International Conference on Data Mining(ICDM), December 2012.
 * Steffen Rendle, Zeno Gantner, Christoph Freudenthaler, and Lars Schmidt-Thieme. 2011. Fast context-aware recommendations with factorization machines. In Proceedings of the 34th international ACM SIGIR conference on Research and development in Information Retrieval (SIGIR '11). ACM, New York, NY, USA, 635-644.
 * Written by Danny Bickson, CMU
 */

#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab.hpp>
#include <graphlab/engine/gl3engine.hpp>
#include <Eigen/Dense>
#include "eigen_serialization.hpp"
#include <graphlab/macros_def.hpp>

typedef Eigen::VectorXd vec_type;

#define ALS_COORD_MAP_REDUCE 0
//when using negative node id range, we are not allowed to use
//0 and 1 so we add 2.
const static int SAFE_NEG_OFFSET=2;
const static int regnormal = 0;
static bool debug;
int max_iter = 0;

bool isuser(uint node){
  return ((int)node) >= 0;
}

/**
 * \ingroup toolkit_matrix_factorization
 *
 * \brief the vertex data type which contains the latent pvec.
 *
 * Each row and each column in the matrix corresponds to a different
 * vertex in the SGD graph.  Associated with each vertex is a pvec
 * (vector) of latent parameters that represent that vertex.  The goal
 * of the SGD algorithm is to find the values for these latent
 * parameters such that the non-zero entries in the matrix can be
 * predicted by taking the dot product of the row and column pvecs.
 */
struct vertex_data {
  /**
   * \brief A shared "constant" that specifies the number of latent
   * values to use.
   */
  static size_t NLATENT;
  /** \brief The latent pvec for this vertex */
  vec_type pvec;
  int t; //index inside the latent feature vector

  int nupdates;
  /**
   * \brief Simple default constructor which randomizes the vertex
   *  data
   */
  vertex_data() : nupdates(0),t(0) { if (debug) pvec = vec_type::Ones(NLATENT); else randomize(); }
  /** \brief Randomizes the latent pvec */
  void randomize() { pvec.resize(NLATENT); pvec.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const {
    arc << nupdates << pvec << t;
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) {
    arc >> nupdates >> pvec >> t;
  }
}; // end of vertex data

std::size_t hash_value(vertex_data const& b) {
  return b.nupdates;
}


/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data sgdo stores the most recent error estimate.
 */
struct edge_data : public graphlab::IS_POD_TYPE {
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
  edge_data(float obs = 0, data_role_type role = PREDICT) :
    obs(obs), role(role) { }

}; // end of edge data

std::size_t hash_value(edge_data const& b) {
  return boost::hash_value(b.obs);
}


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;
typedef graphlab::gl3engine<graph_type> engine_type;

bool isuser_node(const graph_type::vertex_type& vertex){
  return isuser(vertex.id());
}



/**
 * \brief The graph loader function is a line parser used for
 * distributed graph construction.
 */
inline bool graph_loader(graph_type& graph,
                         const std::string& filename,
                         const std::string& line) {
  ASSERT_FALSE(line.empty());
  // Determine the role of the data
  edge_data::data_role_type role = edge_data::TRAIN;
  if(boost::ends_with(filename,".validate")) role = edge_data::VALIDATE;
  else if(boost::ends_with(filename, ".predict")) role = edge_data::PREDICT;
  // Parse the line
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  strm >> source_id >> target_id;

  // for test files (.predict) no need to read the actual rating value.
  if(role == edge_data::TRAIN || role == edge_data::VALIDATE) {
    strm >> obs;
  }
  target_id = -(graphlab::vertex_id_type(target_id + SAFE_NEG_OFFSET));

  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, edge_data(obs, role));
  return true; // successful load
} // end of graph_loader

double LAMBDA = 0.001;

class gather_type {
public:
  double numerator;
  double denominator;

  gather_type() { 
    numerator = 0;
    denominator = LAMBDA;
  }

  gather_type(double numerator, double denominator) : numerator(numerator),
     denominator(denominator){
  }

  /** \brief Save the values to a binary archive */
  void save(graphlab::oarchive& arc) const { arc << numerator << denominator; }

  /** \brief Read the values from a binary archive */
  void load(graphlab::iarchive& arc) { arc >> numerator >> denominator; }  

  /** 
   * sums up values
   */
  gather_type& operator+=(const gather_type& other) {
    numerator += other.numerator;
    denominator += other.denominator;
    return *this;
  } // end of operator+=

}; // end of gather type



gather_type als_coord_map(const graph_type::vertex_type& center,
                         graph_type::edge_type& edge,
                         const graph_type::vertex_type& other) {
   //compute numerator of equation (5) in ICDM paper above
   //             (A_ij        - w_i^T*h_j  + wit          * h_jt         
   gather_type ret(edge.data().obs - center.data().pvec.dot(other.data().pvec)
                               + center.data().pvec[center.data().t] * other.data().pvec[center.data().t],
   //compute denominator of equation (5) in ICDM paper above
   //h_jt^2
     pow(other.data().pvec[center.data().t], 2));
   return ret;

}

void als_coord_combine(gather_type& v1, const gather_type& v2) {
    v1 += v2;
}


void als_coord_function(engine_type::context_type& context,
                  graph_type::vertex_type& vertex) {
       
   double regularization = LAMBDA;
   for (vertex.data().t=0; vertex.data().t< (int)vertex_data::NLATENT; vertex.data().t++){
     gather_type frac =  context.map_reduce<gather_type>(ALS_COORD_MAP_REDUCE, graphlab::ALL_EDGES);
     assert(frac.denominator > 0);
     double z = (frac.numerator/frac.denominator+regularization);  
     vertex.data().pvec[vertex.data().t] = z;
   }
}


/**
 * \brief Given an edge compute the error associated with that edge
 */
double extract_l2_error(const graph_type::edge_type & edge) {
  double pred =
      edge.source().data().pvec.dot(edge.target().data().pvec);
  double rmse = (edge.data().obs - pred) * (edge.data().obs - pred);
  return rmse;
} // end of extract_l2_error



size_t vertex_data::NLATENT = 20;

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  // Parse command line options -----------------------------------------------
  const std::string description =
      "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string input_dir;
  size_t interval = 0;
  std::string exec_type = "synchronous";
  clopts.attach_option("matrix", input_dir,
                       "The directory containing the matrix file");
  clopts.add_positional("matrix");
  clopts.attach_option("D", vertex_data::NLATENT,
                       "Number of latent parameters to use.");
  clopts.attach_option("interval", interval,
                       "The time in seconds between error reports");
  clopts.attach_option("lambda", LAMBDA,
                       "SGD regularization weight");
  clopts.attach_option("max_iter", max_iter,
                       "number of iterations");
  if(!clopts.parse(argc, argv) || input_dir == "") {
    std::cout << "Error in parsing command line arguments." << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer;
  graph_type graph(dc, clopts);
  graph.load(input_dir, graph_loader);
  dc.cout() << "Loading graph. Finished in "
            << timer.current_time() << std::endl;

  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in "
            << timer.current_time() << std::endl;


  dc.cout()
      << "========== Graph statistics on proc " << dc.procid()
      << " ==============="
      << "\n Num vertices: " << graph.num_vertices()
      << "\n Num edges: " << graph.num_edges()
      << "\n Num replica: " << graph.num_replicas()
      << "\n Replica to vertex ratio: "
      << float(graph.num_replicas())/graph.num_vertices()
      << "\n --------------------------------------------"
      << "\n Num local own vertices: " << graph.num_local_own_vertices()
      << "\n Num local vertices: " << graph.num_local_vertices()
      << "\n Replica to own ratio: "
      << (float)graph.num_local_vertices()/graph.num_local_own_vertices()
      << "\n Num local edges: " << graph.num_local_edges()
      //<< "\n Begin edge id: " << graph.global_eid(0)
        << "\n Edge balance ratio: "
        << float(graph.num_local_edges())/graph.num_edges()
        << std::endl;

  dc.cout() << "Creating engine" << std::endl;

  engine_type engine(dc, graph, clopts);
  engine.register_map_reduce(ALS_COORD_MAP_REDUCE, als_coord_map, als_coord_combine);
  for (int i=0; i< max_iter; i++){
     engine.parfor_all_local_vertices(als_coord_function);
     engine.wait();
     double rmse = graph.map_reduce_edges<double>(extract_l2_error);
     dc.cout() << "RMSE = " << sqrt(rmse / graph.num_edges()) << std::endl;
  }

  const double runtime = timer.current_time();
  dc.cout() << "----------------------------------------------------------"
            << std::endl
            << "Final Runtime (seconds):   " << runtime;

  // Compute the final training error -----------------------------------------
  dc.cout() << "Final error: " << std::endl;


  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
} // end of main



