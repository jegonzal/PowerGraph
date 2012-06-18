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



#ifndef SGD_VERTEX_PROGRAM_HPP
#define SGD_VERTEX_PROGRAM_HPP


/**
 * \file
 * \ingroup toolkit_matrix_pvecization
 *
 * \brief This file describes the vertex program for the alternating
 * least squares (SGD) matrix pvecization algorithm.  See
 * \ref sgd_vertex_program for description of the SGD Algorithm.
 */



#include <Eigen/Dense>

#include <graphlab.hpp>

#include "eigen_serialization.hpp"


typedef Eigen::VectorXd vec_type;
typedef Eigen::MatrixXd mat_type;

/** 
 * \ingroup toolkit_matrix_pvecization
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
  /** \brief The number of times this vertex has been updated. */
  uint32_t nupdates;
  /** \brief The most recent L1 change in the pvec value */
  float residual; //! how much the latent value has changed
  /** \brief The latent pvec for this vertex */
  vec_type pvec;
  float RMSE;
  uint32_t edges;

  /** 
   * \brief Simple default constructor which randomizes the vertex
   *  data 
   */
  vertex_data() : nupdates(0), residual(1), RMSE(0), edges(0) { randomize(); } 
  /** \brief Randomizes the latent pvec */
  void randomize() { pvec.resize(NLATENT); pvec.setRandom(); }
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const { 
    arc << nupdates << residual << pvec;        
  }
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc) { 
    arc >> nupdates >> residual >> pvec;
  }
}; // end of vertex data


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


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * \brief Given a vertex and an edge return the other vertex in the
 * edge.
 */
inline graph_type::vertex_type
get_other_vertex(graph_type::edge_type& edge, 
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex


/**
 * \brief Given an edge compute the error associated with that edge
 */
double extract_l2_error(const graph_type::vertex_type & vertex) {
   if (vertex.data().edges == 0)
     return 0;
   return vertex.data().RMSE / vertex.data().edges;
} // end of extract_l2_error


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
  if(role == edge_data::TRAIN || role == edge_data::VALIDATE) strm >> obs;
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id+1000000, edge_data(obs, role)); 
  return true; // successful load
} // end of graph_loader




/**
 * \brief The gather type used to construct XtX and Xty needed for the SGD
 * update
 *
 * To compute the ALS update we need to compute the sum of 
 * \code
 *  sum: XtX = nbr.pvec.transpose() * nbr.pvec 
 *  sum: Xy  = nbr.pvec * edge.obs
 * \endcode
 * For each of the neighbors of a vertex. 
 *
 * To do this in the Gather-Apply-Scatter model the gather function
 * computes and returns a pair consisting of XtX and Xy which are then
 * added. The gather type represents that tuple and provides the
 * necessary gather_type::operator+= operation.
 *
 */
class gather_type {
public:
  /**
   * \brief Stores the current sum of nbr.pvec.transpose() *
   * nbr.pvec
   */

  /**
   * \brief Stores the current sum of nbr.pvec * edge.obs
   */
  vec_type pvec;

  /** \brief basic default constructor */
  gather_type() { }

  /**
   * \brief This constructor computes XtX and Xy and stores the result
   * in XtX and Xy
   */
  gather_type(const vec_type& X, const double y) {
    pvec = X;
  } // end of constructor for gather type

  /** \brief Save the values to a binary archive */
  void save(graphlab::oarchive& arc) const { arc << pvec; }

  /** \brief Read the values from a binary archive */
  void load(graphlab::iarchive& arc) { arc >> pvec; }  

  /** 
   * \brief Computes XtX += other.XtX and Xy += other.Xy updating this
   * tuples value
   */
  gather_type& operator+=(const gather_type& other) {
    this->pvec += other.pvec;
    return *this;
  } // end of operator+=

}; // end of gather type



/**
 * SGD vertex program type
 */ 
class sgd_vertex_program : 
  public graphlab::ivertex_program<graph_type, gather_type,
                                   graphlab::messages::sum_priority>,
  public graphlab::IS_POD_TYPE {
public:
  /** The convergence tolerance */
  static double TOLERANCE;
  static double LAMBDA;
  static double GAMMA;
  static size_t MAX_UPDATES;

  /** The set of edges to gather along */
  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 

  /** The gather function computes XtX and Xy */
  gather_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const {
    //if(edge.data().role == edge_data::TRAIN) {
   if (vertex.num_in_edges() == 0){
      vertex_type other_vertex(get_other_vertex(edge, vertex));
      vertex_type my_vertex(vertex);
      vertex_data & my_data = my_vertex.data();
      const double pred = vertex.data().pvec.dot(other_vertex.data().pvec);
      const float err = (edge.data().obs - pred);
      my_vertex.data().RMSE += err*err;
      my_vertex.data().edges++;
      if (edge.data().role == edge_data::TRAIN){
        my_vertex.data().pvec += GAMMA*(err*vertex.data().pvec - LAMBDA*other_vertex.data().pvec);
        other_vertex.data().pvec += GAMMA*(err*other_vertex.data().pvec - LAMBDA*vertex.data().pvec);
      }
    } 
    return gather_type();
  } // end of gather function

  /** apply collects the sum of XtX and Xy */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& sum) {
    // Get and reset the vertex data
    vertex_data& vdata = vertex.data(); 
    // Determine the number of neighbors.  Each vertex has only in or
    // out edges depending on which side of the graph it is located
    //if(sum.Xy.size() == 0) { vdata.residual = 0; ++vdata.nupdates; return; }
    //mat_type XtX = sum.XtX;
    //vec_type Xy = sum.Xy;
    // Add regularization
    //for(int i = 0; i < XtX.rows(); ++i) XtX(i,i) += LAMBDA; // /nneighbors;
    // Solve the least squares problem using eigen ----------------------------
    //const vec_type old_pvec = vdata.pvec;
    //vdata.pvec = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xy);
    // Compute the residual change in the pvec pvec -----------------------
    //vdata.residual = (vdata.pvec - old_pvec).cwiseAbs().sum() / XtX.rows();
    ++vdata.nupdates;
  } // end of apply
  
  /** The edges to scatter along */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  /** Scatter reschedules neighbors */  
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {
    edge_data& edata = edge.data();
    if(edata.role == edge_data::TRAIN) {
      const vertex_type other_vertex = get_other_vertex(edge, vertex);
      //const double pred = vdata.pvec.dot(other_vdata.pvec);
      //const float error = std::fabs(edata.obs - pred);
      //const double priority = (error * vdata.residual); 
      // Reschedule neighbors ------------------------------------------------
      if(other_vertex.data().nupdates < MAX_UPDATES) 
        context.signal(other_vertex, 1);
    }
  } // end of scatter function


  /**
   * \brief Signal all vertices on one side of the bipartite graph
   */
  static graphlab::empty signal_left(icontext_type& context,
                                     vertex_type& vertex) {
    if(vertex.num_out_edges() > 0) context.signal(vertex);
    return graphlab::empty();
  } // end of signal_left 

}; // end of sgd vertex program


struct error_aggregator : public graphlab::IS_POD_TYPE {
  typedef sgd_vertex_program::icontext_type icontext_type;
  typedef graph_type::edge_type edge_type;
  double train_error, validation_error;
  size_t ntrain, nvalidation;
  error_aggregator() : 
    train_error(0), validation_error(0), ntrain(0), nvalidation(0) { }
  error_aggregator& operator+=(const error_aggregator& other) {
    train_error += other.train_error;
    assert(!std::isnan(train_error));
    validation_error += other.validation_error;
    ntrain += other.ntrain;
    nvalidation += other.nvalidation;
    return *this;
  }
  static error_aggregator map(icontext_type& context, const graph_type::vertex_type& vertex) {
    error_aggregator agg;
    agg.train_error = extract_l2_error(vertex); 
    assert(!std::isnan(agg.train_error));
    agg.ntrain = vertex.data().edges;
    return agg;
  }
  static void finalize(icontext_type& context, const error_aggregator& agg) {
    ASSERT_GT(agg.ntrain, 0);
    const double train_error = std::sqrt(agg.train_error / agg.ntrain);
    assert(!std::isnan(train_error));
    context.cout() << context.elapsed_seconds() << "\t" << train_error;
    if(agg.nvalidation > 0) {
      const double validation_error = 
        std::sqrt(agg.validation_error / agg.nvalidation);
      context.cout() << "\t" << validation_error; 
    }
    context.cout() << std::endl;
  }
}; // end of error aggregator


struct prediction_saver {
  typedef graph_type::vertex_type vertex_type;
  typedef graph_type::edge_type   edge_type;
  std::string save_vertex(const vertex_type& vertex) const {
    return ""; //nop
  }
  std::string save_edge(const edge_type& edge) const {
    std::stringstream strm;
    const double prediction = 
      edge.source().data().pvec.dot(edge.target().data().pvec);
    strm << edge.source().id() << '\t' 
         << edge.target().id() << '\t'
         << prediction << '\n';
    return strm.str();
  }
}; // end of prediction_saver




#endif
