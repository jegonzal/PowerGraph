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



#ifndef ALS_VERTEX_PROGRAM_HPP
#define ALS_VERTEX_PROGRAM_HPP


/**
 * \file
 * \ingroup toolkit_matrix_factorization
 *
 * \brief This file describes the vertex program for the alternating
 * least squares (ALS) matrix factorization algorithm.  See
 * \ref als_vertex_program for description of the ALS Algorithm.
 */



#include <Eigen/Dense>

#include <graphlab.hpp>





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
  /** \brief The latent factor for this vertex */
  Eigen::VectorXd latent;
  /** 
   * \brief Simple default constructor which randomizes the vertex
   *  data 
   */
  vertex_data(); 
  /** \brief Randomizes the latent factor */
  void randomize();
  /** \brief Save the vertex data to a binary archive */
  void save(graphlab::oarchive& arc) const;
  /** \brief Load the vertex data from a binary archive */
  void load(graphlab::iarchive& arc);
}; // end of vertex data

/**
 * \brief The edge data stores the entry in the matrix.
 *
 * In addition the edge data also stores the most recent error estimate.
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  float obs, error;
  edge_data(const float& obs = 0);
}; // end of edge data


/**
 * \brief The graph type is defined in terms of the vertex and edge
 * data.
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;




/**
 * \brief The gather type used to construct XtX and Xty needed for the ALS
 * update
 *
 * To compute the ALS update we need to compute the sum of 
 * \code
 *  sum: XtX = nbr.latent.transpose() * nbr.latent 
 *  sum: Xy  = nbr.latent * edge.obs
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
   * \brief Stores the current sum of nbr.latent.transpose() *
   * nbr.latent
   */
  Eigen::MatrixXd XtX;

  /**
   * \brief Stores the current sum of nbr.latent * edge.obs
   */
  Eigen::VectorXd Xy;

  /** \brief basic default constructor */
  gather_type() { }

  /**
   * \brief This constructor computes XtX and Xy and stores the result
   * in XtX and Xy
   */
  gather_type(const Eigen::VectorXd& X, const double y);
  /** 
   * \brief Computes XtX += other.XtX and Xy += other.Xy updating this
   * tuples value
   */
  gather_type& operator+=(const gather_type& other);

  /** \brief Save the values to a binary archive */
  void save(graphlab::oarchive& arc) const;

  /** \brief Read the values from a binary archive */
  void load(graphlab::iarchive& arc);
}; // end of gather type


/**
 * ALS vertex program type
 */ 
class als_vertex_program : 
  public graphlab::ivertex_program<graph_type, gather_type,
                                   graphlab::messages::sum_priority>,
  public graphlab::IS_POD_TYPE {
public:
  /** The convergence tolerance */
  static double TOLERANCE;
  static double LAMBDA;
  static size_t MAX_UPDATES;

  /** The set of edges to gather along */
  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const;

  /** The gather function computes XtX and Xy */
  gather_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const;

  /** apply collects the sum of XtX and Xy */
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& sum);
  
  /** The edges to scatter along */
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const;

  /** Scatter reschedules neighbors */  
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const;

  /**
   * \brief Signal all vertices on one side of the bipartite graph
   */
  static graphlab::empty signal_left(icontext_type& context,
                                     vertex_type& vertex);
  
  

private:
  /** Since the edges are undirected we use this helper function to
      get the vertex on the other side of the edge */
  vertex_type get_other_vertex(edge_type& edge, const vertex_type& vertex) const;

}; // end of als vertex program





//=============================================================================
// Graph operations 


/**
 *  the extract error function is used to compute the error on an edge
 */
double extract_error(graph_type::edge_type edge);


/**
 * The graph loader function is a line parser used for distributed
 * graph construction.
 */
bool graph_loader(graph_type& graph, 
                  const std::string& filename,
                  const std::string& line);


#endif
