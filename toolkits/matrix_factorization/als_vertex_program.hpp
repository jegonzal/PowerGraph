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


#include <Eigen/Dense>

#include <graphlab.hpp>


/** Vertex and edge data types **/
struct vertex_data {
  static size_t NLATENT;
  uint32_t nupdates; //! the number of times the vertex was updated
  float residual; //! how much the latent value has changed
  Eigen::VectorXd latent; //! vector of learned values 
  vertex_data(); 
  void randomize();
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);
}; // end of vertex data

/**
 * The edge data is just an observation float
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  float obs;
  float error;
  uint32_t nupdates;
  edge_data(const float& obs = 0);
}; // end of edge data


/**
 * The graph type
 */ 
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * The gather type used to construct XtX and Xty needed for the ALS
 * update
 */
class gather_type {
public:
  Eigen::MatrixXd XtX;
  Eigen::VectorXd Xy;
  gather_type() { }
  gather_type(const Eigen::VectorXd& X, const double y);
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);
  gather_type& operator+=(const gather_type& other);
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

  /** Initialize the vertex data */
  void init(icontext_type& context, vertex_type& vertex);

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
