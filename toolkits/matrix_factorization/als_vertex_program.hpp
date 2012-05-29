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
  float residual; 
  float neighborhood_total; //! sum of values on edges
  uint32_t nupdates; //! the number of times the vertex was updated
  Eigen::VectorXd latent; //! vector of learned values 
  vertex_data() : residual(std::numeric_limits<float>::max()), 
                  neighborhood_total(0), nupdates(0) { }
  void randomize();
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);
}; // end of vertex data

/**
 * The edge data is just an observation float
 */
struct edge_data : public graphlab::IS_POD_TYPE {
  float obs, error;
  edge_data(const float& obs = 0) :
    obs(obs), error(std::numeric_limits<float>::max()) { }
}; // end of edge data


/**
 * The graph type
 */ 
graphlab::distributed_graph<vertex_data, edge_data> graph_type;


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
class als_vertex_program : public ivertex_program<vertex_data,
                                                  edge_data,
                                                  gather_type> {
public:
  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of gather_edges 
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  gather_type gather(icontext_type& context, 
                     const vertex_type& vertex, 
                     edge_type& edge) const;

private:
  vertex_type get_other_vertex(edge_type& edge, 
                               const vertex_type& vertex) const {
    return vertex.id() == edge.source().id()? edge.target() : edge.source();
  }; // end of other_vertex



}; // end of als vertex program


#endif
