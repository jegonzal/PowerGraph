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

#include "eigen_serialization.hpp"
#include "als_vertex_program.hpp"

//=============================================================================
// Vertex data
size_t vertex_data::NLATENT = 20;

vertex_data::vertex_data() : nupdates(0), residual(1) { }

void vertex_data::randomize()  { latent.resize(NLATENT); latent.setRandom(); }

void vertex_data::save(graphlab::oarchive& arc) const { 
  arc << nupdates << residual << latent;        
}
void vertex_data::load(graphlab::iarchive& arc) { 
  arc >> nupdates >> residual >> latent;
}

//=============================================================================
// Edge data

edge_data::edge_data(const float& obs) :
  obs(obs), error(std::numeric_limits<float>::max()) { }




//=============================================================================
// Gather type
gather_type::gather_type(const Eigen::VectorXd& X, const double y) :
  XtX(X.size(), X.size()), Xy(X.size()) {
  XtX.triangularView<Eigen::Upper>() = X * X.transpose();
  Xy = X * y;
} // end of constructor for gather type

void gather_type::save(graphlab::oarchive& arc) const { arc << XtX << Xy; }
void gather_type::load(graphlab::iarchive& arc) { arc >> XtX >> Xy; }  
gather_type& gather_type::operator+=(const gather_type& other) {
  XtX.triangularView<Eigen::Upper>() += other.XtX;  
  Xy += other.Xy;
  return *this;
} // end of operator+=



//=============================================================================
// Vertex program

double als_vertex_program::TOLERANCE = 1e-3;
double als_vertex_program::LAMBDA = 0.01;
size_t als_vertex_program::MAX_UPDATES = -1;

void als_vertex_program::
init(icontext_type& context, vertex_type& vertex) {
  vertex.data().randomize();
  // Schedule self if on the left side of the bipartite graph
  if(vertex.num_out_edges() > 0) context.signal(vertex);
} // end of init


graphlab::edge_dir_type als_vertex_program::
gather_edges(icontext_type& context, const vertex_type& vertex) const { 
  return graphlab::ALL_EDGES; 
}; // end of gather_edges 


gather_type als_vertex_program::
gather(icontext_type& context, const vertex_type& vertex, 
       edge_type& edge) const {
  const vertex_type other_vertex = get_other_vertex(edge, vertex);
  return gather_type(other_vertex.data().latent, edge.data().obs);
} // end of gather function


void als_vertex_program::
apply(icontext_type& context, vertex_type& vertex, const gather_type& sum) {
  // Get and reset the vertex data
  vertex_data& vdata = vertex.data(); 
  // Determine the number of neighbors.  Each vertex has only in or
  // out edges depending on which side of the graph it is located
  const size_t nneighbors = vertex.num_in_edges() + vertex.num_out_edges();
  if(nneighbors == 0) { vdata.residual = 0; ++vdata.nupdates; return; }
  Eigen::MatrixXd XtX = sum.XtX;
  Eigen::VectorXd Xy = sum.Xy;
  // Add regularization
  for(int i = 0; i < XtX.rows(); ++i) XtX(i,i) += LAMBDA; // /nneighbors;
  // Solve the least squares problem using eigen ----------------------------
  const Eigen::VectorXd old_latent = vdata.latent;
  vdata.latent = XtX.selfadjointView<Eigen::Upper>().ldlt().solve(Xy);
  // Compute the residual change in the latent factor -----------------------
  vdata.residual = (vdata.latent - old_latent).cwiseAbs().sum() / XtX.rows();
  ++vdata.nupdates;
} // end of apply


graphlab::edge_dir_type als_vertex_program::
scatter_edges(icontext_type& context, const vertex_type& vertex) const { 
  return graphlab::ALL_EDGES; 
}; // end of scatter edges


void als_vertex_program::
scatter(icontext_type& context, const vertex_type& vertex, 
        edge_type& edge) const {
  const vertex_type other_vertex = get_other_vertex(edge, vertex);
  edge_data& edata = edge.data();
  const vertex_data& vdata = vertex.data();
  const vertex_data& other_vdata = other_vertex.data();
  const double pred = vdata.latent.dot(other_vdata.latent);
  const float error = std::fabs(edata.obs - pred);
  edata.error = error;
  assert(!std::isinf(error));
  assert(!std::isnan(error));
  const double priority = (error * vdata.residual); 
  assert(!std::isinf(priority));
  assert(!std::isnan(priority));
  // Reschedule neighbors ------------------------------------------------
  if( priority > TOLERANCE && other_vdata.nupdates < MAX_UPDATES) 
    context.signal(other_vertex, priority);
} // end of scatter function


als_vertex_program::vertex_type als_vertex_program::
get_other_vertex(edge_type& edge, const vertex_type& vertex) const {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex


//=============================================================================
// Graph operations 
double extract_error(graph_type::edge_type edge) {
  const double pred = 
    edge.source().data().latent.dot(edge.target().data().latent);
  return (edge.data().obs - pred) * (edge.data().obs - pred);
} // end of extract_train_error


bool graph_loader(graph_type& graph, const std::string& filename,
                  const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  std::stringstream strm(line);
  graph_type::vertex_id_type source_id(-1), target_id(-1);
  float obs(0);
  strm >> source_id >> target_id >> obs;
  // since this is a bipartite graph I need a method to number the
  // left and right vertices differently.  To accomplish I make sure
  // all vertices have non-zero ids and then negate the right vertex.
  source_id++; target_id++;
  target_id = -target_id;
  // Create an edge and add it to the graph
  graph.add_edge(source_id, target_id, obs); 
  return true; // successful load
} // end of graph_loader


