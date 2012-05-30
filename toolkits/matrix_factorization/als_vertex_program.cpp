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

void vertex_data::randomize()  { latent.resize(NLATENT); latent.setRandom(); }

void vertex_data::save(graphlab::oarchive& arc) const { 
  arc << residual << neighborhood_total << nupdates << latent;        
}
void vertex_data::load(graphlab::iarchive& arc) { 
  arc >> residual >> neighborhood_total >> nupdates >> latent;
}

gather_type::gather_type(const Eigen::VectorXd& X, const double y) :
  XtX(X.size(), X.size()), Xy(X.size()) {
  XtX.triangularView<Eigen::Upper>() = X * X.transpose();
  Xy = X * y;
} // end of constructor for gather type


void gather_type::save(graphlab::oarchive& arc) const { arc << XtX << Xy; }
void gather_type::load(graphlab::iarchive& arc) { arc >> XtX >> Xy; }  
gather_type& gather_type::operator+=(const gather_type& other) {
  XtX.triangularView<Eigen::Upper>() += other.XtX;
  Xty += other.Xty;
  return *this;
} // end of operator+=


gather_type als_vertex_program::
gather(icontext_type& context, const vertex_type& vertex, 
       edge_type& edge) const {
  const vertex_type other_vertex = get_other_vertex(edge, vertex);
  edge_data& edata = edge.data();
  return gather_type(vertex.data().latent, edge.data().obs);
}


