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



graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::VectorXd& vec) {
  typedef Eigen::VectorXd::Index index_type;
  typedef Eigen::VectorXd::Scalar scalar_type;
  const index_type size = vec.size();
  arc << size;
  graphlab::serialize(arc, vec.data(), size * sizeof(scalar_type));
  return arc;
} // end of save vector

graphlab::iarchive& operator>>(graphlab::iarchive& arc, Eigen::VectorXd& vec) {
  typedef Eigen::VectorXd::Index index_type;
  typedef Eigen::VectorXd::Scalar scalar_type;
  index_type size = 0;
  arc >> size;
  vec.resize(size);
  graphlab::deserialize(arc, vec.data(), size * sizeof(scalar_type));
  return arc;
} // end of save vector


graphlab::oarchive& operator<<(graphlab::oarchive& arc, const Eigen::MatrixXd& mat) {
  typedef Eigen::MatrixXd::Index index_type;
  typedef Eigen::MatrixXd::Scalar scalar_type;
  const index_type rows = mat.rows();
  const index_type cols = mat.cols();
  arc << rows << cols;
  graphlab::serialize(arc, mat.data(), rows*cols*sizeof(scalar_type));
  return arc;
} // end of save matrix

graphlab::iarchive& operator>>(graphlab::iarchive& arc,  Eigen::MatrixXd& mat) {
  typedef Eigen::MatrixXd::Index index_type; 
  typedef Eigen::MatrixXd::Scalar scalar_type;
  index_type rows=0, cols=0;
  arc >> rows >> cols;
  mat.resize(rows,cols);
  graphlab::deserialize(arc, mat.data(), rows*cols*sizeof(scalar_type));
  return arc;
} // end of load matrix
