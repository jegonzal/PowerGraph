/* Copyright (c) 2009 Carnegie Mellon University. 
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
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 */





#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>


#include <graphlab.hpp>

#include "matrix_loader.hpp"





#include <graphlab/macros_def.hpp>
/**
 * Define linear algbebra types
 */
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;




/** Vertex and edge data types **/
struct vertex_data {
  vec latent; //! vector of learned values 
  float rmse; //!root of mean square error
  //constructor
  vertex_data() : rmse(0) { }
  void save(graphlab::oarchive& arc) const { 
    arc << latent.size();
    for(int i = 0; i < latent.size(); ++i) arc << latent(i);
    arc << rmse;
  }   
  void load(graphlab::iarchive& arc) {
    size_t size; 
    arc >> size; latent.resize(size);
    for(int i = 0; i < latent.size(); ++i) arc >> latent(i);
    arc >> rmse;
  }
}; // end of vertex data

struct edge_data {
  float observation, weight;
  edge_data() : observation(0), weight(1)  { } 
  void save(graphlab::oarchive& archive) const { 
    archive << observation << weight; }     
  void load(graphlab::iarchive& archive) { 
    archive >> observation >> weight; }
}; // end of edge data
typedef graphlab::graph<vertex_data, edge_data> graph_type;


/***
 * UPDATE FUNCTION
 */
class als_update : 
  public graphlab::iupdate_functor<graph_type, als_update > {
  float residual;
public:
  als_update(float residual = 0) : residual(residual) { }
  double priority() const { return residual; }
  void operator+=(const als_update& other) { residual += other.residual; }
  void operator()(icontext_type& context) {    
    vertex_data& vdata = context.vertex_data();
    vdata.rmse = 0;    // First reset the error;
    const edge_list_type out_eids = context.out_edge_ids();
    const edge_list_type in_eids  = context.in_edge_ids();
    // If there are no neighbors just return
    if(out_eids.size() + in_eids.size() == 0) return;
    // Get the number of latent dimensions
    const size_t& nlatent = context.get_global_const<size_t>("nlatent");
    // Get the local temporaries (cached for efficiency) ----------------------
    mat XtX(nlatent, nlatent); XtX.setZero();
    vec Xty(nlatent); Xty.setZero();
    // Get the non-zero edge list
    const bool is_in_eids = in_eids.size() > 0;
    const edge_list_type eids = is_in_eids? in_eids : out_eids;    
    // Compute X'X and X'y (weighted) -----------------------------------------
    foreach(const edge_id_type eid, eids) {
      // get the neighbor id
      const vertex_id_type neighbor_id = is_in_eids? 
        context.source(eid) : context.target(eid);
      const vertex_data& neighbor = 
        context.const_vertex_data(neighbor_id);
      const edge_data& edata = context.const_edge_data(eid);
      // Update the X'X and X'y
      Xty += neighbor.latent * (edata.observation * edata.weight);
      XtX += neighbor.latent.transpose() * (neighbor.latent * edata.weight);
    }
    // Add regularization
    const double& lambda = context.get_global_const<double>("lambda");
    for(int i = 0; i < nlatent; ++i) XtX(i,i) += lambda;
    // Solve the least squares problem using eigen ----------------------------
    const vec old_latent = vdata.latent;
    vdata.latent = XtX.ldlt().solve(Xty);
    // Update the rmse --------------------------------------------------------
    foreach(const edge_id_type eid, eids) {
      // get the neighbor id
      const vertex_id_type neighbor_id = is_in_eids? 
        context.source(eid) : context.target(eid);
      const vertex_data& neighbor = 
        context.const_vertex_data(neighbor_id);
      const edge_data& edata = context.const_edge_data(eid);
      const double pred = vdata.latent.dot(neighbor.latent);
      vdata.rmse += edata.weight * 
        (edata.observation - pred) * (edata.observation - pred);
    }
    // Reschedule neighbors ---------------------------------------------------
    const float change = 
      (old_latent - vdata.latent).dot(old_latent - vdata.latent);
    if(change > context.get_global_const<float>("tolerance")) {
      context.schedule_neighbors(context.vertex_id(), als_update(change));
    }    
  } // end of operator()
}; // end of class user_movie_nodes_update_function



int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_DEBUG);
  global_logger().set_log_to_console(true);
  
  // Parse command line options -----------------------------------------------
  const std::string description = 
    "Compute the ALS factorization of a matrix.";
  graphlab::command_line_options clopts(description);
  std::string matrix_file;
  std::string format = "matrixmarket";
  double tolerance = 1e-2;
  double holdout = 0.1;
  double lambda = 1;
  size_t nlatent = 10;
  clopts.attach_option("matrix",
                       &matrix_file, matrix_file,
                       "The file containing the matrix. If none is provided"
                       "then a toy matrix will be created");
  clopts.add_positional("matrix");
  clopts.attach_option("format",
                       &format, format,
                       "The matrix file format: {matrixmarket, binary}");
  clopts.attach_option("D",
                       &nlatent, nlatent,
                       "Number of latent parameters to use.");
  clopts.attach_option("lambda", &lambda, lambda, "ALS regularization weight"); 
  clopts.attach_option("tol",
                       &tolerance, tolerance,
                       "residual termination threshold");
  clopts.attach_option("holdout",
                       &holdout, holdout,
                       "The proportion of the data to use for testing.");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  // Setup the GraphLab execution core and load graph -------------------------
  graphlab::core<graph_type, als_update> core;
  core.set_options(clopts); // attach the command line options to the core
  typedef matrix_entry<graph_type> matrix_entry_type;
  std::vector< matrix_entry_type > test_set;
  const bool success = 
    load_graph(matrix_file, format, holdout, core.graph(), test_set);
  if(!success) {
    std::cout << "Error in reading file: " << matrix_file << std::endl;
  }

  // Set global variables -----------------------------------------------------
  core.add_global_const("tolerance", tolerance);
  core.add_global_const("nlatent", nlatent);
  core.add_global_const("lambda", lambda);
  // Run the PageRank ---------------------------------------------------------
  core.schedule_all(als_update(10000));
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  std::cout << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;
  
  // Output Results -----------------------------------------------------------
  double error = 0, weight = 0;
  foreach(const matrix_entry_type& entry, test_set) {
    const vertex_data& v1 = core.graph().vertex_data(entry.source);
    const vertex_data& v2 = core.graph().vertex_data(entry.target);
    const double prediction = v1.latent.dot(v2.latent);
    error = entry.edata.weight * 
      std::pow((entry.edata.observation - prediction),2);
    weight += entry.edata.weight;
  }
  std::cout << "Error: " << std::sqrt(error/weight) << std::endl;


  return EXIT_SUCCESS;
} // end of main
