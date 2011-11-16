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
  double squared_error; //!root of mean square error
  double residual;
  size_t nupdates; //! the number of times the vertex was updated
  //constructor
  vertex_data() : 
    squared_error(std::numeric_limits<double>::max()), 
    residual(std::numeric_limits<double>::max()), 
    nupdates(0) { }
  void save(graphlab::oarchive& arc) const { }
  void load(graphlab::iarchive& arc) { }
}; // end of vertex data

struct edge_data {
  float observation, weight;
  edge_data(const float observation = 0, const float weight = 1) : 
    observation(observation), weight(weight)  { } 
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
  double residual;
public:
  als_update(double residual = 0) : residual(residual) { }
  double priority() const { return residual; }
  void operator+=(const als_update& other) { residual += other.residual; }
  void operator()(icontext_type& context) {
    vertex_data& vdata = context.vertex_data(); 
    vdata.squared_error = 0; vdata.residual = 0; ++vdata.nupdates;
    const edge_list_type out_eids = context.out_edge_ids();
    const edge_list_type in_eids  = context.in_edge_ids();
    // If there are no neighbors just return
    if(out_eids.size() + in_eids.size() == 0) return;
    // Get the number of latent dimensions
    const size_t& nlatent = context.get_global_const<size_t>("nlatent");
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
      const vertex_data& neighbor = context.const_vertex_data(neighbor_id);
      const edge_data& edata = context.const_edge_data(eid);
      // Update the X'X and X'y (eigen calls are too slow)
      // Xty += neighbor.latent*(edata.observation*edata.weight);
      // XtX += (neighbor.latent*neighbor.latent.transpose()) * edata.weight;
      for(size_t i = 0; i < nlatent; ++i) {
        Xty(i) += neighbor.latent(i)*(edata.observation*edata.weight);
        for(size_t j = 0; j < nlatent; ++j) 
          XtX(i,j) += neighbor.latent(i)*neighbor.latent(j)*edata.weight;
      }
    }
    // Add regularization
    const double& lambda = context.get_global_const<double>("lambda");
    for(size_t i = 0; i < nlatent; ++i) XtX(i,i) += (lambda);
    // Solve the least squares problem using eigen ----------------------------
    const vec old_latent = vdata.latent;
    vdata.latent = XtX.ldlt().solve(Xty);
    // Compute the residual change in the latent factor -----------------------
    vdata.residual = 0;
    for(size_t i = 0; i < nlatent; ++i)
      vdata.residual += std::fabs(old_latent(i) - vdata.latent(i));
    vdata.residual /= nlatent;
    // Update the rmse and reschedule neighbors -------------------------------
    const double tolerance = context.get_global_const<double>("tolerance");
    vdata.squared_error = 0;
    foreach(const edge_id_type eid, eids) {
      // get the neighbor id
      const vertex_id_type neighbor_id = is_in_eids? 
        context.source(eid) : context.target(eid);
      const vertex_data& neighbor = context.const_vertex_data(neighbor_id);
      const edge_data& edata = context.const_edge_data(eid);
      const double pred = vdata.latent.dot(neighbor.latent);
      const double error = std::fabs(edata.observation - pred);
      vdata.squared_error += error*error;
      // Reschedule neighbors ------------------------------------------------
      if( error > tolerance && vdata.residual > tolerance) 
        context.schedule(neighbor_id, als_update(residual));
    }
  } // end of operator()
}; // end of class user_movie_nodes_update_function




class accumulator :
  public graphlab::iaccumulator<graph_type, als_update, accumulator> {
private:
  double rmse, max_rmse, residual, max_residual;
  size_t min_updates, max_updates, total_updates;
public:
  accumulator() : 
    rmse(0), max_rmse(0), residual(0), max_residual(0),
    min_updates(-1), max_updates(0), total_updates(0) { }
  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.vertex_data();
    if((context.in_edge_ids().size() + 
        context.out_edge_ids().size()) == 0)
      return;
    rmse += vdata.squared_error;
    const size_t num_edges = context.in_edge_ids().size() +
      context.out_edge_ids().size();
    max_rmse = 
      std::max(max_rmse, 
               std::sqrt(vdata.squared_error/num_edges));
    residual += vdata.residual;
    max_residual = std::max(max_residual, vdata.residual);
    min_updates = std::min(min_updates, vdata.nupdates);
    max_updates = std::max(max_updates, vdata.nupdates);
    total_updates += vdata.nupdates;
  }
  void operator+=(const accumulator& other) { 
    rmse += other.rmse; 
    max_rmse = std::max(max_rmse, other.max_rmse);
    residual += other.residual;
    max_residual = std::max(max_residual, other.residual);
    min_updates = std::min(min_updates, other.min_updates);
    max_updates = std::max(max_updates, other.max_updates);
    total_updates += other.total_updates;
  }
  void finalize(iglobal_context_type& context) {
    std::cout 
      << std::setw(10) << sqrt(rmse /(2*context.num_edges())) << '\t'
      << std::setw(10) << max_rmse << '\t'
      << std::setw(10) << (residual / context.num_vertices()) << '\t'
      << std::setw(10) << max_residual << '\t'
      << std::setw(10) << max_updates << '\t'
      << std::setw(10) << min_updates << '\t'
      << std::setw(10) << (double(total_updates) / context.num_vertices()) 
      << std::endl;
  }
}; // end of  accumulator




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
  size_t nlatent = 10;
  double lambda = 0.065;
  size_t freq = 100000;
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
  clopts.attach_option("freq",
                       &freq, freq,
                       "The number of updates between rmse calculations");
  clopts.set_scheduler_type("sweep");
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing command line arguments." << std::endl;
    return EXIT_FAILURE;
  }

  // Setup the GraphLab execution core and load graph -------------------------
  graphlab::core<graph_type, als_update> core;
  core.set_options(clopts); // attach the command line options to the core
  typedef matrix_entry<graph_type> matrix_entry_type;
  matrix_descriptor desc;
  std::vector< matrix_entry_type > test_set;
  std::cout << "Loading graph-------------------------------" << std::endl;
  const bool success = 
    load_graph(matrix_file, format, holdout, desc, core.graph(), test_set);
  if(!success) {
    std::cout << "Error in reading file: " << matrix_file << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Randomizing initial latent factors----------" << std::endl;
  initialize_vertex_data(nlatent, core.graph());

  std::cout << "Finished initializing all vertices." << std::endl;

  // Set global variables -----------------------------------------------------
  core.add_global_const("tolerance", tolerance);
  core.add_global_const("nlatent", nlatent);
  core.add_global_const("lambda", lambda);
  core.add_sync("rmse", accumulator(), freq);


  // Run the PageRank ---------------------------------------------------------
  core.schedule_all(als_update(10000));
  std::cout 
    << std::setw(10) << "Avg RMSE" << '\t'
    << std::setw(10) << "Max RMSE" << '\t'
    << std::setw(10) << "Residual" << '\t'
    << std::setw(10) << "Max Resid" << '\t'
    << std::setw(10) << "Max Update" << '\t'
    << std::setw(10) << "Min Update" << '\t'
    << std::setw(10) << "Avg Update" << std::endl;
  const double runtime = core.start();  // Run the engine
  std::cout << "Graphlab finished, runtime: " << runtime 
            << " seconds." << std::endl;
  std::cout << "Updates executed: " << core.last_update_count() 
            << std::endl;
  std::cout << "Update Rate (updates/second): " 
            << core.last_update_count() / runtime
            << std::endl;
  
  // Output Results -----------------------------------------------------------
  double squared_error = 0, weight = 0;
  foreach(const matrix_entry_type& entry, test_set) {
    const vertex_data& v1 = core.graph().vertex_data(entry.source);
    const vertex_data& v2 = core.graph().vertex_data(entry.target);
    const double prediction = v1.latent.dot(v2.latent);
    const double error = entry.edata.observation - prediction;
    squared_error = entry.edata.weight * error * error;
    weight += entry.edata.weight;
  }
  std::cout << "Test Error: " << std::sqrt(squared_error/weight) << std::endl;


  return EXIT_SUCCESS;
} // end of main
