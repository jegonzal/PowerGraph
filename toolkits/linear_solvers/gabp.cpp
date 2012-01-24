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
 * GRAPHLAB implementation of Gaussiabn Belief Propagation Code See
 * algrithm description and explanation in: Danny Bickson, Gaussian
 * Belief Propagation: Theory and Application. Ph.D. Thesis. The
 * Hebrew University of Jerusalem. Submitted October 2008.
 * http://arxiv.org/abs/0811.2518 By Danny Bickson, CMU. Send any bug
 * fixes/reports to bickson@cs.cmu.edu Code adapted to GraphLab by
 * Joey Gonzalez, CMU July 2010
 *
 * Functionality: The code solves the linear system Ax = b using
 * Gaussian Belief Propagation. (A is either square matrix or
 * skinny). A assumed to be full column rank.  Algorithm is described
 * in Algorithm 1, page 14 of the above Phd Thesis.
 *
 * If you are using this code, you should cite the above reference. Thanks!
 */

#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <cmath>
#include <cstdio>
#include <limits>
#include <iostream>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;

#include <graphlab/macros_def.hpp>

bool debug = false;
bool support_null_variance = false;
double regularization = 0;

enum gabp_fields{
  GABP_REAL_X = 1,
  GABP_Y = 2
};

enum gabp_output_felds{
  GABP_CUR_MEAN = 1, 
  GABP_CUR_PREC = 2
};

/** Vertex data types for Gaussian B{**/
struct vertex_data {
  real_type prior_mean;  //prior mean (b_i / A_ii)
  real_type prior_prec;   //prior precision P_ii = A_ii
  real_type real;   //the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
  real_type cur_mean; //intermediate value to store the mean \mu_i
  real_type cur_prec; // //current precision value P_i
  real_type prev_mean; //  mean value from previous round (for convergence detection)
  real_type prev_prec; //precision value from previous round (for convergence detection)

  vertex_data(){ 
     prev_mean = 1000;
     prior_mean = prior_prec = real = cur_mean = cur_prec = prev_prec = 0;
   };
  
  void add_self_edge(double value) { prior_prec = value + regularization; }

  void set_val(double value, int field_type) { 
     if (field_type == GABP_Y){
        if (!support_null_variance)
	  assert(prior_prec > 0);
        prior_mean = value; 
     }
     else if (field_type == GABP_REAL_X)
       real = value;
     else assert(false);
  }
  double get_output(int field_type){ 
     if (field_type == GABP_CUR_PREC)
       return cur_prec;
     else if (field_type == GABP_CUR_MEAN)
       return cur_mean;
     else assert(false);    
  }

};


struct edge_data {
  real_type weight; //edge value (nnz entry in matrix A)
  real_type mean; //message \mu_ij
  real_type prec; //message P_ij
  edge_data(real_type weight) : weight(weight), mean(0), prec(0) { }
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;

/***
 * GABP UPDATE FUNCTION
 */
struct gabp_update :
  public graphlab::iupdate_functor<graph_type, gabp_update> {
  void operator()(icontext_type& context) {
    /* GET current vertex data */
    vertex_data& vdata = context.vertex_data();
    const edge_list_type out_edges = context.out_edges();
    const edge_list_type in_edges = context.in_edges();
    
    //store last round values
    vdata.prev_mean = vdata.cur_mean;
    vdata.prev_prec = vdata.cur_prec;

    //initialize accumlated values
    double mu_i = vdata.prior_mean;
    real_type J_i = vdata.prior_prec + regularization;
    if (!support_null_variance) assert(J_i != 0);

    /* CALCULATE new value */
    if (debug) {
      std::cout << "entering node " << context.vertex_id()
                << " P=" << vdata.prior_prec
                << " u=" << vdata.prior_mean
                << std::endl;
    }
    
    //accumlate all messages (the inner summation in section 4 of Algorithm 1)
    for(size_t i = 0; i < in_edges.size(); i++) {
      const edge_data& edata = context.edge_data(in_edges[i]);
      mu_i += edata.mean;
      J_i +=  edata.prec;
    }
    
    if (debug) {
      std::cout << context.vertex_id() << ") summing up all messages "
                << mu_i << " " << J_i << std::endl;
    }
    
    // optional support for null variances
    if (support_null_variance && J_i == 0){
      vdata.cur_mean = mu_i;
      vdata.cur_prec = 0;
    } else {
      assert(J_i != 0);
      vdata.cur_mean = mu_i / J_i;
      assert(vdata.cur_mean != NAN);
      vdata.cur_prec = J_i;
    }
    assert(vdata.cur_mean != NAN);
    
    /* SEND new value and schedule neighbors */
    for(size_t i = 0; i < in_edges.size(); ++i) {
      assert(in_edges[i].source() == out_edges[i].target());
      edge_data& in_edge = context.edge_data(in_edges[i]);
      edge_data& out_edge = context.edge_data(out_edges[i]);
      //graphlab::vertex_id_type target = context.target(outedgeid[i]);
      
      //substruct the sum of message sent from node j
      const real_type mu_i_j = mu_i - in_edge.mean;
      const real_type J_i_j  = J_i - in_edge.prec;
      
      if (!support_null_variance)  assert(J_i_j != 0);
      assert(out_edge.weight != 0);
      
      if (support_null_variance && J_i_j == 0){
        out_edge.mean = 0;
        out_edge.prec = 0;
      } else {
        //compute the update rule (Section 4, Algorithm 1)
        out_edge.mean = -(out_edge.weight * mu_i_j / J_i_j);
        out_edge.prec = -((out_edge.weight * out_edge.weight) / J_i_j);//matrix is assumed symmetric!
      }
      context.schedule(context.vertex_id(), *this);
      if (debug) {
        std::cout << "Sending to " << out_edges[i].target() << ' '
                  << out_edge.mean << " "
                  << out_edge.prec << " wdge weight "
                  << out_edge.weight << std::endl;
      }
    }
  } // end of operator()
}; // end of update_functor





class aggregator :
  public graphlab::iaggregator<graph_type, gabp_update, aggregator> {
private:
  real_type real_norm, relative_norm;
public:
  aggregator() : 
    real_norm(0), 
    relative_norm(0) { }
  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    real_norm += std::pow(vdata.cur_mean - vdata.real,2);
    relative_norm += std::pow(vdata.cur_mean - vdata.prev_mean, 2);
    if (debug)
	std::cout << "Real_norm: " << real_norm << "relative norm: " <<relative_norm << std::endl;
  }
  void operator+=(const aggregator& other) { 
    real_norm += other.real_norm; 
    relative_norm += other.relative_norm;
  }
  void finalize(iglobal_context_type& context) {
    // here we can output something as a progress monitor
    std::cout << "Real Norm:     " << real_norm << std::endl
              << "Relative Norm: " << relative_norm << std::endl;
    // write the final result into the shared data table
    context.set_global("REAL_NORM", real_norm);
    context.set_global("RELATIVE_NORM", relative_norm);
    const real_type threshold = context.get_global<real_type>("THRESHOLD");
    if(relative_norm < threshold) 
      context.terminate();
  }
}; // end of  aggregator





int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile, yfile, xfile;
  std::string format = "matrixmarket";
  real_type threshold = 1e-5;
  size_t sync_interval = 10000;
  int unittest = 0;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix A input file");
  clopts.add_positional("data");
  clopts.attach_option("yfile", &yfile, yfile,
                       "vector y input file");
  clopts.attach_option("xfile", &xfile, xfile,
                       "vector x input file (optional)");
  clopts.attach_option("threshold", &threshold, threshold, "termination threshold.");
  clopts.add_positional("threshold");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("syncinterval", 
                       &sync_interval, sync_interval, 
                       "sync interval (number of update functions before convergen detection");
  clopts.attach_option("regularization", &regularization, regularization, 
	               "regularization added to the main diagonal");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("support_null_variance", &support_null_variance, 
		       support_null_variance, "support null precision");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab Linear solver library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Gaussian Belief Propagation, "
    << "Jacobi method, Conjugate Gradient" << std::endl;



  // Create a core
  graphlab::core<graph_type, gabp_update> core;
  core.set_options(clopts); // Set the engine options

  //unit testing
  if (unittest == 1){
 /* A = [6.4228    2.0845    2.1617
        2.0845    4.7798    1.2418
        2.1617    1.2418    5.5250];
    y = [0.9649    0.1576    0.9706]';
    x = [0.1211    -0.0565   0.1410 ]';
    prec = [    0.1996 0.2474 0.2116]';
 */
    datafile = "A_gabp"; yfile = "y"; xfile = "x_gabp"; sync_interval = 120;
  }

  std::cout << "Load matrix A" << std::endl;
  bipartite_graph_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_vector(yfile, format, matrix_info, core.graph(), GABP_Y, false);
  std::cout << "Load x values" << std::endl;
  load_vector(xfile, format, matrix_info, core.graph(), GABP_REAL_X, true);
  
  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(gabp_update());

  if (sync_interval < core.graph().num_vertices()){
    sync_interval = core.graph().num_vertices(); 
    logstream(LOG_WARNING) << "Sync interval is lower than the number of nodes: setting sync interval to " 
			   << sync_interval << std::endl;
  }
   

  aggregator acum;
  core.add_aggregator("sync", acum, sync_interval);
  core.add_global("REAL_NORM", double(0));
  core.add_global("RELATIVE_NORM", double(0));
  core.add_global("THRESHOLD", threshold); 
 
  double runtime= core.start();
  // POST-PROCESSING *****
 
  std::cout << "GaBP finished in " << runtime << std::endl;

  vec x = fill_output(&core.graph(), matrix_info, GABP_CUR_MEAN);
  write_output_vector(datafile + ".curmean.out", format, x, false);

  vec prec = fill_output(&core.graph(), matrix_info, GABP_CUR_PREC);
  write_output_vector(datafile + ".curprec.out", format, prec, false);


  if (unittest == 1){
    double real_norm = core.get_global<double>("REAL_NORM");
    double relative_norm = core.get_global<double>("RELATIVE_NORM");
    assert(real_norm < 1e-30);
    assert(relative_norm < 1e-30);
  }

   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
#endif
