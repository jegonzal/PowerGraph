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
bool final_residual = true;
bool debug_conv_fix = false;
int fix_conv = 0;
double regularization = 0;

enum gabp_fields{
  GABP_Y = 0,
  GABP_PRIOR_PREC = 1,
  GABP_REAL_X = 2,
  GABP_CUR_MEAN = 3,
  GABP_CUR_PREC = 4,
  GABP_PREV_MEAN = 5,
  GABP_PREV_PREC = 6,
  GABP_OLD_Y = 7,
  GABP_MAX_FIELD = 8
};

int data_size = GABP_MAX_FIELD;

/** Vertex data types for Gaussian BP**/
struct vertex_data {
  vec pvec;
  double A_ii;
 /* real_type prior_mean;  //prior mean (b_i / A_ii)
  real_type prior_prec;   //prior precision P_ii = A_ii
  real_type real;   //the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
  real_type cur_mean; //intermediate value to store the mean \mu_i
  real_type cur_prec; // //current precision value P_i
  real_type prev_mean; //  mean value from previous round (for convergence detection)
  real_type prev_prec; //precision value from previous round (for convergence detection)
*/
  vertex_data(){ 
     pvec = zeros(GABP_MAX_FIELD);
     pvec[GABP_PREV_MEAN] = 1000;
   };
  
  void add_self_edge(double value) { A_ii = value; }

  void set_val(double value, int field_type) { 
    if (field_type == GABP_PRIOR_PREC)
      A_ii = value;
    else
      pvec[field_type] = value;
  }
  double get_output(int field_type){ 
    return pvec[field_type];
  }

};


struct edge_data {
  real_type weight; //edge value (nnz entry in matrix A)
  real_type mean; //message \mu_ij
  real_type prec; //message P_ij
  edge_data(real_type weight) : weight(weight), mean(0), prec(0) { }
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;

#include "../shared/math.hpp" //has to be included after vertex_data and edge_data and graph_type
#include "../shared/printouts.hpp"

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
    vdata.pvec[GABP_PREV_MEAN] = vdata.pvec[GABP_CUR_MEAN];
    vdata.pvec[GABP_PREV_PREC] = vdata.pvec[GABP_CUR_PREC];

    //initialize accumlated values
    double mu_i = vdata.pvec[GABP_Y];
    real_type J_i = vdata.A_ii + regularization;
    if (!support_null_variance) assert(J_i != 0);

    /* CALCULATE new value */
    if (debug) {
      std::cout << "entering node " << context.vertex_id()
                << " P=" << vdata.A_ii + regularization
                << " u=" << vdata.pvec[GABP_Y]
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
      vdata.pvec[GABP_CUR_MEAN] = mu_i;
      vdata.pvec[GABP_CUR_PREC] = 0;
    } else {
      assert(J_i != 0);
      vdata.pvec[GABP_CUR_MEAN] = mu_i / J_i;
      assert(vdata.pvec[GABP_CUR_MEAN] != NAN);
      vdata.pvec[GABP_CUR_PREC] = J_i;
    }
    assert(vdata.pvec[GABP_CUR_MEAN] != NAN);
    
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
  real_type relative_norm;
public:
  aggregator() : 
    relative_norm(0) { }
  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    relative_norm += std::pow(vdata.pvec[GABP_CUR_MEAN] - vdata.pvec[GABP_PREV_MEAN], 2);
    if (debug)
	std::cout << "relative norm: " <<relative_norm << std::endl;
  }
  void operator+=(const aggregator& other) { 
    relative_norm += other.relative_norm;
  }
  void finalize(iglobal_context_type& context) {
    // here we can output something as a progress monitor
              logstream(LOG_INFO) << "Relative Norm: " << relative_norm << std::endl;
    // write the final result into the shared data table
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
  clopts.attach_option("debug", &debug, debug, "Verbose mode");
  clopts.attach_option("debug_conv_fix", &debug_conv_fix, debug_conv_fix, "Verbose mode for convergence fix");
  clopts.attach_option("syncinterval", 
                       &sync_interval, sync_interval, 
                       "sync interval (number of update functions before convergen detection");
  clopts.attach_option("regularization", &regularization, regularization, 
	               "regularization added to the main diagonal");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("support_null_variance", &support_null_variance, 
		       support_null_variance, "support null precision");
  clopts.attach_option("final_residual", &final_residual, final_residual, "calc residual at the end (norm(Ax-b))");
  clopts.attach_option("fix_conv", &fix_conv, fix_conv, "Fix convergence, using XX outer loop iterations");
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
    clopts.set_scheduler_type("sweep(max_iterations=100,ordering=ascending,strict=true)");   
  }
  //./gabp --data=A_gabp --yfile=y --fix_conv=10 --regularization=1 --scheduler="sweep(max_iterations=100,ordering=ascending,strict=true)"
  else if (unittest == 2){
    datafile = "A_gabp"; yfile = "y"; xfile = "x_gabp"; sync_interval = 120; fix_conv = 10; regularization = 1;
    clopts.set_scheduler_type("sweep(max_iterations=100,ordering=ascending,strict=true)");   
  }

  // Create a core
  graphlab::core<graph_type, gabp_update> core;
  core.set_options(clopts); // Set the engine options


  std::cout << "Load matrix A" << std::endl;
  bipartite_graph_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_vector(yfile, format, matrix_info, core.graph(), GABP_Y, false);
  if (xfile.size() > 0){
    std::cout << "Load x values" << std::endl;
    load_vector(xfile, format, matrix_info, core.graph(), GABP_REAL_X, true);
  }  

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(gabp_update());

  if (sync_interval < core.graph().num_vertices()){
    sync_interval = core.graph().num_vertices(); 
    logstream(LOG_WARNING) << "Sync interval is lower than the number of nodes: setting sync interval to " 
			   << sync_interval << std::endl;
  }
   

  aggregator acum;
  core.add_aggregator("sync", acum, sync_interval);
  core.add_global("RELATIVE_NORM", double(0));
  core.add_global("THRESHOLD", threshold); 

  math_info mi;
  graphlab::core<graph_type, Axb> tmp_core;
  tmp_core.graph() = core.graph();
  init_math(&tmp_core.graph(), &tmp_core, matrix_info);
  DistMat A(matrix_info);
  DistVec b(matrix_info, GABP_Y,true, "b");
  DistVec p(matrix_info, GABP_PREV_MEAN,true, "p");
  DistVec x(matrix_info, GABP_CUR_MEAN,true, "x");
  DistVec old_b(matrix_info, GABP_OLD_Y, true, "old_b");
  DistDouble pnorm;
  vec xj;
   
  int outer_iterations = 1;
  if (fix_conv > 0){
    outer_iterations = fix_conv;
    xj = zeros(matrix_info.num_nodes(false));
    old_b = b;
  }

  for (int i=0; i< outer_iterations; i++) {
    if (fix_conv > 0){
       double old_regularization = regularization;
       regularization = 0;
       p = xj;
       b = old_b - A*p;
       if (debug_conv_fix){
         PRINT_VEC(b);
         PRINT_VEC(xj);
       }
       regularization = old_regularization;
       pnorm = norm(b);
       logstream(LOG_INFO) << "Convergence fix norm of gradient is: " << pnorm.toDouble() << std::endl;
       if (pnorm < threshold)
          break;
       core.graph() = tmp_core.graph();
    }
    double runtime= core.start();
    std::cout << "GaBP instrance " << i << " finished in " << runtime << std::endl;
    tmp_core.graph() = core.graph();
    if (fix_conv > 0)
      xj += x.to_vec();
    else xj = x.to_vec();
    if (debug_conv_fix){
      PRINT_VEC(xj);
      PRINT_VEC(x); 
    }
   }

  if (final_residual){
    if (fix_conv){
       x = xj;
      regularization = 0;
      b = old_b;
    } 
    p = A*x - b;
    DistDouble ret = norm(p);
    logstream(LOG_INFO) << "Solution converged to residual: " << ret.toDouble() << std::endl;
    if (unittest == 1)
      assert(ret.toDouble() < 1e-15);
    else if (unittest == 2)
      assert(ret.toDouble() < 1e-5);
  }
 
  write_output_vector(datafile + ".curmean.out", format, xj, false, "GraphLab linear solver library. vector x, as computed by GaBP, includes the solution x = inv(A)*y");

  vec prec = fill_output(&core.graph(), matrix_info, GABP_CUR_PREC);
  write_output_vector(datafile + ".curprec.out", format, prec, false, "GraphLab linear solver library, vector prec, as computed by GaBp, includes an approximation of diag(inv(A))");
   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
#endif
