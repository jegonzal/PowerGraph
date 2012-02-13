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
 * Functionality: The code solves the linear system Ax = b using
 * The Jacobi algorithm. (A is a square matrix). 
 * A assumed to be full column rank.  Algorithm is described
 * http://en.wikipedia.org/wiki/Jacobi_method
 * Written by Danny Bickson 
 */


#include <cmath>
#include <cstdio>
#include <limits>
#include <iostream>
#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
using namespace graphlab;

enum jacobi_fields{
  JACOBI_X = 0,
  JACOBI_REAL_X = 1,
  JACOBI_Y = 2,
  JACOBI_PREV_X = 3,
  JACOBI_PREC = 4
};

int data_size = 5;
bool debug = false;
double regularization = 0;
bool final_residual = true;
int max_iter = 10;


struct vertex_data {
  vec pvec;
  double A_ii;
  //real_type y, Aii;
  //real_type pvec[JACOBI_X], pvec[JACOBI_REAL_X], pvec[JACOBI_PREV_X];
  vertex_data(): A_ii(1) { //: y(0), Aii(1), pvec[JACOBI_X](0), pvec[JACOBI_REAL_X](0), 
                 // pvec[JACOBI_PREV_X](-1) 
     pvec = zeros(data_size);
     pvec[JACOBI_PREV_X] = -1;
  }
  void add_self_edge(double value) { A_ii = value + regularization; }

  void set_val(double value, int field_type) { 
     pvec[field_type] = value;
  }  
  double get_output(int field_type){ return pvec[field_type]; }
}; // end of vertex_data

struct edge_data {
  real_type weight;
  edge_data(double weight = 0) : weight(weight) { }
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
#include "../shared/math.hpp"


class aggregator :
  public graphlab::iaggregator<graph_type, Axb, aggregator> {
private:
  double real_norm, relative_norm;
public:
  aggregator() : 
    real_norm(0), 
    relative_norm(0) { }
  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    assert(!std::isnan(real_norm));
    real_norm += std::pow(vdata.pvec[JACOBI_X] - vdata.pvec[JACOBI_REAL_X],2);
    assert(!std::isnan(real_norm));
    relative_norm += std::pow(vdata.pvec[JACOBI_X] - vdata.pvec[JACOBI_PREV_X], 2);
    if (debug)
	std::cout << "Real_norm: " << real_norm << "relative norm: " <<relative_norm << std::endl;
  }
  void operator+=(const aggregator& other) { 
    assert(!std::isnan(real_norm));
    real_norm += other.real_norm; 
    assert(!std::isnan(real_norm));
    relative_norm += other.relative_norm;
  }
  void finalize(iglobal_context_type& context) {
    // here we can output something as a progress monitor
    std::cout << "Real Norm:     " << real_norm << std::endl
              << "Relative Norm: " << relative_norm << std::endl;
    // write the final result into the shared data table
    assert(!std::isnan(real_norm));
    context.set_global("REAL_NORM", real_norm);
    context.set_global("RELATIVE_NORM", relative_norm);
    const real_type threshold = context.get_global<real_type>("THRESHOLD");
    if(relative_norm < threshold) 
      context.terminate();
  }
}; // end of  aggregator

void verify_values(int unittest, double residual){
   if (unittest == 1)
     assert(residual < 1e-5);
}


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
  clopts.attach_option("final_residual", &final_residual, final_residual, "calc residual at the end (norm(Ax-b))");
  clopts.attach_option("max_iter", &max_iter, max_iter, "max number of iterations");
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
  graphlab::core<graph_type, Axb> core;
  core.set_options(clopts); // Set the engine options

  //unit testing
  if (unittest == 1){
/*A= [  1.8147    0.9134    0.2785
       0.9058    1.6324    0.5469
       0.1270    0.0975    1.9575 ];

  y= [ 0.9649    0.1576    0.9706 ]';
  x= [ 0.6803   -0.4396    0.4736 ]';
*/
    datafile = "A"; yfile = "y"; xfile = "x"; sync_interval = 120;
  }

  std::cout << "Load matrix A" << std::endl;
  bipartite_graph_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_vector(yfile, format, matrix_info, core.graph(), JACOBI_Y, false);
  std::cout << "Load x values" << std::endl;
  load_vector(xfile, format, matrix_info, core.graph(), JACOBI_REAL_X, true);
  
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

  init_math(&core.graph(), &core, matrix_info);
  DistMat A(matrix_info);
  DistVec b(matrix_info, JACOBI_Y,true, "b");
  DistVec x(matrix_info, JACOBI_X,true, "x", JACOBI_PREV_X);
  DistVec A_ii(matrix_info, JACOBI_PREC, true, "A_ii");
  A_ii = diag(A);
  //A.set_use_diag(false);
  for (int i=0; i< info.total(); i++)
    core.graph().vertex_data(i).A_ii = 0;

  for (int i=0; i < max_iter; i++){
    x = (b - A*x)/A_ii;
  }

 
  std::cout << "Jacobi finished in " << runtime << std::endl;
  std::cout << "\t Updates: " << core.last_update_count() << " per node: " 
     << core.last_update_count() / core.graph().num_vertices() << std::endl
    << "\t Rate:     " << (core.last_update_count()/runtime) << std::endl;


  if (final_residual){
    graphlab::core<graph_type, Axb> tmp_core;
    tmp_core.graph() = core.graph();
    init_math(&tmp_core.graph(), &tmp_core, matrix_info);
    DistMat A(matrix_info);
    DistVec b(matrix_info, JACOBI_Y,true,"b");
    DistVec x(matrix_info, JACOBI_X,true,"x");
    DistVec p(matrix_info, JACOBI_PREV_X, true, "p");
    for (int i=0; i< info.total(); i++)
      tmp_core.graph().vertex_data(i).A_ii = A_ii[i];

    p = A*x -b;
    DistDouble ret = norm(p);
    logstream(LOG_INFO) << "Solution converged to residual: " << ret.toDouble() << std::endl;
    if (unittest > 0){
     verify_values(unittest, ret.toDouble());
    }
  }
 
  vec ret = fill_output(&core.graph(), matrix_info, JACOBI_X);
  write_output_vector(datafile + "x.out", format, ret, false);


   return EXIT_SUCCESS;
}



