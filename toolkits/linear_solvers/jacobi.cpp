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
double regularization = 0;

enum jacobi_fields{
  JACOBI_REAL_X = 1,
  JACOBI_Y = 2
};

enum jacobi_output_fields{
  JACOBI_X = 1
};

struct vertex_data {
  real_type y, Aii;
  real_type pred_x, real_x, prev_x;
  vertex_data() : y(0), Aii(1), pred_x(0), real_x(0), 
                  prev_x(-1) {
    if(debug) std::cout << "hello" << std::endl;
  }
  void add_self_edge(double value) { Aii = value + regularization; }

  void set_val(double value, int field_type) { 
     if (field_type == JACOBI_REAL_X) 
       y = value; 
     else if (field_type == JACOBI_Y)
       real_x = value;
  }
  //only one output for jacobi - solution x
  double get_output(int field_type){ return pred_x; }
}; // end of vertex_data

struct edge_data {
  real_type weight;
  edge_data(double weight = 0) : weight(weight) { }
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;

/***
 * JACOBI UPDATE FUNCTION
 * x_i = (b_i - \sum_j A_ij * x_j)/A_ii
 */
struct jacobi_update :
  public graphlab::iupdate_functor<graph_type, jacobi_update> {
  void operator()(icontext_type& context) {
    /* GET current vertex data */
    vertex_data& vdata = context.vertex_data();
    const edge_list_type out_edges = context.out_edges();

    //store last round values
    vdata.prev_x = vdata.pred_x;

    //initialize accumlated values in x_i
    real_type x_i = vdata.y;
    assert(!std::isnan(x_i));
    const real_type A_ii = vdata.Aii;

    if (debug) 
      std::cout << "entering node " << context.vertex_id() 
                << " A_ii=" << vdata.Aii 
                << " prev_x=" << vdata.prev_x 
                << " y: " << vdata.y << std::endl;
  
    for(size_t i = 0; i < out_edges.size(); ++i) {
      edge_data& out_edge = context.edge_data(out_edges[i]);
      const vertex_data & other = 
        context.const_vertex_data(out_edges[i].target());
      double val = out_edge.weight * other.pred_x;
      x_i = x_i - val;
      assert(!std::isnan(x_i));
    }
    assert(A_ii != 0);
    vdata.pred_x = x_i / A_ii;
    assert(!std::isnan(vdata.pred_x));   
 
    if (debug)
      std::cout << context.vertex_id()<< ") x_i: " << x_i << std::endl;

    assert(!std::isnan(x_i));
    context.schedule(context.vertex_id(), *this);
  }
}; // end of update_functor





class accumulator :
  public graphlab::iaccumulator<graph_type, jacobi_update, accumulator> {
private:
  double real_norm, relative_norm;
public:
  accumulator() : 
    real_norm(0), 
    relative_norm(0) { }
  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    assert(!std::isnan(real_norm));
    real_norm += std::pow(vdata.pred_x - vdata.real_x,2);
    assert(!std::isnan(real_norm));
    relative_norm += std::pow(vdata.pred_x - vdata.prev_x, 2);
    if (debug)
	std::cout << "Real_norm: " << real_norm << "relative norm: " <<relative_norm << std::endl;
  }
  void operator+=(const accumulator& other) { 
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
}; // end of  accumulator





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
  graphlab::core<graph_type, jacobi_update> core;
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
  matrix_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_vector(yfile, format, matrix_info, core.graph(), JACOBI_REAL_X, false);
  std::cout << "Load x values" << std::endl;
  load_vector(xfile, format, matrix_info, core.graph(), JACOBI_Y, true);
  

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(jacobi_update());

  if (sync_interval < core.graph().num_vertices()){
    sync_interval = core.graph().num_vertices(); 
    logstream(LOG_WARNING) << "Sync interval is lower than the number of nodes: setting sync interval to " 
			   << sync_interval << std::endl;
  }
 
  accumulator acum;
  core.add_sync("sync", acum, sync_interval);
  core.add_global("REAL_NORM", double(0));
  core.add_global("RELATIVE_NORM", double(0));
  core.add_global("THRESHOLD", threshold); 
 
  double runtime= core.start();
  // POST-PROCESSING *****
 
  std::cout << "Jacobi finished in " << runtime << std::endl;
  std::cout << "\t Updates: " << core.last_update_count() << " per node: " 
     << core.last_update_count() / core.graph().num_vertices() << std::endl
    << "\t Rate:     " << (core.last_update_count()/runtime) << std::endl;

  vec ret = fill_output(&core.graph(), matrix_info, JACOBI_X);

  write_output_vector(datafile + "x.out", format, ret);


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
