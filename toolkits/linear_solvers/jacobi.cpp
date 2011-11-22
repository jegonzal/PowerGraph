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

#include <graphlab/macros_def.hpp>

bool debug = false;

typedef double real_type;

struct vertex_data {
  real_type y, Aii;
  real_type pred_x, real_x, prev_x;
  vertex_data() : y(0), Aii(1), pred_x(0), real_x(0), 
                  prev_x(std::numeric_limits<real_type>::max()) {
    if(debug) std::cout << "hello" << std::endl;
  }
  void add_self_edge(double value) { Aii = value; }
  void update_value(double value) { y = value; }
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
  void operator()(context_type& context) {
    /* GET current vertex data */
    vertex_data& vdata = context.vertex_data();
    edge_list_type outedgeid = context.out_edge_ids();

    //store last round values
    vdata.prev_mean = vdata.cur_mean;

    //initialize accumlated values in x_i
    real_type& x_i = vdata.x;
    const real_type& A_ii = vdata.Aii;
    assert(A_ii != 0);

    if (debug) 
      std::cout << "entering node " << context.vertex() 
                << " P=" << vdata.prior_prec 
                << " u=" << vdata.prior_mean << std::endl;
  
    for(size_t i = 0; i < outedgeid.size(); ++i) {
      edge_data& out_edge = context.edge_data(outedgeid[i]);
      const vertex_data & other = 
        context.neighbor_vertex_data(context.target(outedgeid[i]));
      x_i -= out_edge.weight * other.cur_mean;
    }
    x_i /= A_ii;
    if (debug)
      std::cout << context.vertex()<< ") x_i: " << x_i << std::endl;
  }
} // end of update_functor





class accumulator :
  public graphlab::iaccumulator<graph_type, jacobi_update, accumulator> {
private:
  real_type real_norm, relative_norm;
public:
  accumulator() : 
    real_norm(std::numeric_limits<real_type>::max()), 
    relative_norm(std::numeric_limits<real_type>::max()) { }
  void operator()(icontext_type& context) {
    const vertex_data& vdata = context.const_vertex_data();
    real_norm += std::pow(vdata.pred_x - vdata.real_x,2);
    relative_norm += std::pow(vdata.pred_x - vdata.prev_x, 2);
  }
  void operator+=(const accumulator& other) { 
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
    if(real_norm < threshold) context.terminate();
  }
}; // end of  accumulator






int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile;
  std::string format = "mm";
  real_type threshold = 1e-5;
  size_t sync_interval = 10000;
  clopts.attach_option("data", &datafile, 
                       "Binary input file (as created by the save_c_gl.m script)");
  clopts.add_positional("data");
  clopts.attach_option("threshold", &threshold, threshold, "termination threshold.");
  clopts.add_positional("threshold");
  clopts.attach_option("format", &format, format, "matrix format");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("syncinterval", 
                       &sync_interval, sync_interval, 
                       "sync interval (number of update functions before convergen detection");

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

  std::cout << "Load Graph" << std::endl;
  matrix_descriptor matrix_info;
  load_graph(datafile, format, matrix_info, core.graph());
  std::cout << "Load Y values" << std::endl;
  load_y_values(yfile, format, matrix_info, core.graph());

  std::cout << "Schedule all vertices" << std::endl;
  core.schedule_all(jacobi_update());

 // Initialize the shared data --------------------------------------
  // Set syncs
  //
  switch(config.algorithm){
     case JACOBI:
     case GaBP:
      if (config.syncinterval > 0){
       core.set_sync(REAL_NORM_KEY,
                gl_types::glshared_sync_ops::sum<double, get_real_norm>,
                apply_func_real,
                double(0),  config.syncinterval,
                gl_types::glshared_merge_ops::sum<double>);

        core.set_sync(RELATIVE_NORM_KEY,
                gl_types::glshared_sync_ops::sum<double, get_relative_norm>,
                apply_func_relative,
                double(0),  config.syncinterval,
                gl_types::glshared_merge_ops::sum<double>);

  	core.engine().add_terminator(termination_condition);
     }
  }
  // Create an atomic entry to track iterations (as necessary)
  ITERATION_KEY.set(0);
  // Set all cosntants
  THRESHOLD_KEY.set(config.threshold);
  SUPPORT_NULL_VARIANCE_KEY.set(config.support_null_variance);
  DEBUG_KEY.set(config.debug);
  MAX_ITER_KEY.set(config.iter);

  // START GRAPHLAB *****
  double runtime;
  double diff = 0;

  switch(config.algorithm){
      case GaBP:
      case JACOBI:
        runtime= core.start();
        break;

      case CONJUGATE_GRADIENT:
        runtime = cg(&core,ps.means,diff,config);
        break;
  }
  // POST-PROCESSING *****
  std::cout << runmodesnames[config.algorithm] << " finished in " << runtime << std::endl;

  fill_output(&core.graph());

  write_output();

   return diff; 
  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (ps.counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], ps.counter[i]); 
   }

   verify_unittest_result(diff);
  
   return EXIT_SUCCESS;
}



#include <graphlab/macros_undef.hpp>
#endif
