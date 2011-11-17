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

#define NO_CG_SUPPRT //comment this flag if you would like to compile cg code

#include <cmath>
#include <cstdio>
#include "linear.h"

const char * countername[] = {"EDGE_TRAVERSAL", "NODE_TRAVERSAL", "RECOMPUTE_EXP_AX_LOGREG", "GAMP_MULT_A", "GAMP_MULT_AT", "GAMP_VEC_A", "GAMP_VEC_AT"};
const char* runmodesnames[]= {"GaBP", "Jacobi", "Conjugate Gradient", "GaBP inverse", "Least Squares", "Shotgun Lasso", "Shotgun Logreg", "Generalized Approximate Message Passing (GAMP)"};


graphlab::glshared<double> REAL_NORM_KEY;
graphlab::glshared<double> RELATIVE_NORM_KEY;
graphlab::glshared<size_t> ITERATION_KEY;
graphlab::glshared<double> THRESHOLD_KEY;
graphlab::glshared<bool> SUPPORT_NULL_VARIANCE_KEY;
graphlab::glshared<bool> ROUND_ROBIN_KEY;
graphlab::glshared<bool> DEBUG_KEY;
graphlab::glshared<size_t> MAX_ITER_KEY;
graphlab::glshared<int> MATRIX_WIDTH_KEY;




#include "gabp.hpp"
#include "jacobi.hpp"
#include "cg.hpp"
#include "math.hpp"
#include "syncs.hpp"
#include "gabp_inv.hpp"
#include "io.hpp"
#include "advanced_config.h"


#include <graphlab/macros_def.hpp>

//function declerations

void unittest(command_line_options & clopts);
void verify_unittest_result(double diff);
void compute_logreg(gl_types_shotgun::core & glcore); 
void solveLasso(gl_types_shotgun::core & glcore); 
void gamp_main(gl_types_gamp::core &glcore);


// global variables

advanced_config config;
problem_setup ps;





template <typename coretype>
double start_inv(graphlab::command_line_options &clopts, advanced_config &config){

  assert(config.algorithm == GaBP_INV);
  // Create a core
  coretype core;
  core.set_engine_options(clopts); // Set the engine options

  // Load the graph --------------------------------------------------
  load_data<gl_types_inv::graph, vertex_data_inv, edge_data_inv>(&core.graph());

  // CREATE INITIAL TASKS ******
  // TOOD is this the correct starting priority?
  double initial_priority = 1.0;
  core.add_task_to_all(gabp_update_inv_function, initial_priority); 
  
  // Create an atomic entry to track iterations (as necessary)
  ITERATION_KEY.set(0);
  // Set all cosntants
  THRESHOLD_KEY.set(config.threshold);
  SUPPORT_NULL_VARIANCE_KEY.set(config.support_null_variance);
  DEBUG_KEY.set(config.debug);
  MAX_ITER_KEY.set(config.iter);
  MATRIX_WIDTH_KEY.set(ps.n);


  // START GRAPHLAB *****
  double runtime = core.start();
  // POST-PROCESSING *****
  std::cout << runmodesnames[config.algorithm] << " finished in " << runtime << std::endl;

  mat ret = zeros(core.graph().num_vertices(), core.graph().num_vertices());
  for (size_t i = 0; i < core.graph().num_vertices(); i++){
    const vertex_data_inv& vdata = core.graph().vertex_data(i);
    
    for (int j=0; j< (int)core.graph().num_vertices(); j++){
      set_val(ret, i,j,vdata.cur_mean[j]);
    }   

   //diff += pow(vdata.real - vdata.cur_mean,2);
       //x[i] = vdata.cur_mean;
       //prec[i] = vdata.cur_prec;
       //TODO
  }

    if (config.debug){
      std::cout<<ret<<std::endl;
      if (config.unittest == 6){
         mat A = init_mat("1.7011078828 0.4972882085 1.01358835; 0.4972882085 2.0077549 1.09088855; 1.01358835 1.09088855 2.690904500", 3, 3);
         std::cout<<A*ret<<std::endl;
      }
         
   }

  fill_output(&core.graph());

  write_output();

  //TODO, measure qulity of solutio
  return 0;
}

template <typename coretype>
double start_shotgun(graphlab::command_line_options &clopts, advanced_config &config){

  assert(config.algorithm == SHOTGUN_LOGREG || config.algorithm == SHOTGUN_LASSO);
  // Create a core
  coretype core;
  core.set_engine_options(clopts); // Set the engine options

  // Load the graph --------------------------------------------------
  load_data<gl_types_shotgun::graph, vertex_data_shotgun, edge_data_shotgun>(&core.graph());

#ifdef HAS_ITPP
  if (config.algorithm == SHOTGUN_LOGREG)
  compute_logreg(core);
  else solveLasso(core);
#else
  logstream(LOG_ERROR) << "itpp must be installed for running this algorithm!" << std::endl;
#endif

  // START GRAPHLAB *****
  //double runtime = core.start();
  // POST-PROCESSING *****
  std::cout << runmodesnames[config.algorithm] << " finished in " << runtime << std::endl;

  fill_output(&core.graph());

  write_output();

  return 0;
}


template <typename coretype>
double start_gamp(graphlab::command_line_options &clopts, advanced_config &config){

  assert(config.algorithm == GAMP);
  // Create a core
  coretype core;
  core.set_engine_options(clopts); // Set the engine options

  // Load the graph --------------------------------------------------
  load_data<gl_types_gamp::graph, vertex_data_gamp, edge_data_gamp>(&core.graph());

#ifdef HAS_ITPP
  gamp_main(core);
#else
  logstream(LOG_ERROR) << "itpp must be installed for running this algorithm!" << std::endl;
#endif

  std::cout << runmodesnames[config.algorithm] << " finished in " << runtime << std::endl;

  return 0;
}



template <typename coretype>
double start(graphlab::command_line_options &clopts, advanced_config &config){


  assert(config.algorithm != GaBP_INV);
  // Create a core
  coretype core;
  core.set_engine_options(clopts); // Set the engine options

  load_data<graph_type, vertex_data, edge_data>(&core.graph());
  // CREATE INITIAL TASKS ******
  // TOOD is this the correct starting priority?
  double initial_priority = 1.0;
  switch(config.algorithm){
    case GaBP: 
	  core.add_task_to_all(gabp_update_function, initial_priority); break;

    case JACOBI:
          core.add_task_to_all(jacobi_update_function, initial_priority); break;

    case CONJUGATE_GRADIENT:
          //deliberately empty, will be done later
          break;
   
    default:
         logstream(LOG_ERROR) << "Unknown algorithm" << std::endl;
         clopts.print_description(); 
         return EXIT_FAILURE;
  } 

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
}



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

#ifdef ITPP
  logstream(LOG_WARNING) << "it++ detected." << std::endl;
#elif defined(HAS_EIGEN)
  logstream(LOG_WARNING) << "Eigen detected." << std::endl;
#endif
  logstream(LOG_INFO) << "GraphLab Linear solver library code by Danny Bickson, CMU" << std::endl <<
                         "Send comments and bug reports to danny.bickson@gmail.com" << std::endl <<
                         "Currently implemented algorithms are: Gaussian Belief Propagation, Jacobi method, Conjugate Gradient" << std::endl;

  config.init_command_line_options(clopts); 
// Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  unittest(clopts);
  if (!config.unittest)
    ps.verify_setup(clopts);

  //run graphlab
  double diff = 0;
  switch(config.algorithm){
    case GaBP:
    case JACOBI:
    case CONJUGATE_GRADIENT:
       diff=start<gl_types::core>(clopts,config);
       break;

    case GaBP_INV:
       diff=start_inv<gl_types_inv::core>(clopts,config);
       break;

    case SHOTGUN_LOGREG:
    case SHOTGUN_LASSO:
       diff=start_shotgun<gl_types_shotgun::core>(clopts, config);
       break;

    case GAMP:
       diff=start_gamp<gl_types_gamp::core>(clopts, config);
       break;

    default: assert(false);
  }
 
  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (ps.counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], ps.counter[i]); 
   }

   verify_unittest_result(diff);
  
   return EXIT_SUCCESS;
}

#include <graphlab/macros_undef.hpp>

