#include "linear.h"

#include "advanced_config.h"

extern problem_setup ps;
extern advanced_config ac;

void advanced_config::init_command_line_options(graphlab::command_line_options & clopts){
  // Setup additional command line arguments for the GABP program
  threshold = 1e-5;
  support_null_variance = false;
  iter = 0;
  syncinterval = 10000;
  unittest = 0;

  clopts.attach_option("algorithm", &algorithm, "Algorithm 0=Gaussian BP, 1= Jacobi");
  clopts.add_positional("algorithm");
 
  clopts.attach_option("data", &datafile, "Binary input file (as created by the save_c_gl.m script)");
  clopts.add_positional("data");
  clopts.attach_option("threshold", &threshold, threshold, "termination threshold.");
  clopts.add_positional("threshold");
  clopts.attach_option("nullvar", &support_null_variance, support_null_variance,
                       "(optional) support invalid covariance matrices - with null variance.");
  clopts.attach_option("square", &square, square, "is the matrix square? ");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("syncinterval", &syncinterval, syncinterval, "sync interval (number of update functions before convergen detection");
  clopts.attach_option("supportgraphlabcf", &supportgraphlabcf, supportgraphlabcf, "input is given in GraphLab collaborative filtering format");
  clopts.attach_option("float", &isfloat, isfloat, "input file is given in float format");
  clopts.attach_option("cg_resid", &cg_resid, cg_resid, "compute cg residual progress ");
  clopts.attach_option("zero", &zero, zero, "support sparse matrix entry containing zero val ");
  clopts.attach_option("unittest", &unittest, unittest, "unit testing ( allowed values: 1/2)");
  clopts.attach_option("supportnullvariance", &support_null_variance, support_null_variance, "GABP: support zero precision (infinite variance)"); 
  clopts.attach_option("regularization", &regularization, regularization, "regularization"); 
  clopts.attach_option("oldformat", &oldformat, oldformat, "support for old edge file format: [int, int, double]"); 
  clopts.attach_option("shotgun_cost", &display_cost, display_cost, "shotgun: display cost function on each iteration"); 
  clopts.attach_option("max_iter", &iter, iter, "maximal number of iterations");
  clopts.attach_option("lambda", &shotgun_lambda, shotgun_lambda, "shotgun: lambda");
  clopts.attach_option("matrixmarket", &matrixmarket, matrixmarket, "Input data is in matrix market format");
}

void problem_setup::verify_setup(graphlab::command_line_options& clopts){
  // Ensure that a data file is provided
  if (!clopts.is_set("data")) {
    logstream(LOG_ERROR)<<"No data file provided!" << std::endl;
    clopts.print_description();
    exit(1);
  }
  
 // Ensure that an algorithm is provided
  if(!clopts.is_set("algorithm")) {
    logstream(LOG_ERROR)<< "No algorithm provided! Choose from: 0) GaBP 1) Jacobi 2) CG" << std::endl;
    clopts.print_description();
    exit(1);
  }


 
}
