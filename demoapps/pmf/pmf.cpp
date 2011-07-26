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


#include <fstream>
#include <cmath>
#include <cstdio>

#include "graphlab.hpp"
#include "itppvecutils.hpp"
#include "pmf.h"
#include "prob.hpp"
#include "bptf.hpp"
#include "sgd.hpp"
#include "lanczos.hpp"
#include "nmf.hpp"
#include "als.hpp"
#include "tensor.hpp"
#include "unittest.hpp"
#include "io.hpp"
#include "stats.hpp"
#include "implicit.hpp"
#include "lasso.hpp"
#include "cosamp.hpp"

#ifdef GL_SVD_PP
#include "svdpp.hpp"
#endif

#include <graphlab/macros_def.hpp>
/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU
 See documentation in header file pmf.h

*/

using namespace graphlab;
using namespace itpp;
using namespace std;

bool BPTF = true; //is this a sampling MCMC algo?
bool tensor = true; //is this tensor or a matrix
bool debug = false; //debug mode
bool ZERO = false; //support edges with zero weight
bool loadfactors = false; //start from a previous solution instead of a random point
bool savefactors = true; //save solution to file
bool loadgraph = false; //load graph from a binary saved file
bool savegraph = false; //save input file into graph binary file
bool stats = false; //print out statistics and exit
bool regnormal = false; //regular normalization
bool aggregatevalidation = false; //use validation dataset as training data
bool outputvalidation = false; //experimental: output validation results of kdd format
bool binaryoutput = false; //export the factors U,V,T to a binary file
bool printhighprecision = false; //print RMSE output with high precision

double scaling = 1.0; //aggregate time values into bins? (default =1, no aggregation)
double truncating = 0.0; // truncate unused time bins (optional, default = 0, no truncation)
double scalerating = 1.0; //scale the rating by dividing to the scalerating factor (optional)
double pT = 1; //regularization for tensor time nodes
runmodes algorithm; //type of algorithm
int ialgo = 0;
timer gt;
std::string infile;
int iiter = 1;//count number of time zero node run
mat U,V,T; //for storing the output

/* Variables for PMF */
int M,N,K,L;//training size: users, movies, times, number of edges
int Le = 0; //number of ratings in validation dataset 
int Lt = 0;//number of rating in test data set
mat dp;


bool FLOAT=false; //is data in float format
double LAMBDA=1;//regularization weight

/* variables of BPTF */
int delayalpha = 0; //delay alpha sampling (optional, for BPTF)
int BURN_IN =10; //burn-in priod (for MCMC sampling - optional)

/* Variables for SVD++ */
float svdpp_step_dec = 0.9;//step decrement size for SVD++
double globalMean[3] = {0}; //store global mean of matrix/tensor entries

/* Variables for SGD */
float sgd_gamma = 1e-2; //step size
float sgd_lambda = 0.3; //starting step size
float sgd_step_dec = 0.9; //step decrement size

/* variables for SVD */
int svd_iter = 10; //number of iterations (which is the number of extracted eigenvectors)

/* implicit ratings variables (see reference 10 in pmf.h) */
float implicitratingweight = 1;
float implicitratingvalue = 0;
string implicitratingtype = "none";
float implicitratingpercentage = 0;

/* sparsity enforcing priors (see reference 11 in pmf.h) */
int lasso_max_iter = 10;
#define DEFAULT_SPARSITY 0.8
double user_sparsity = DEFAULT_SPARSITY;
double movie_sparsity= DEFAULT_SPARSITY;
//performance counters
#define MAX_COUNTER 20
double counter[MAX_COUNTER];
int unittest = 0;
vertex_data * times = NULL;

gl_types::iengine * engine;
graph_type* g;
graph_type validation_graph;
graph_type test_graph;


/* Function declerations */ 
void load_pmf_graph(const char* filename, graph_type * g, testtype flag,gl_types::core & glcore);    
void calc_T(int id);    
double calc_obj(double res);
void last_iter();
void export_kdd_format(graph_type * _g, testtype type, bool dosave);
void calc_stats(testtype type);
void lanczos(gl_types::core & glcore);
void init_lanczos();
void nmf_init();




/**
 * printout RMSE statistics after each iteration
 */
void last_iter(){
  printf("Entering last iter with %d\n", iiter);

  double res,res2;
  double rmse = (algorithm != STOCHASTIC_GRADIENT_DESCENT && algorithm != NMF) ? agg_rmse_by_movie(res) : agg_rmse_by_user(res);
  //rmse=0;
  printf(printhighprecision ? 
        "%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n":
        "%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n"
        , gt.current_time(), runmodesname[algorithm], iiter,calc_obj(res),  rmse, calc_rmse_wrapper(&validation_graph, true, res2));
  iiter++;

  if (BPTF)
    last_iter_bptf(res);        
}


void add_tasks(gl_types::core & glcore){

  std::vector<vertex_id_t> um;
  for (int i=0; i< M+N; i++)
    um.push_back(i);
 
 // add update function for user and movie nodes (tensor dims 1+2) 
  switch (algorithm){
     case ALS_TENSOR_MULT:
     case ALS_MATRIX:
     case ALS_SPARSE_USR_FACTOR:
     case ALS_SPARSE_USR_MOVIE_FACTORS:
     case ALS_SPARSE_MOVIE_FACTOR:
     case BPTF_TENSOR:
     case BPTF_TENSOR_MULT:
     case BPTF_MATRIX:
     case WEIGHTED_ALS:
       glcore.add_tasks(um, user_movie_nodes_update_function, 1);
       break;

     case STOCHASTIC_GRADIENT_DESCENT:
       glcore.add_tasks(um, sgd_update_function, 1);
       break;
    
     case SVD_PLUS_PLUS: 
       glcore.add_tasks(um, svd_plus_plus_update_function, 1);
       break;

     case LANCZOS:
     case NMF:
       //lanczos is unique since it has more than one update function
       //lanczos code is done later
       break;
 }

  // add update function for time nodes (dim 3)
  if (tensor){
    std::vector<vertex_id_t> tv;
    for (int i=M+N; i< M+N+K; i++)
      tv.push_back(i);
    glcore.add_tasks(tv, time_node_update_function, 1);
  }


}


void init(){

  if (BPTF)
     init_self_pot(); 

  switch(algorithm){
   case SVD_PLUS_PLUS:
     svd_init(); break;

   case LANCZOS: 
     init_lanczos(); break;
   
   case NMF:
      nmf_init(); break;

   case ALS_MATRIX:
   case ALS_TENSOR_MULT:
   case ALS_SPARSE_USR_FACTOR:
   case ALS_SPARSE_USR_MOVIE_FACTORS:
   case ALS_SPARSE_MOVIE_FACTOR:
   case WEIGHTED_ALS:
   case BPTF_TENSOR_MULT:
   case BPTF_MATRIX:
   case BPTF_TENSOR:
   case STOCHASTIC_GRADIENT_DESCENT:
      init_pmf(); break;
  }

}


void run_graphlab(gl_types::core &glcore,timer & gt ){
        glcore.start();
        // calculate final RMSE
        double res, train_rmse =  agg_rmse_by_movie(res), res2;
        double obj = calc_obj(res);
        double validation_rmse = calc_rmse_wrapper(&validation_graph, true, res2);
        printf(printhighprecision ? 
              "Final result. Obj=%g, TRAIN RMSE= %0.12f VALIDATION RMSE= %0.12f.\n":
              "Final result. Obj=%g, TRAIN RMSE= %0.4f VALIDATION RMSE= %0.4f.\n"
              , obj,  train_rmse, validation_rmse);
        double runtime = gt.current_time();
        printf("Finished in %lf seconds\n", runtime);
        if (unittest > 0){
            verify_result(obj, train_rmse, validation_rmse);
        }

}



/** 
 * ==== SETUP AND START
 */
void start(int argc, char ** argv) {
      
  command_line_options clopts;
  clopts.attach_option("infile", &infile, infile, "Input matrix/tensor");
  clopts.add_positional("infile");

  clopts.attach_option("algorithm", &ialgo, ialgo, "algorithm");
  clopts.add_positional("algorithm");
  
  clopts.attach_option("debug", &debug, debug, "Display debug output. (optional)");
  clopts.attach_option("float", &FLOAT, FLOAT, "is data in float format?");
  clopts.attach_option("D", &D, D, "number of features (dimension of computed weight vector)");
  clopts.attach_option("zero", &ZERO, ZERO, "support zero edges");  
  clopts.attach_option("scaling", &scaling, scaling, "time scaling factor (optional)");  
  clopts.attach_option("truncating", &truncating, truncating, "time truncation factor (optional)");  
  clopts.attach_option("savegraph", &savegraph, savegraph, "save graphs to file");  
  clopts.attach_option("loadgraph", &loadgraph, loadgraph, "load graphs to file");  
  clopts.attach_option("savefactors", &savefactors, savefactors, "save factors to file");  
  clopts.attach_option("loadfactors", &loadfactors, loadfactors, "load factors from file");  
  clopts.attach_option("stats", &stats, stats, "compute graph statistics");  
  clopts.attach_option("binaryoutput", &binaryoutput, binaryoutput, "export U,V,T to a binary file"); 
 
  //BPTF related switches
  clopts.attach_option("alpha", &alpha, alpha, "BPTF alpha (noise parameter)");  
  clopts.attach_option("burn_in", &BURN_IN, BURN_IN, "BPTF burn-in period");
  clopts.attach_option("delayalpha", &delayalpha, delayalpha, "BPTF start sampling alpha (noise level) the delayalpha round ");  
  
  //ALS related switches
  clopts.attach_option("regnormal", &regnormal, regnormal, "ALS - use identical normalization for each variable? (default is weighted regularization by the number of edges");  
  clopts.attach_option("lambda", &LAMBDA, LAMBDA, "ALS regularization weight");  
  
  clopts.attach_option("scalerating", &scalerating, scalerating, "scale rating value ");  
  clopts.attach_option("aggregatevalidation", &aggregatevalidation, aggregatevalidation, "aggregate training and validation into one dataset ");  
  clopts.attach_option("maxval", &maxval, maxval, "maximal allowed value in matrix/tensor");
  clopts.attach_option("minval", &minval, minval, "minimal allowed value in matrix/tensor");
  clopts.attach_option("outputvalidation", &outputvalidation, outputvalidation, "output prediction on vadlidation data in kdd format");
  clopts.attach_option("unittest", &unittest, unittest, "unit testing. ");

  //SVD++ related switches
  clopts.attach_option("svdpp_step_dec", &svdpp_step_dec, svdpp_step_dec, "SVD++ step decrement ");
 
  //SGD related switches
  clopts.attach_option("sgd_lambda", &sgd_lambda, sgd_lambda, "SGD step size");
  clopts.attach_option("sgd_gamma", &sgd_gamma, sgd_gamma, "SGD starting step size");
  clopts.attach_option("sgd_step_dec", &sgd_step_dec, sgd_step_dec, "SGD step decrement");
 
  //SVD related switches
  clopts.attach_option("svd_iter", &svd_iter, svd_iter, "SVD iteration number"); 
  clopts.attach_option("printhighprecision", &printhighprecision, printhighprecision, "print RMSE output with high precision");

  //implicit rating options (see reference 10 in pmf.h)
  clopts.attach_option("implicitratingtype", &implicitratingtype, implicitratingtype, "type can be: user/item/uniform");
  clopts.attach_option("implicitratingpercentage", &implicitratingpercentage, implicitratingpercentage, " precentage of implicit added edges (0-100)");
  clopts.attach_option("implicitratingvalue", &implicitratingvalue, implicitratingvalue, "value for implicit negative ratings");
  clopts.attach_option("implicitratingweight", &implicitratingweight, implicitratingweight, "weight/time for implicit negative ratings");

  //sparsity enforcing priors (see reference 11 in pmf.h)
  clopts.attach_option("user_sparsity", &user_sparsity, user_sparsity, "user sparsity [0.5->1) (for L1 regularization for sparsity enforcing priors - run modes 10,11");
  clopts.attach_option("movie_sparsity", &movie_sparsity, movie_sparsity, "movie sparsity [0.5->1) (for L1 regularization for sparsity enforcing priors - run modes 11,12");
  clopts.attach_option("lasso_max_iter", &lasso_max_iter, lasso_max_iter, "max iter for lasso sparsity (run modes 10-12)");

  assert(clopts.parse(argc, argv));
  
  if (unittest > 0)
     unit_testing(unittest,clopts);

  algorithm = (runmodes)ialgo;
  printf("Setting run mode %s\n", runmodesname[algorithm]);

  switch(algorithm){
  // iterative matrix factorization using alternating least squares
  // or SVD ++
  case ALS_MATRIX:
  case ALS_SPARSE_USR_FACTOR:
  case ALS_SPARSE_USR_MOVIE_FACTORS:
  case ALS_SPARSE_MOVIE_FACTOR:
  case WEIGHTED_ALS:
  case SVD_PLUS_PLUS:
  case STOCHASTIC_GRADIENT_DESCENT:
  case LANCZOS:
  case NMF:
    tensor = false; BPTF = false;
    break;

    // MCMC tensor factorization
  case BPTF_TENSOR:
    // tensor factorization , allow for multiple edges between user and movie in different times
  case BPTF_TENSOR_MULT:
    tensor = true; BPTF = true;
   break;
    //MCMC matrix factorization
  case BPTF_MATRIX:
    tensor = false; BPTF = true;
    break;
   // tensor factorization
  case ALS_TENSOR_MULT:
    tensor = true; BPTF = false;
    break;
  default:
    assert(0);
  }

  logger(LOG_INFO, "%s starting\n",runmodesname[algorithm]);

//INPUT SANITY CHECKS
#ifdef GL_NO_MCMC
  if (BPTF){
    logstream(LOG_ERROR) << "Can not run MCMC method with GL_NO_MCMC flag. Please comment flag on pmf.h and recompile\n";
    exit(1); 
 }
#endif

#ifdef GL_NO_MULT_EDGES
  if (algorithm == ALS_TENSOR_MULT || algorithm == BPTF_TENSOR_MULT){
    logstream(LOG_ERROR) << "Can not have support for multiple edges with GL_NO_MULT_EDGES flag. Please comment flag on pmf.h and recompile\n";
   exit(1);
  }
#endif  
#ifdef GL_SVD_PP
  if (algorithm != SVD_PLUS_PLUS){
    logstream(LOG_ERROR) << "Can not run required algorithm with GL_SVD_PP flag. Please comment flag on pmf.h and recompile\n";
    exit(1);
  }
#else
  if (algorithm == SVD_PLUS_PLUS){
    logstream(LOG_ERROR) << "Can not run required algorithm without GL_SVD_PP flag. Please define flag on pmf.h and recompile\n";
    exit(1);
  }
#endif

  if (delayalpha != 0 && (algorithm != BPTF_TENSOR_MULT && algorithm != BPTF_TENSOR))
	logstream(LOG_WARNING) << "Delaying alpha (sampling of noise level) is ignored in non-MCMC methods" << std::endl;

  if (BURN_IN != 10 && (algorithm != BPTF_TENSOR_MULT && algorithm != BPTF_TENSOR && algorithm != BPTF_MATRIX))
	logstream(LOG_WARNING) << "Markov chain burn in period is ignored in non-MCMC methods" << std::endl;

  if (user_sparsity < 0.5 || user_sparsity >= 1){
	logstream(LOG_ERROR) << "user_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
        exit(1);
  }
  if (movie_sparsity < 0.5 || movie_sparsity >= 1){
	logstream(LOG_ERROR) << "movie_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
        exit(1);
  }
   gl_types::core glcore;
  //read the training data
  printf("loading data file %s\n", infile.c_str());
  if (!loadgraph){
    g=&glcore.graph();
    load_pmf_graph(infile.c_str(), g, TRAINING, glcore);

  //read the vlidation data (optional)
    printf("loading data file %s\n", (infile+"e").c_str());
    load_pmf_graph((infile+"e").c_str(),&validation_graph, VALIDATION, glcore);

  //read the test data (optional)
    printf("loading data file %s\n", (infile+"t").c_str());
    load_pmf_graph((infile+"t").c_str(),&test_graph, TEST, glcore);


    if (savegraph){
	printf("Saving .graph files\n");
	char filename[256];
        sprintf(filename, "%s%d.graph", infile.c_str(), D);
        std::ofstream fout(filename, std::fstream::binary);
        graphlab::oarchive oarc(fout);
	oarc << M << N << K << L << Le << Lt << D;
        oarc << *g << validation_graph << test_graph;
        printf("Done!\n");
        fout.close();
	exit(0);
    }

  } else {
    char filename[256];
    sprintf(filename, "%s%d.graph", infile.c_str(), D);
    std::ifstream fin(filename, std::fstream::binary);
    graphlab::iarchive iarc(fin);
    iarc >> M >> N >> K >> L >> Le >> Lt >> D;
    printf("Loading graph from file\n");
    iarc >> glcore.graph() >> validation_graph >> test_graph;
    g=&glcore.graph();
    printf("Matrix size is: USERS %dx MOVIES %dx TIME BINS %d D=%d\n", M, N, K, D);   
    printf("Creating %d edges (observed ratings)...\n", L);
  }
  

  if (loadfactors){
     import_uvt_from_file();
  }


  if (stats){
    calc_stats(TRAINING);
    calc_stats(VALIDATION);
    calc_stats(TEST);
    exit(0);
  }

  if (algorithm == ALS_TENSOR_MULT || algorithm == ALS_MATRIX || algorithm == ALS_SPARSE_USR_FACTOR || algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || algorithm == ALS_SPARSE_MOVIE_FACTOR){
    printf("setting regularization weight to %g\n", LAMBDA);
    pU=pV=LAMBDA;
  }
  glcore.set_engine_options(clopts); 

  if (tensor)
    dp = GenDiffMat(K)*pT;
  if (debug)
    std::cout<<dp<<std::endl;

  
  add_tasks(glcore);

  
  printf("%s for %s (%d, %d, %d):%d.  D=%d\n", runmodesname[algorithm], tensor?"tensor":"matrix", M, N, K, L, D);
  
  init();

   if (infile == "netflix" || infile == "netflix-r"){
       minval = 1; maxval = 5;
    }
   else if ((infile == "kddcup" || infile == "kddcup2") && maxval == DEF_MAX_VAL){
       minval = 0; maxval = 100;
   }

   if (algorithm != LANCZOS){
     double res, res2;
     double rmse =  calc_rmse_wrapper(g, false, res);
     printf(printhighprecision ? 
           "complete. Objective=%g, TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n" :
           "complete. Objective=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n" 
           , calc_obj(res), rmse, calc_rmse(&validation_graph, true, res2));
  }
 
 
  if (BPTF){
    //sample hyper priors and noise level
    if (delayalpha < iiter)
    	sample_alpha(L);
    sample_U();
    sample_V();
    if (tensor) 
      sample_T();
  }

  g->finalize();  
  gt.start();

  /**** START GRAPHLAB AND RUN UNTIL COMPLETION *****/
    switch(algorithm){
      case ALS_TENSOR_MULT:
      case ALS_MATRIX:
      case ALS_SPARSE_USR_FACTOR:
      case ALS_SPARSE_USR_MOVIE_FACTORS:
      case ALS_SPARSE_MOVIE_FACTOR:
      case WEIGHTED_ALS:
      case BPTF_TENSOR_MULT:
      case BPTF_TENSOR:
      case BPTF_MATRIX:
      case SVD_PLUS_PLUS:
      case STOCHASTIC_GRADIENT_DESCENT:
         run_graphlab(glcore, gt);
         break;
     
     case LANCZOS:
        lanczos(glcore); break;

     case NMF:
        nmf(&glcore); break;
  }

 if (algorithm != LANCZOS){
    /**** OUTPUT KDD FORMAT *****/
    if (infile == "kddcup" || infile == "kddcup2"){
      if (outputvalidation) //experimental: output prediction of validation data
 	     export_kdd_format(&validation_graph, VALIDATION, true);
      else //output prediction of test data, as required by KDD 
	     export_kdd_format(&test_graph, TEST, true);
    }
 }

  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], counter[i]); 
  }

  //write output matrices U,V,T to file
  if (binaryoutput)
     export_uvt_to_binary_file();
  else // it++ output
   export_uvt_to_itpp_file();
}



//main function 
int main(int argc,  char *argv[]) {

  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logstream(LOG_INFO)<< "PMF/BPTF/ALS/SVD++/SGD/SVD Code written By Danny Bickson, CMU\nSend bug reports and comments to danny.bickson@gmail.com\n";
#ifdef GL_NO_MULT_EDGES
  logstream(LOG_WARNING)<<"Code compiled with GL_NO_MULT_EDGES flag - this mode does not support multiple edges between user and movie in different times\n";
#endif
#ifdef GL_NO_MCMC
  logstream(LOG_WARNING)<<"Code compiled with GL_NO_MCMC flag - this mode does not support MCMC methods.\n";
#endif
#ifdef GL_SVD_PP
  logstream(LOG_WARNING)<<"Code compiled with GL_SVD_PP flag - this mode only supports SVD++ run.\n";
#endif

   start(argc, argv);
}
#include <graphlab/macros_undef.hpp>
