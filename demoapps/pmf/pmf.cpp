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
#include <stdlib.h>

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
#include "../gabp/advanced_config.h"

#ifdef GL_SVD_PP
#include "svdpp.hpp"
#endif

#include <graphlab/macros_def.hpp>

const char * runmodesname[] = {"ALS_MATRIX (Alternating least squares)", "BPTF_MATRIX (Bayesian Prob. Matrix Factorization)", "BPTF_TENSOR (Bayesian Prob. Tensor Factorization)", "BPTF_TENSOR_MULT", "ALS_TENSOR_MULT", "SVD++", "SGD (Stochastic Gradient Descent)", "SVD (Singular Value Decomposition via LANCZOS)", "NMF (non-negative factorization)", "Weighted alternating least squares", "Alternating least squares with sparse user factor matrix", "Alternating least squares with doubly sparse (user/movie) factor matrices", "Alternating least squares with sparse movie factor matrix"};

const char * countername[] = {"EDGE_TRAVERSAL", "BPTF_SAMPLE_STEP", "CALC_RMSE_Q", "ALS_LEAST_SQUARES", \
  "BPTF_TIME_EDGES", "BPTF_LEAST_SQUARES", "CALC_OBJ", "BPTF_MVN_RNDEX", "BPTF_LEAST_SQUARES2", "SVD_MULT_A", "SVD_MULT_A_TRANSPOSE"};

const char * testtypename[] = {"TRAINING", "VALIDATION", "TEST"};





using namespace graphlab;
using namespace itpp;
using namespace std;

advanced_config ac;
problem_setup ps;





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

double get_rmse(const vertex_data & v){
    return v.rmse;
 }

#ifndef GL_SVD_PP
void svd_init(){};
void svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler){};
#endif



//methods to compute the Root mean square error (RMSE)     
float predict(const vec& x1, const vec& x2, const edge_data * edge, float rating, float & prediction){
	prediction = dot(x1, x2);	
   //return the squared error
   prediction = std::min((double)prediction, ac.maxval);
   prediction = std::max((double)prediction, ac.minval);
	float sq_err = powf(prediction - rating, 2);
   return sq_err;
}

double predict(const vertex_data& user, const vertex_data &movie, const edge_data * edge, float rating, float & prediction){
#ifdef GL_SVD_PP	
   return svd_predict(user, movie, edge, rating, prediction);
#else
   return predict(user.pvec, movie.pvec, edge, rating, prediction);
#endif
} 


float predict(const vertex_data& v1, const vertex_data& v2, const edge_data * edge, const vertex_data *v3, float rating, float &prediction){
	if (v3 == NULL) //matrix	
		return predict(v1,v2,edge, rating,prediction);

	prediction = 0;
	for (int i=0; i< v1.pvec.size(); i++){
	   prediction += (v1.pvec[i] * v2.pvec[i] * v3->pvec.get(i));
	}
   prediction = std::min((double)prediction, ac.maxval);
   prediction = std::max((double)prediction, ac.minval);
   float sq_err = powf(prediction - rating, 2);
   return sq_err;
   
}

  //constructor
  vertex_data::vertex_data(){
    pvec = zeros(ac.D);
    rmse = 0;
    num_edges = 0;
#ifdef GL_SVD_PP
    bias =0;
    weight = zeros(ac.D);
#endif
  }

  void vertex_data::save(graphlab::oarchive& archive) const {  
    ////TODO archive << pvec;
    archive << rmse << num_edges; 
#ifdef GL_SVD_PP
    archive << bias << weight;
#endif
  }  
   
  void vertex_data::load(graphlab::iarchive& archive) {  
     //TODO archive >> pvec;
     archive >> rmse >> num_edges;  
#ifdef GL_SVD_PP
     archive >> bias >> weight;
#endif
  }


/**
 * printout RMSE statistics after each iteration
 */
void last_iter(){
  printf("Entering last iter with %d\n", ps.iiter);

  double res,res2;
  double rmse = (ps.algorithm != STOCHASTIC_GRADIENT_DESCENT && ps.algorithm != NMF) ? agg_rmse_by_movie(res) : agg_rmse_by_user(res);
  //rmse=0;
  printf(ac.printhighprecision ? 
        "%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n":
        "%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n"
        , ps.gt.current_time(), runmodesname[ps.algorithm], ps.iiter,calc_obj(res),  rmse, calc_rmse_wrapper(&ps.validation_graph, true, res2));
  ps.iiter++;

  if (ps.BPTF)
    last_iter_bptf(res);        
}


void add_tasks(gl_types::core & glcore){

  std::vector<vertex_id_t> um;
  for (int i=0; i< ps.M+ps.N; i++)
    um.push_back(i);
 
 // add update function for user and movie nodes (tensor dims 1+2) 
  switch (ps.algorithm){
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
  if (ps.tensor){
    std::vector<vertex_id_t> tv;
    for (int i=ps.M+ps.N; i< ps.M+ps.N+ps.K; i++)
      tv.push_back(i);
    glcore.add_tasks(tv, time_node_update_function, 1);
  }


}


void init(){

  if (ps.BPTF)
     init_self_pot(); 

  switch(ps.algorithm){
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
        double validation_rmse = calc_rmse_wrapper(&ps.validation_graph, true, res2);
        printf(ac.printhighprecision ? 
              "Final result. Obj=%g, TRAIN RMSE= %0.12f VALIDATION RMSE= %0.12f.\n":
              "Final result. Obj=%g, TRAIN RMSE= %0.4f VALIDATION RMSE= %0.4f.\n"
              , obj,  train_rmse, validation_rmse);
        double runtime = gt.current_time();
        printf("Finished in %lf seconds\n", runtime);
        if (ac.unittest > 0){
            verify_result(obj, train_rmse, validation_rmse);
        }

}



/** 
 * ==== SETUP AND START
 */
void start(int argc, const char * argv[]) {
   
  command_line_options clopts;
  ac.init_command_line_options(clopts);
  gl_types::core glcore;
  if (ps.glcore == NULL)
    ps.glcore = &glcore;

  if (ac.mainfunc){ //if called from main(), parse command line arguments
    assert(clopts.parse(argc, argv));

   if (ac.unittest > 0)
      unit_testing(ac.unittest,clopts);
  }
  
  ps.algorithm = (runmodes)ac.algorithm;
  printf("Setting run mode %s\n", runmodesname[ps.algorithm]);

  if (ac.scheduler == "round_robin"){
    char schedulerstring[256];
    sprintf(schedulerstring, "round_robin(max_iterations=%d,block_size=1)", ac.iter);
    clopts.set_scheduler_type(schedulerstring);
    assert(ac.iter > 0);
  }
  ps.verify_setup();
  ps.glcore->set_engine_options(clopts); 

  logger(LOG_INFO, "%s starting\n",runmodesname[ps.algorithm]);
  //read the training data
  printf("loading data file %s\n", ac.datafile.c_str());
  if (!ac.manualgraphsetup){
  if (!ac.loadgraph){
    ps.g=&ps.glcore->graph();
    load_pmf_graph(ac.datafile.c_str(), ps.g, TRAINING,* ps.glcore);

  //read the vlidation data (optional)
    printf("loading data file %s\n", (ac.datafile+"e").c_str());
    load_pmf_graph((ac.datafile+"e").c_str(),&ps.validation_graph, VALIDATION, *ps.glcore);

  //read the test data (optional)
    printf("loading data file %s\n", (ac.datafile+"t").c_str());
    load_pmf_graph((ac.datafile+"t").c_str(),&ps.test_graph, TEST, *ps.glcore);


    if (ac.savegraph){
	printf("Saving .graph files\n");
	char filename[256];
        sprintf(filename, "%s%d.graph", ac.datafile.c_str(), ac.D);
        std::ofstream fout(filename, std::fstream::binary);
        graphlab::oarchive oarc(fout);
	oarc << ps.M << ps.N << ps.K << ps.L << ps.Le << ps.Lt << ac.D;
        oarc << *ps.g << ps.validation_graph << ps.test_graph;
        printf("Done!\n");
        fout.close();
	exit(0);
    }

  } else {
    char filename[256];
    sprintf(filename, "%s%d.graph", ac.datafile.c_str(), ac.D);
    std::ifstream fin(filename, std::fstream::binary);
    graphlab::iarchive iarc(fin);
    iarc >> ps.M >> ps.N >> ps.K >> ps.L >> ps.Le >> ps.Lt >> ac.D;
    printf("Loading graph from file\n");
    iarc >> ps.glcore->graph() >> ps.validation_graph >> ps.test_graph;
    ps.g=&ps.glcore->graph();
    printf("Matrix size is: USERS %dx MOVIES %dx TIME BINS %d D=%d\n", ps.M, ps.N, ps.K, ac.D);   
    printf("Creating %d edges (observed ratings)...\n", ps.L);
  }
  }

  if (ac.loadfactors){
     import_uvt_from_file();
  }


  if (ac.stats){
    calc_stats(TRAINING);
    calc_stats(VALIDATION);
    calc_stats(TEST);
    exit(0);
  }

  if (ps.algorithm == ALS_TENSOR_MULT || ps.algorithm == ALS_MATRIX || ps.algorithm == ALS_SPARSE_USR_FACTOR || ps.algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || ps.algorithm == ALS_SPARSE_MOVIE_FACTOR || ps.algorithm == WEIGHTED_ALS){
    printf("setting regularization weight to %g\n", ac.als_lambda);
    pU=pV=ac.als_lambda;
  }

   if (ps.tensor){
    ps.dp = GenDiffMat(ps.K)*ps.pT;
    if (ac.debug)
      std::cout<<ps.dp<<std::endl;
  }
  
  add_tasks(*ps.glcore);

  
  printf("%s for %s (%d, %d, %d):%d.  D=%d\n", runmodesname[ac.algorithm], ps.tensor?"tensor":"matrix", ps.M, ps.N, ps.K, ps.L, ac.D);
  
  init();

   if (ac.datafile == "netflix" || ac.datafile == "netflix-r"){
       ac.minval = 1; ac.maxval = 5;
    }
   else if ((ac.datafile == "kddcup" || ac.datafile == "kddcup2") && ac.maxval == DEF_MAX_VAL){
       ac.minval = 0; ac.maxval = 100;
   }

   if (ps.algorithm != LANCZOS){
     double res, res2;
     double rmse =  calc_rmse_wrapper(ps.g, false, res);
     printf(ac.printhighprecision ? 
           "complete. Objective=%g, TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n" :
           "complete. Objective=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n" 
           , calc_obj(res), rmse, calc_rmse(&ps.validation_graph, true, res2));
  }
 
 
  if (ps.BPTF){
    //sample hyper priors and noise level
    if (ac.bptf_delay_alpha < ps.iiter)
    	sample_alpha(ps.L);
    sample_U();
    sample_V();
    if (ps.tensor) 
      sample_T();
  }

  ps.g->finalize();  
  ps.gt.start();

  /**** START GRAPHLAB AND RUN UNTIL COMPLETION *****/
    switch(ps.algorithm){
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
         run_graphlab(*ps.glcore, ps.gt);
         break;
     
     case LANCZOS:
        lanczos(*ps.glcore); break;

     case NMF:
        nmf(ps.glcore); break;
  }

 if (ps.algorithm != LANCZOS){
    /**** OUTPUT KDD FORMAT *****/
    if (ac.datafile == "kddcup" || ac.datafile == "kddcup2"){
      if (ac.outputvalidation) //experimental: output prediction of validation data
 	     export_kdd_format(&ps.validation_graph, VALIDATION, true);
      else //output prediction of test data, as required by KDD 
	     export_kdd_format(&ps.test_graph, TEST, true);
    }
 }

  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (ps.counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], ps.counter[i]); 
  }

  //write output matrices U,V,T to file
  if (ac.binaryoutput)
     export_uvt_to_binary_file();
  else if (ac.matrixmarket)
     export_uvt_to_matrixmarket();
  else // it++ output
   export_uvt_to_itpp_file();
}



void do_main(int argc, const char *argv[]){
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
