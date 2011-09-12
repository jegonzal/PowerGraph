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
#include "clustering.h"
#include "kmeans.hpp"
#include "unittest.hpp"
#include "io.hpp"
#include "../gabp/advanced_config.h"
#include "lda.h"
#ifdef OMP_SUPPORT
#include "omp.h"
#endif

#include <graphlab/macros_def.hpp>
/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU
 See documentation in header file pmf.h

*/

using namespace graphlab;
using namespace itpp;
using namespace std;

advanced_config ac;
problem_setup ps;

const char * runmodesname[] = {"K-Means", "K-Means++", "Fuzzy K-Means", "Latent Dirichlet Allocation"};
const char * inittypenames[]= {"RANDOM", "ROUND_ROBIN", "KMEANS++", "RANDOM_CLUSTER"};
const char * countername[] = {"DISTANCE_CALCULTION", "LDA_NEWTON_METHOD", "LDA_ACCUM_BETA", "LDA_LIKELIHOOD", "LDA_NORMALIZE"};


/* Function declerations */ 
void load_graph(const char* filename, graph_type * g,gl_types::core & glcore);    
void last_iter();
void initialize_clusters();
void dumpcluster();
void tfidf_weighting();
void plus_mul(vec& v1, sparse_vec &v2, double factor);


  void vertex_data::save(graphlab::oarchive& archive) const {  
    ////TODO archive << pvec;
  }  
   
  void vertex_data::load(graphlab::iarchive& archive) {  
     //TODO archive >> pvec;
  }

double get_distance(const vertex_data & v){
    return v.min_distance;
 }

double calc_cost(){
     double cost = 0;
     for (int i=0; i< ps.M; i++){
       const vertex_data & data = ps.g->vertex_data(i);
       cost += data.min_distance;
     }
   return cost;

}




 

int calc_cluster_centers(){
   int total = 0;
   if (ps.algorithm == K_MEANS){
#ifdef OMP_SUPPORT
#pragma omp parallel for
#endif    
     for (int i=0; i< ps.K; i++){
         ps.clusts.cluster_vec[i].location = ps.clusts.cluster_vec[i].cur_sum_of_points / ps.clusts.cluster_vec[i].num_assigned_points;
         ps.clusts.cluster_vec[i].sum_sqr = sum_sqr(ps.clusts.cluster_vec[i].location);
     }
     for (int i=0; i< ps.K; i++)
         total += ps.clusts.cluster_vec[i].num_assigned_points;
   }
  else if (ps.algorithm == K_MEANS_FUZZY){
#ifdef OMP_SUPPORT
#pragma omp parallel for
#endif
    for (int i=0; i< ps.K; i++){
       if (ps.iiter > 0){
         ps.clusts.cluster_vec[i].location = zeros(ps.N);
         ps.clusts.cluster_vec[i].num_assigned_points = ps.M;
       //ps.clusts.cluster_vec[i].cur_sum_of_points = zeros(ps.K);
         double sum_u_i_j = 0;
         for (int j=0; j< ps.M; j++){
            vertex_data & data = ps.g->vertex_data(j);
	    plus_mul(ps.clusts.cluster_vec[i].location, data.datapoint, data.distances[i]);
            sum_u_i_j += data.distances[i];
         }
         assert(sum_u_i_j > 0); 
         ps.clusts.cluster_vec[i].location /= sum_u_i_j;
       }
       ps.clusts.cluster_vec[i].sum_sqr = sum_sqr(ps.clusts.cluster_vec[i].location);
   }
 }
  return total;

}

void add_tasks(){

  std::vector<vertex_id_t> um;
  for (int i=0; i< ps.M; i++)
    um.push_back(i);
 
  switch (ps.algorithm){
     case K_MEANS:
     case K_MEANS_PLUS_PLUS:
     case K_MEANS_FUZZY:
       ps.glcore->add_tasks(um, kmeans_update_function, 1);
       break;
;
     case LDA:
       ps.glcore->add_tasks(um, lda_em_update_function,1);  
       break;
 }

}



void init_clusters(){
  assert(ps.N >0);
  
   for (int i=0; i< ac.K; i++){
        cluster a ;
        a.location = zeros(ps.N); 
        a.cur_sum_of_points = zeros(ps.N);
        ps.clusts.cluster_vec.push_back(a);
   }


}
	

void init_random_cluster(){
   for (int i=0; i < ac.K; i++){
       int tries = 0;
       while(true){
        ::plus(ps.clusts.cluster_vec[i].location,  ps.g->vertex_data(randi(0, ps.M-1)).datapoint);
  	 if (sum(abs(ps.clusts.cluster_vec[i].location))>0)	
           break;
         tries++;
	 if (tries > 100){
	    logstream(LOG_ERROR)<<"Failed to assign non-zero cluster head"<<std::endl;
	    exit(1);
	 }
       }
   }
}

void init(){
#ifdef OMP_SUPPORT
   omp_set_num_threads(ac.omp_support ? ac.ncpus: 1);
   logstream(LOG_INFO) << "setting the number of omp threads to: " << (ac.omp_support ? ac.ncpus: 1) << std::endl;
#endif 

  switch(ps.algorithm){
   case K_MEANS:
   case K_MEANS_FUZZY:
    init_clusters();
    break;

   case K_MEANS_PLUS_PLUS:
   case LDA:
    break;

  }
}


void run_graphlab(){
   for (int i=0; i< ac.iter; i++){
    ps.glcore->start();
    last_iter();
    add_tasks();
  }
}



/** 
 * ==== SETUP AND START
 */
void start(int argc, const char * argv[]) {
  
  ps.gt.start();
  command_line_options clopts;
  ac.init_command_line_options(clopts);
  gl_types::core glcore;
  if (ps.glcore == NULL)
    ps.glcore = &glcore;

  if (ac.mainfunc){ //if called from main(), parse command line arguments
    assert(clopts.parse(argc, argv));
    ac.ncpus = clopts.get_ncpus();

   if (ac.unittest > 0)
      unit_testing(ac.unittest,clopts);
  }
  
  ps.algorithm = (runmodes)ac.algorithm;
  printf("Setting run mode %s\n", runmodesname[ps.algorithm]);


  ps.verify_setup();
  
  ps.glcore->set_engine_options(clopts); 

  logger(LOG_INFO, "%s starting\n",runmodesname[ps.algorithm]);
  //read the training data
  printf("loading data file %s\n", ac.datafile.c_str());
  if (!ac.manualgraphsetup){
  if (!ac.loadgraph){
    ps.g=&ps.glcore->graph();
    load_graph(ac.datafile.c_str(), ps.g,* ps.glcore);

    if (ps.init_type == INIT_RANDOM_CLUSTER)
       init_random_cluster();

    if (ac.savegraph){
	printf("Saving .graph files\n");
	char filename[256];
        sprintf(filename, "%s%d.graph", ac.datafile.c_str(), ac.D);
        std::ofstream fout(filename, std::fstream::binary);
        graphlab::oarchive oarc(fout);
	oarc << ps.M << ps.N << ps.K << ps.L;
        oarc << *ps.g;
        printf("Done!\n");
        fout.close();
	exit(0);
    }

  } else {
    char filename[256];
    sprintf(filename, "%s%d.graph", ac.datafile.c_str(), ac.D);
    std::ifstream fin(filename, std::fstream::binary);
    graphlab::iarchive iarc(fin);
    iarc >> ps.M >> ps.N >> ps.K >> ps.L;
    printf("Loading graph from file\n");
    iarc >> ps.glcore->graph();
    ps.g=&ps.glcore->graph();
    printf("Matrix size is: ROWS %dx COLS %dx CLUSTERS %d ", ps.M, ps.N, ps.K);   
    printf("Creating %d nnz data entires...\n", ps.L);
  }
  }

  if (ac.loadfactors){
     import_from_file();
  }


  if (ac.tfidf)
    tfidf_weighting();

  if (ac.stats){
    calc_stats();
    exit(0);
  }
  

  add_tasks();
 
 
  printf("%s for (%d, %d, %d):%d.\n", runmodesname[ac.algorithm], ps.M, ps.N, ps.K, ps.L);
  

  
  ps.iiter--; 
 
  ps.g->finalize();  

  /**** START GRAPHLAB AND RUN UNTIL COMPLETION *****/
    switch(ps.algorithm){
      case K_MEANS:
      case K_MEANS_FUZZY:
         calc_cluster_centers();
         last_iter();
         run_graphlab();
         break;

      case K_MEANS_PLUS_PLUS:
         initialize_clusters();
     	 break;

      case LDA:
        lda_main();
	break;
  }


  if (ac.clusterdump)
     dumpcluster();

  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (ps.counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], ps.counter[i]); 
  }

  //write output matrices U,V,T to file
  if (ac.binaryoutput)
     export_to_binary_file();
  else if (ac.matrixmarket)
     export_to_matrixmarket();
  else // it++ output
   export_to_itpp_file();
}



void do_main(int argc, const char *argv[]){
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logstream(LOG_INFO)<< "Clustering Code (K-Means/Fuzzy K-means/K-Means++/LDA) written By Danny Bickson, CMU\nSend bug reports and comments to danny.bickson@gmail.com\n";

#ifdef OMP_SUPPORT
  logstream(LOG_INFO)<<"Program compiled with OMP support\n";
#endif
   start(argc, argv);
}


#include <graphlab/macros_undef.hpp>
