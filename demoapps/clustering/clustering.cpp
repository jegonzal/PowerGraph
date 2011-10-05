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
using namespace std;

advanced_config ac;
problem_setup ps;

const char * runmodesname[] = {"K-Means", "K-Means++", "Fuzzy K-Means", "Latent Dirichlet Allocation"};
const char * inittypenames[]= {"RANDOM", "ROUND_ROBIN", "KMEANS++", "RANDOM_CLUSTER"};
const char * countername[] = {"DISTANCE_CALCULTION", "LDA_NEWTON_METHOD", "LDA_ACCUM_BETA", "LDA_LIKELIHOOD", "LDA_NORMALIZE"};


/* Function declerations */ 
void load_graph(const char* filename, graph_type * g,gl_types::core & glcore);    
void last_iter();
void initialize_clusters(gl_types::core & glcore);
void initialize_clusters(gl_types_kcores::core & glcore){ assert(false); }
void dumpcluster();
void tfidf_weighting();
void plus_mul(vec& v1, sparse_vec &v2, double factor);
void kcores_update_function(gl_types_kcores::iscope & scope, gl_types_kcores::icallback & scheduler);
void kcores_main();

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
       const vertex_data & data = ps.g<graph_type>()->vertex_data(i);
       cost += data.min_distance;
     }
     ps.cost = cost;
   return cost;

}




 

int calc_cluster_centers(){
   int total = 0;
   if (ps.algorithm == K_MEANS){
#ifdef OMP_SUPPORT
#pragma omp parallel for
#endif    
     for (int i=0; i< ps.K; i++){
         if (ps.clusts.cluster_vec[i].num_assigned_points == 0){
	   if (ps.iiter == 0)
	       assert(false);
          }
          else {
              ps.clusts.cluster_vec[i].location = ps.clusts.cluster_vec[i].cur_sum_of_points / ps.clusts.cluster_vec[i].num_assigned_points;
              ps.clusts.cluster_vec[i].sum_sqr = sum_sqr(ps.clusts.cluster_vec[i].location);
          }
          if (ac.debug)
            std::cout<<"New cluster center is: " << ps.clusts.cluster_vec[i].location << std::endl;
     }
     for (int i=0; i< ps.K; i++)
         total += ps.clusts.cluster_vec[i].num_assigned_points;
   }
  else if (ps.algorithm == K_MEANS_FUZZY){
#ifdef OMP_SUPPORT
#pragma omp parallel for
#endif
   /* see algo description in http://www.cs.princeton.edu/courses/archive/fall08/cos436/Duda/C/fk_means.htm . The below code replaces m_i with fuzzy mean of all examples in the cluster */
   /* for each cluster 1...k */
    for (int i=0; i< ps.K; i++){
      //if (ps.iiter > 0){
         ps.clusts.cluster_vec[i].location = zeros(ps.N);
         ps.clusts.cluster_vec[i].num_assigned_points = ps.M;
         double sum_u_i_j = 0;
         /* for each pot 1..M*/
         for (int j=0; j< ps.M; j++){
            vertex_data & data = ps.g<graph_type>()->vertex_data(j);
            /* m_i = \sum u(j,i)^2 x_j */
	    plus_mul(ps.clusts.cluster_vec[i].location, data.datapoint, powf(data.distances[i], 1));
            sum_u_i_j += powf(data.distances[i],1);
         }
         assert(sum_u_i_j > 0); 
         ps.clusts.cluster_vec[i].location /= sum_u_i_j;
         if (ac.debug)
           std::cout<<" cluster " << i << " is now on: " << ps.clusts.cluster_vec[i].location << std::endl;
     //}    
     ps.clusts.cluster_vec[i].sum_sqr = sum_sqr(ps.clusts.cluster_vec[i].location);
    }
  }

  if (ac.algorithm == K_MEANS && ac.init_mode != INIT_KMEANS_PLUS_PLUS)
      assert(total == ps.total_assigned);
  return total;

}

void add_tasks(gl_types::core & glcore){

  std::vector<vertex_id_t> um;
  for (int i=0; i< ps.M; i++)
    um.push_back(i);
 
  switch (ps.algorithm){
     case K_MEANS:
     case K_MEANS_PLUS_PLUS:
     case K_MEANS_FUZZY:
       glcore.add_tasks(um, kmeans_update_function, 1);
       break;
;
     case LDA:
       glcore.add_tasks(um, lda_em_update_function,1);  
       break;

     default: assert(false);
  }
}

void add_tasks(gl_types_kcores::core & glcore){

  std::vector<vertex_id_t> um;
  for (int i=0; i< ps.M; i++)
    um.push_back(i);
 
  switch (ps.algorithm){
     case KSHELL_DECOMPOSITION:
       glcore.add_tasks(um, kcores_update_function, 1);
       break;

     default: assert(false);
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
        ::plus(ps.clusts.cluster_vec[i].location,  ps.g<graph_type>()->vertex_data(randi(0, ps.M-1)).datapoint);
  	 if (sum(abs(ps.clusts.cluster_vec[i].location))>0){	
            if (ac.debug)
              std::cout<<"Selected random cluster: " << i << " to be: " << ps.clusts.cluster_vec[i].location << std::endl;
            break;
         }
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


   case KSHELL_DECOMPOSITION:
     break; 
  }
}


void run_graphlab(gl_types::core & glcore){
   for (int i=0; i< ac.iter; i++){
    glcore.start();
    last_iter();
    add_tasks(glcore);
  }
}

void run_graphlab(gl_types_kcores::core & glcore){
  assert(false);
}




/** 
 * ==== SETUP AND START
 */
template<typename core, typename graph_type>
void start(command_line_options & clopts) {
  
  core glcore;
  glcore.set_engine_options(clopts); 
  ps.set_core(&glcore);
  ps.verify_setup();
  

  logger(LOG_INFO, "%s starting\n",runmodesname[ps.algorithm]);
  //read the training data
  printf("loading data file %s\n", ac.datafile.c_str());
  if (!ac.manualgraphsetup){
  if (!ac.loadgraph){
    ps.set_graph(&glcore.graph());
    load_graph<graph_type>(ac.datafile.c_str(), &glcore.graph());

    if (ps.init_type == INIT_RANDOM_CLUSTER)
       init_random_cluster();

    if (ac.savegraph){
	printf("Saving .graph files\n");
	char filename[256];
        sprintf(filename, "%s%d.graph", ac.datafile.c_str(), ac.D);
        std::ofstream fout(filename, std::fstream::binary);
        graphlab::oarchive oarc(fout);
	oarc << ps.M << ps.N << ps.K << ps.L;
        oarc << *ps.g<graph_type>();
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
    iarc >> *ps.g<graph_type>();
    printf("Matrix size is: ROWS %dx COLS %dx CLUSTERS %d ", ps.M, ps.N, ps.K);   
    printf("Creating %d nnz data entires...\n", ps.L);
  }
  }

  if (ac.loadfactors){
     //import_from_file();
  }


  if (ac.tfidf)
    tfidf_weighting();

  if (ac.stats){
    calc_stats();
    exit(0);
  }
  

  add_tasks(glcore);
 
 
  printf("%s for (%d, %d, %d):%d.\n", runmodesname[ac.algorithm], ps.M, ps.N, ps.K, ps.L);
  

  
  ps.iiter--; 
 
  ps.g<graph_type>()->finalize();  

  /**** START GRAPHLAB AND RUN UNTIL COMPLETION *****/
    switch(ps.algorithm){
      case K_MEANS:
      case K_MEANS_FUZZY:
         if (ps.algorithm == K_MEANS_FUZZY){
            init_fuzzy_kmeans();
         }
         else 
           calc_cluster_centers();
         last_iter();
         run_graphlab(glcore);
         break;

      case K_MEANS_PLUS_PLUS:
         initialize_clusters(glcore);
     	 break;

      case LDA:
        lda_main();
	break;
  
      case KSHELL_DECOMPOSITION:
        kcores_main();
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
     export_to_binary_file<graph_type>();
  else if (ac.matrixmarket)
     export_to_matrixmarket<graph_type>();
  else // it++ output
   export_to_itpp_file<graph_type>();
}



int do_main(int argc, const char *argv[]){
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logstream(LOG_INFO)<< "Clustering Code (K-Means/Fuzzy K-means/K-Means++/LDA) written By Danny Bickson, CMU\nSend bug reports and comments to danny.bickson@gmail.com\n";

#ifdef OMP_SUPPORT
  logstream(LOG_INFO)<<"Program compiled with OMP support\n";
#endif

  ps.gt.start();
  command_line_options clopts;
  ac.init_command_line_options(clopts);

  if (ac.mainfunc){ //if called from main(), parse command line arguments
    if (!clopts.parse(argc, argv))
       return EXIT_FAILURE;
    ac.ncpus = clopts.get_ncpus();

   if (ac.unittest > 0)
      unit_testing(ac.unittest,clopts);
  }
  
  ps.algorithm = (runmodes)ac.algorithm;
  printf("Setting run mode %s\n", runmodesname[ps.algorithm]);
    switch(ps.algorithm){
    case K_MEANS:
    case K_MEANS_PLUS_PLUS:
    case K_MEANS_FUZZY:
    case LDA: 
       start<gl_types::core, graph_type>(clopts);
       break;

    case KSHELL_DECOMPOSITION:
       start<gl_types_kcores::core, graph_type_kcores>(clopts);
       break;
  }  


  return EXIT_SUCCESS;
}


#include <graphlab/macros_undef.hpp>
