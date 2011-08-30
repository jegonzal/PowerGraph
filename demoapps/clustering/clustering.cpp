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



/* Function declerations */ 
void load_graph(const char* filename, graph_type * g,gl_types::core & glcore);    
void last_iter();


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
 

/**
 * printout cost after each iteration
 */
void last_iter(){
  printf("Entering last iter with %d\n", ps.iiter);

  double cost = calc_cost();
  printf("%g) Iter %s %d  Cost=%g\n",
        ps.gt.current_time(), 
	runmodesname[ps.algorithm], 
	ps.iiter,
        cost);
  ps.iiter++;

}


void add_tasks(gl_types::core & glcore){

  std::vector<vertex_id_t> um;
  for (int i=0; i< ps.M; i++)
    um.push_back(i);
 
  switch (ps.algorithm){
     case K_MEANS:
       glcore.add_tasks(um, kmeans_update_function, 1);
       break;
 }

}


void init(){


  switch(ps.algorithm){
   case K_MEANS:
     //TODO; 
     break;

  }
}


void run_graphlab(gl_types::core &glcore,timer & gt ){
    glcore.start();
    last_iter();
    ps.iiter--;
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
    load_graph(ac.datafile.c_str(), ps.g,* ps.glcore);

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


  if (ac.stats){
    calc_stats();
    exit(0);
  }

  add_tasks(*ps.glcore);

  
  printf("%s for (%d, %d, %d):%d.\n", runmodesname[ac.algorithm], ps.M, ps.N, ps.K, ps.L);
  
  init();


  last_iter();
  ps.iiter--; 
 
  ps.g->finalize();  
  ps.gt.start();

  /**** START GRAPHLAB AND RUN UNTIL COMPLETION *****/
    switch(ps.algorithm){
      case K_MEANS:
         run_graphlab(*ps.glcore, ps.gt);
         break;
     
  }

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
  logstream(LOG_INFO)<< "Clustering Code (K-Means) written By Danny Bickson, CMU\nSend bug reports and comments to danny.bickson@gmail.com\n";

   start(argc, argv);
}


#include <graphlab/macros_undef.hpp>
