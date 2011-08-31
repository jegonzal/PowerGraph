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
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 */


#ifndef _KMEANS_HPP
#define _KMEANS_HPP

#include "clustering.h"
#include "../gabp/advanced_config.h"
#include "distance.h"
#include <graphlab/macros_def.hpp>

extern advanced_config ac;
extern problem_setup ps;


void last_iter();
double calc_cost();
void update_clusters();
vec plus(vec & v1, const sparse_vec & v2);

 /***
 * UPDATE FUNCTION
 */
void kmeans_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& vdata = scope.vertex_data();

  int id = scope.vertex();
  bool toprint = ac.debug /*&& (id == 0 || (id == ps.M-1))*/; 
  
 /* print statistics */
  if (toprint){
    printf("entering data point %u, current cluster %d\n",  id, vdata.current_cluster);   
  }

  const mat &clusters = CLUSTER_LOCATIONS.get_val();

  double min_dist = 1e100;
  int pos = -1;
  graphlab::timer t; t.start();
  for (int i=0; i< ps.K; i++){
     vec row = clusters.get_row(i);
     double dist = calc_distance(vdata.datapoint, row);
     if (min_dist > dist){
         min_dist = dist;
         pos = i;
     }
  }  
  ps.counter[DISTANCE_CALCULATION] += t.current_time();
  assert(pos != -1);


  if (toprint && pos != vdata.current_cluster){
    std::cout <<id<<" assigned to cluster: " << pos << std::endl;
  }
      
  vdata.current_cluster = pos;
  vdata.min_distance = min_dist;

  if (id == ps.M-1)
	last_iter();
}


/**
 * printout cost after each iteration
 */
void last_iter(){
  printf("Entering last iter with %d\n", ps.iiter);

  update_clusters();
  double cost = calc_cost();
  printf("%g) Iter %s %d  Cost=%g\n",
        ps.gt.current_time(), 
	runmodesname[ps.algorithm], 
	ps.iiter,
        cost);
  ps.iiter++;

}



void update_clusters(){
 

   mat means = zeros(ps.K, ps.N);
   vec count_in_cluster = zeros(ps.K);
   for (int i=0; i< ps.M; i++){
       const vertex_data & data = ps.g->vertex_data(i);
       assert(data.current_cluster >=0 && data.current_cluster < ps.K);
       count_in_cluster[data.current_cluster]++;
       vec row = means.get_row(data.current_cluster);
       row = plus(row, data.datapoint);
       means.set_row(data.current_cluster, row);
   }
   for (int i=0; i< ps.K; i++){
      vec row = means.get_row(i);
      if (count_in_cluster[i] > 0)
        row /= count_in_cluster[i];
      means.set_row(i, row);
   }

   CLUSTER_LOCATIONS.set(means);

   if (ac.debug){
      std::cout<<"counter in clusters: " << count_in_cluster<< std::endl;
      std::cout<<"cluster locations: " << means << std::endl;
   }
}
 


#include <graphlab/macros_undef.hpp>
#endif //_ALS_HPP
