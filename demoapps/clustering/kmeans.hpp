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
void plus( vec &v1,  sparse_vec &v2);
void minus( vec &v1, sparse_vec &v2);
int calc_cluster_centers();


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

  if (!vdata.reported) //this matrix row have no non-zero entries, and thus ignored
     return;
  //const mat &clusters = CLUSTER_LOCATIONS.get_val();

  double min_dist = 1e100;
  int pos = -1;
  graphlab::timer t; t.start();
  for (int i=0; i< ps.K; i++){
     //vec row = clusters.get_row(i);
     vec & row = ps.clusts.cluster_vec[i].location;
     double dist = calc_distance(vdata.datapoint, row, ps.clusts.cluster_vec[i].sum_sqr);
     if (min_dist > dist){
         min_dist = dist;
         pos = i;
     }
  }  
  ps.counter[DISTANCE_CALCULATION] += t.current_time();
  assert(pos != -1);


  if (pos != vdata.current_cluster){
    vdata.hot=true;
    if (toprint)
      std::cout <<id<<" assigned to cluster: " << pos << std::endl;
  }
  vdata.prev_cluster = vdata.current_cluster; 
  vdata.current_cluster = pos;
  vdata.min_distance = min_dist;

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
   for (int i=0; i< ps.M; i++){
       vertex_data & data = ps.g->vertex_data(i);
       if (!data.hot)
         continue;
       else {
         assert(data.current_cluster >=0 && data.current_cluster < ps.K);
         assert(data.prev_cluster >=0 && data.prev_cluster < ps.K);
         assert(data.prev_cluster != data.current_cluster);
         plus(ps.clusts.cluster_vec[data.current_cluster].cur_sum_of_points , data.datapoint);    
         minus(ps.clusts.cluster_vec[data.prev_cluster].cur_sum_of_points , data.datapoint);    
         ps.clusts.cluster_vec[data.current_cluster].num_assigned_points++;
         ps.clusts.cluster_vec[data.prev_cluster].num_assigned_points--;
         data.hot = false;  
         if (ac.debug)
           std::cout<<"in hot node: " << i << std::endl;
       }
   }
   int total_assigned =calc_cluster_centers();
   assert(total_assigned == ps.total_assigned);

}
 


#include <graphlab/macros_undef.hpp>
#endif //_ALS_HPP
