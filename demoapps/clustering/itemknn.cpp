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



#include "clustering.h"
#include "../gabp/advanced_config.h"
#include "distance.h"

extern advanced_config ac;
extern problem_setup ps;
extern const char * runmodesname[];
double sum_sqr(sparse_vec & v);

vec wrap_answer(const vec& distances, const ivec& indices, int num){
   vec ret = zeros(num*2);
   for (int i=0; i< num; i++){
      ret[2*i] = distances[i];
      ret[2*i+1] = indices[i];
   }
   return ret;
}

void init_knn(){

   if (ac.distance_measure != EUCLIDEAN)
     return;

   graph_type * training = ps.g<graph_type>(TRAINING);
   for (int i=0; i< ps.M; i++){
      vertex_data & data = training->vertex_data(i);
      data.min_distance = sum_sqr(data.datapoint);
   }

}

 /***
 * UPDATE FUNCTION
 */
void knn_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& vdata = scope.vertex_data();

  int id = scope.vertex();
  bool toprint = ac.debug /*&& (id == 0 || (id == ps.M-1))*/; 
  
 /* print statistics */
  if (toprint){
    printf("Item Knn: entering data point %u, current cluster %d\n",  id, vdata.current_cluster);  
    print(vdata.datapoint); 
  }

  if (!vdata.reported) //this matrix row have no non-zero entries, and thus ignored
     return;

  graphlab::timer t; t.start();

  vec distances = zeros(ps.M);
  for (int i=0; i<ps.M; i++){
      vertex_data & other = ps.g<graph_type>(TRAINING)->vertex_data(i);
      distances[i] = calc_distance(vdata.datapoint, other.datapoint, other.min_distance);
  }

  ivec indices = sort_index(distances);
  sort(distances); 
  vdata.distances = wrap_answer(distances, indices, ps.K);
  if (toprint)
    printf("Closest is: %d with distance %g\n", (int)vdata.distances[1], vdata.distances[0]);
}

void copy_assignments(mat &a, const vec& distances, int i){
   for (int j=0; j<ps.K; j++){
     set_val(a,j, i, distances[j*2+1]);
   }
}

void copy_distances(mat &a, const vec& distances, int i){
   for (int j=0; j<ps.K; j++){
     set_val(a,j, i, distances[j*2]);
   }
}


void prepare_output(){

  graph_type * validation = ps.g<graph_type>(VALIDATION);
  ps.output_assignements = zeros(ps.K, validation->num_vertices());
  ps.output_clusters = zeros(ps.K, validation->num_vertices());
  for (int i=0; i< (int)validation->num_vertices(); i++){
     const vertex_data& data = validation->vertex_data(i);
     copy_assignments(ps.output_assignements, data.distances, i);
     copy_distances(ps.output_clusters, data.distances, i); 
  }
}


 
void knn_main(){

    init_knn();
    ps.gt.start();
    ps.glcore->start();

    prepare_output();
};


