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
//double sum_sqr(sparse_vec & v);
const int matlab_offset = 1;

flt_dbl_vec wrap_answer(const flt_dbl_vec& distances, const ivec& indices, int num){
   flt_dbl_vec ret = zeros(num*2);
   for (int i=0; i< num; i++){
      ret[2*i] = distances[i];
      ret[2*i+1] = indices[i];
   }
   return ret;
}

void init_knn(){

   if (ac.distance_measure != EUCLIDEAN && ac.distance_measure != COSINE && ac.distance_measure != TANIMOTO)
     return;

   int start = (ps.algorithm == USER_KNN ? 0 : ps.M);
   int end = (ps.algorithm == USER_KNN ? ps.M : ps.M+ps.N);
  
   graph_type * training = ps.g<graph_type>(TRAINING);
   for (int i=start; i< end; i++){
      vertex_data & data = training->vertex_data(i);
      data.min_distance = sum_sqr(data.datapoint);
      assert( data.reported == ( nnz(data.datapoint) > 0));
   }

   int startv = (ps.algorithm == USER_KNN ? 0 : ps.M_validation);
   int endv = (ps.algorithm == USER_KNN ? ps.M_validation : ps.M_validation+ps.N_validation);

   graph_type * validation = ps.g<graph_type>(VALIDATION);
   for (int i=startv; i< endv; i++){
      vertex_data & data = validation->vertex_data(i);
      data.min_distance = sum_sqr(data.datapoint);
   }
}


void stats(){

   flt_dbl min=1e100, max=0, avg=0;
   int cnt = 0;
   graph_type * validation = ps.g<graph_type>(VALIDATION);
   int startv = (ps.algorithm == USER_KNN ? 0 : ps.M_validation);
   int endv = (ps.algorithm == USER_KNN ? ps.M_validation : ps.M_validation+ps.N_validation);

   for (int i=startv; i< endv; i++){
     vertex_data& data = validation->vertex_data(i);
     if (data.distances.size() > 0){
       min = std::min(min, data.distances[0]);
       max = std::max(max, data.distances[0]);
       if (std::isnan(data.distances[0]))
          printf("bug: nan on %d\n", i);
       else {
         avg += data.distances[0];    
         cnt++;
       }
     }
   }

  printf("Distance statistics: min %g max %g avg %g\n", min, max, avg/cnt);
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
  graph_type *training = ps.g<graph_type>(TRAINING);

   int start = (ps.algorithm == USER_KNN ? 0 : ps.M);
   int end = (ps.algorithm == USER_KNN ? ps.M : ps.M+ps.N);
  int howmany = (end-start)*ac.knn_sample_percent;
  flt_dbl_vec distances = zeros(howmany);
  ivec indices = ivec(howmany);
   if (ac.knn_sample_percent == 1.0){
     for (int i=start; i< end; i++){
        vertex_data & other = training->vertex_data(i);
        distances[i-start] = calc_distance(vdata.datapoint, other.datapoint, other.min_distance, vdata.min_distance);
        indices[i-start] = i;
     }
  }
  else for (int i=0; i<howmany; i++){
        int random_other = ::randi(start, end-1);
        while(random_other == id){
	  random_other == ::randi(start, end-1);
        }
        vertex_data & other = training->vertex_data(random_other);
        distances[i] = calc_distance(vdata.datapoint, other.datapoint, other.min_distance, vdata.min_distance);
        indices[i] = random_other;
  }

  ivec indices_sorted = sort_index2(distances, indices);
  sort(distances); 
  vdata.distances = wrap_answer(distances, indices_sorted, ps.K);
  if (toprint)
    printf("Closest is: %d with distance %g\n", (int)vdata.distances[1], vdata.distances[0]);


  if (id % 100 == 0)
    printf("handling validation row %d\n", id);
}

void copy_assignments(flt_dbl_mat &a, const flt_dbl_vec& distances, int i, graph_type* validation){
   bool reported = validation->vertex_data(i).reported;
   for (int j=0; j<ac.K; j++){
     set_val(a,j, i, reported? distances[j*2+1]+matlab_offset : -1);
   }
}

void copy_distances(flt_dbl_mat &a, const flt_dbl_vec& distances, int i, graph_type* validation){
   bool reported = validation->vertex_data(i).reported;
   for (int j=0; j<ac.K; j++){
     set_val(a,j, i, reported? distances[j*2] : -1);
   }
}


void prepare_output(){

  graph_type * validation = ps.g<graph_type>(VALIDATION);
  ps.output_assignements = zeros(ps.K, validation->num_vertices());
  ps.output_clusters = zeros(ps.K, validation->num_vertices());
  for (int i=0; i< (int)validation->num_vertices(); i++){
     const vertex_data& data = validation->vertex_data(i);
     copy_assignments(ps.output_assignements, data.distances, i, validation);
     copy_distances(ps.output_clusters, data.distances, i, validation); 
  }
  ps.output_clusters = transpose(ps.output_clusters);
  ps.output_assignements = transpose(ps.output_assignements);
}


 
void knn_main(){

    init_knn();
    ps.gt.start();
    ps.glcore->start();
    stats();
    prepare_output();
};


