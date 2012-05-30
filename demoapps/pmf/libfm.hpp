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
 */

/*
 * This file contains an implementation of the mehod: Steffen Rendle (2010): Factorization Machines, in Proceedings of the 10th IEEE International Conference on Data Mining (ICDM 2010), Sydney, Australia.
 * Original implementation by Qiang Yan, Chinese Academy of Science.
 * Adapted to GraphLab by Danny Bickson, CMU
 *
 */

#ifndef __LIBFM_HPP
#define __LIBFM_HPP

#include "graphlab.hpp"
#include "als.hpp"
#include "svdpp.hpp"
#include "stats.hpp"
#include <graphlab/macros_def.hpp>

double reg0 = 0.1;
/**
 *
 *  Implementation of the bias SVD algorithm
 *
 * */
extern advanced_config ac;
extern problem_setup ps;

using namespace graphlab;

std::vector<vertex_data> last_items;

struct vertex_data_libfm{
   float * bias;
   double * v;
   int *last_item; 
   float * rmse;

   vertex_data_libfm(const vertex_data & vdata){
     v = (double*)&vdata.pvec[0];
     bias = (float*)&vdata.bias;
     last_item = (int*)&vdata.num_edges;
     rmse = (float*)& vdata.rmse;
   }
   
   vertex_data_libfm & operator=(vertex_data & data){
     v = (double*)&data.pvec[0];
     bias = (float*)&data.bias;
     last_item = (int*)&data.num_edges;
     rmse = (float*)&data.rmse;
		 return * this;
    }   
};



void libfm_update_function(gl_types_mcmc::iscope &scope, 
			 gl_types_mcmc::icallback &scheduler) {
   assert(false);
} 
void libfm_update_function(gl_types_mult_edge::iscope &scope, 
			 gl_types_mult_edge::icallback &scheduler) {
   assert(false);
}
template<typename graph_type>
void init_libfm(graph_type* _g){
  assert(false);
}

template<>
void init_libfm(graph_type *_g){
   fprintf(stderr, "%s %d factors\n",runmodesname[ps.algorithm], ac.D);
   double factor = 0.1/sqrt(ac.D);
#pragma omp parallel for
   for (int i=0; i<ps.M+ps.N; i++){
       vertex_data & vdata = _g->vertex_data(i);
       vdata.pvec = zeros(ac.D);
       vdata.bias = 0;
       vertex_data_libfm data(_g->vertex_data(i));
       for (int j=0; j < ac.D; j++){
          data.v[j] = (ac.debug ? 0.1 : (::randu()*factor));
       }
       if (i >= ps.M){
         vertex_data data;
         data.pvec = (ac.debug ? ones(ac.D)*0.1 : ::randu(ac.D) * factor);
         last_items.push_back(data);
       }
       else { //user node. find the last rated item and store it
         gl_types::edge_list outs = _g->out_edge_ids(i);
         int max_time = 0;
         vdata.num_edges = -1;
         foreach(graphlab::edge_id_t oedgeid, outs) {
           edge_data & edge = _g->edge_data(oedgeid);
           if (edge.time >= max_time){
             max_time = edge.time;
             vdata.num_edges = _g->target(oedgeid) - ps.M;
           }
         }
      }
  } 
   assert((int)last_items.size() == ps.N);

   for (int i=0; i< ac.K; i++){
      vertex_data & data = ps.times[i];
      data.pvec  = (ac.debug ? 0.1*ones(ac.D) : ::randu(ac.D)*factor);
   }
}

float libfm_predict(const vertex_data_libfm & user, const vertex_data_libfm & movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction, vec* sum);

float libfm_predict_new_user(const vertex_data_libfm& user, const vertex_data_libfm& movie, const edge_data * edge, const vertex_data * time, const float rating, float & prediction, vec * sum){
  prediction = ps.globalMean[0] + *movie.bias + time->bias;
  prediction = std::min((double)prediction, ac.maxval);
  prediction = std::max((double)prediction, ac.minval);
  float err = rating - prediction;
  return err*err;
}


void libfm_post_iter(){
  //printf("Entering last iter with %d\n", ps.iiter);

  post_iter_stats<graph_type>();  
  ac.libfm_rate *= ac.libfm_mult_dec;
}
/*
float libfm_predict(const vertex_data_libfm& user, 
                const vertex_data_libfm& movie, 
                const edge_data_mcmc * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction, vec * sum){ assert(false); }
*/

  
float libfm_predict(const vertex_data_libfm& user, 
                const vertex_data_libfm& movie, 
                const edge_data * edge,
                const vertex_data* time,
                const float rating, 
                float & prediction, vec * sum){
 
  if (*user.last_item == -1)
     return libfm_predict_new_user(user, movie, edge, time, rating, prediction, sum);

  vertex_data & last_item = last_items[*user.last_item];
  vec sum_sqr = zeros(ac.D);
  *sum = zeros(ac.D);
  prediction = ps.globalMean[0] + *user.bias + *movie.bias + time->bias + last_item.bias;
  for (int j=0; j< ac.D; j++){
    sum->operator[](j) += user.v[j] + movie.v[j] + time->pvec[j] + last_item.pvec[j];	 
    sum_sqr[j] = pow(user.v[j],2) + pow(movie.v[j],2) + pow(time->pvec[j],2) + pow(last_item.pvec[j],2); 
    prediction += 0.5 * (pow(sum->operator[](j),2) - sum_sqr[j]);
  }
  //truncate prediction to allowed values
  prediction = std::min((double)prediction, ac.maxval);
  prediction = std::max((double)prediction, ac.minval);
  //return the squared error
  float err = rating - prediction;
  assert(!std::isnan(err));
  return err*err; 
 
}
float libfm_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* time,
                const float rating, 
                float & prediction){
  vec sum; 
  return libfm_predict(vertex_data_libfm((vertex_data&)user), vertex_data_libfm((vertex_data&)movie), edge, time, rating, prediction, &sum);
}
 
 /***
 * UPDATE FUNCTION
 */
void libfm_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    


  int id = scope.vertex();
  /* GET current vertex data */
  vertex_data_libfm user(scope.vertex_data());
 
  /* print statistics */
  if (ps.to_print(id)){
    printf("%s: entering user node  %u \n", runmodesname[ps.algorithm], id);   
    debug_print_vec("U", user.v, ac.D);
  }
  *user.rmse = 0;
  gl_types::edge_list outs = scope.out_edge_ids();

  if (outs.size() == 0){
		if (id == ps.M-1){
  	  libfm_post_iter();
    }
		return; //if this user/movie have no ratings do nothing
  }
  

  timer t;
  t.start(); 
  assert(*user.last_item >= 0 && *user.last_item < ps.N);
  vertex_data & last_item = last_items[*user.last_item]; 
  
  foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data_libfm movie(scope.neighbor_vertex_data(scope.target(oedgeid)));
       
				float rui = edge.weight;
				float pui;
        vec sum;
        vertex_data & time = ps.times[(int)edge.time];
        float sqErr = libfm_predict(user, movie, &edge, &time, edge.weight, pui, &sum);
			  float eui = pui - rui;
                
        ps.globalMean[0] -= ac.libfm_rate * (eui + reg0 * ps.globalMean[0]);
        *user.bias -= ac.libfm_rate * (eui + ac.libfm_regw * *user.bias);
        *movie.bias -= ac.libfm_rate * (eui + ac.libfm_regw * *movie.bias);
        time.bias -= ac.libfm_regw * (eui + ac.libfm_regw * time.bias);
        last_item.bias -= ac.libfm_regw * (eui + ac.libfm_regw * last_item.bias);
        
        for(int f = 0; f < ac.D; f++){
                    // user
          float grad = sum[f] - user.v[f];
          user.v[f] -= ac.libfm_rate * (eui * grad + ac.libfm_regv * user.v[f]);
                    // item
          grad = sum[f] - movie.v[f];
          movie.v[f] -= ac.libfm_rate * (eui * grad + ac.libfm_regv * movie.v[f]);
                    // time
          grad = sum[f] - time.pvec[f];
          time.pvec[f] -= ac.libfm_rate * (eui * grad + ac.libfm_regv * time.pvec[f]);
                    // last item
          grad = sum[f] - last_item.pvec[f];
          last_item.pvec[f] -= ac.libfm_rate * (eui * grad + ac.libfm_regv * last_item.pvec[f]);

       }
                
				*user.rmse += sqErr;
		}

    ps.counter[EDGE_TRAVERSAL] += t.current_time();

    if (id == ps.M-1){
  	  libfm_post_iter();
    }


}
void fill_factors_libfm(){
		ps.U = zeros(ps.M,ac.D);
		ps.V = zeros(ps.N,ac.D);
    ps.T = zeros(ps.K,ac.D);
		ps.svdpp_usr_bias = zeros(ps.M);
		ps.svdpp_movie_bias = zeros(ps.N);
    ps.svdpp_time_bias = zeros(ac.K);      

		for (int i=0; i< ps.M+ps.N; i++){ 
      vertex_data &data = (vertex_data&)ps.g<graph_type>(TRAINING)->vertex_data(i);
      if (i < ps.M){
         set_row(ps.U, i, data.pvec);
         ps.svdpp_usr_bias[i] = data.bias;
      }
      else {
        set_row(ps.V, i-ps.M, data.pvec);
        ps.svdpp_movie_bias[i-ps.M] = data.bias;
      }
		}

   for (int i=0; i< ac.K; i++){
      vertex_data &data = ps.times[i];
      assert(data.pvec.size() == ac.D);
      set_row(ps.T, i, data.pvec);
      ps.svdpp_time_bias = data.bias;
   }
}

#include "graphlab/macros_undef.hpp"
#endif //__LIBFM_HPP
