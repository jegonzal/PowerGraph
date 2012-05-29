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

struct vertex_data_libfm{
   float * bias;
   double * v;
   int num_edges; 
   float * rmse;

   vertex_data_libfm(const vertex_data & vdata){
     v = (double*)&vdata.pvec[0];
     bias = (float*)&vdata.bias;
     num_edges = vdata.num_edges;
     rmse = (float*)& vdata.rmse;
   }
   
   vertex_data_libfm & operator=(vertex_data & data){
     v = (double*)&data.pvec[0];
     bias = (float*)&data.bias;
     num_edges = data.num_edges;
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
          data.v[j] = (ac.debug ? 0.1 : (randu()*factor));
       }
   } 
}

float libfm_predict(const vertex_data_libfm & user, const vertex_data_libfm & movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction, vec* sum);

float libfm_predict_new_user(const vertex_data_libfm& user, const vertex_data_libfm& movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction, vec * sum){
  prediction = ps.globalMean[0] + *movie.bias;
  prediction = std::min((double)prediction, ac.maxval);
  prediction = std::max((double)prediction, ac.minval);
  float err = rating - prediction;
  return err*err;
}


void libfm_post_iter(){
  //printf("Entering last iter with %d\n", ps.iiter);

  post_iter_stats<graph_type>();  
  ac.libfm_rate *= ac.libfm_mult_dec;
  ps.iiter++;
}

float libfm_predict(const vertex_data_libfm& user, 
                const vertex_data_libfm& movie, 
                const edge_data_mcmc * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction, vec * sum){ assert(false); }

  
 
float libfm_predict(const vertex_data_libfm& user, 
                const vertex_data_libfm& movie, 
                const edge_data * edge,
                const vertex_data* time,
                const float rating, 
                float & prediction, vec * sum){
 
  if (user.num_edges == 0)
     return libfm_predict_new_user(user, movie, edge, time, rating, prediction, sum);

  vec sum_sqr = zeros(ac.D);
  *sum = zeros(ac.D);
  prediction = ps.globalMean[0] + *user.bias + *movie.bias + time->bias;
  for (int j=0; j< ac.D; j++){
    sum->operator[](j) += user.v[j] + movie.v[j] + time->pvec[j];	 //TODO + last item
    sum_sqr[j] = pow(user.v[j],2) + pow(movie.v[j],2) + pow(time->pvec[j],2); //TODO + last item
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
                float & prediction, vec * sum){
   
  return libfm_predict(vertex_data_libfm((vertex_data&)user), vertex_data_libfm((vertex_data&)movie), edge, time, rating, prediction, sum);
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

  if (user.num_edges == 0){
		if (id == ps.M-1){
  	  libfm_post_iter();
    }
		return; //if this user/movie have no ratings do nothing
  }

  gl_types::edge_list outs = scope.out_edge_ids();
  timer t;
  t.start(); 
            

  int max_time = 0;
  graphlab::edge_id_t lastItemRating = -1;
  foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      if (edge.time > max_time){
          max_time = edge.time;
          lastItemRating = oedgeid;
       }
  } 
  edge_data & last_rating = scope.edge_data(lastItemRating);
 
  foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data_libfm movie(scope.neighbor_vertex_data(scope.target(oedgeid)));
       
				float rui = edge.weight;
				float pui;
        vec sum;
        float sqErr = libfm_predict(user, movie, &last_rating, NULL, edge.weight, pui, &sum);
        vertex_data & time = ps.times[(int)edge.time];
			  float eui = pui - rui;
                
        ps.globalMean[0] -= ac.libfm_rate * (eui + reg0 * ps.globalMean[0]);
        /*        w[user_idx] -= libfm_lrate * (eui + libfm_regw * w[user_idx]);
                w[item_idx] -= libfm_lrate * (eui + libfm_regw * w[item_idx]);
                w[time_idx] -= libfm_lrate * (eui + libfm_regw * w[time_idx]);
                w[last_item_idx] -= libfm_lrate * (eui + libfm_regw * w[last_item_idx]);*/
        *user.bias -= ac.libfm_rate * (eui + ac.libfm_regw * *user.bias);
        *movie.bias -= ac.libfm_rate * (eui + ac.libfm_regw * *movie.bias);
        time.bias -= ac.libfm_regw * (eui + ac.libfm_regw * time.bias);
        //TODO *last_rating.bias -= ac.libfm_regw * (eru + ac.libfm_regw * *last_rating.bias);
        
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
                    //grad = sum[f] - v[last_item_idx][f];
                    //v[last_item_idx][f] -= libfm_rate * (eui * grad + ac.libfm_regv * v[last_item_idx][f]);

                }
                
				*user.rmse += sqErr;
			}

    ps.counter[EDGE_TRAVERSAL] += t.current_time();

    if (id == ps.M-1){
  	  libfm_post_iter();
    }


}

#include "graphlab/macros_undef.hpp"
#endif //__LIBFM_HPP
