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


#ifndef __BIAS_SGD_HPP
#define __BIAS_SGD_HPP

#include "graphlab.hpp"
#include "als.hpp"
#include "svdpp.hpp"
#include <graphlab/macros_def.hpp>


/**
 *
 *  Implementation of the bias SVD algorithm
 *
 * */
extern advanced_config ac;
extern problem_setup ps;

using namespace graphlab;

void bias_sgd_update_function(gl_types_mcmc::iscope &scope, 
			 gl_types_mcmc::icallback &scheduler) {
   assert(false);
} 
void bias_sgd_update_function(gl_types_mult_edge::iscope &scope, 
			 gl_types_mult_edge::icallback &scheduler) {
   assert(false);
}
template<typename graph_type>
void init_biassgd(graph_type* _g){
  assert(false);
}

template<>
void init_biassgd(graph_type *_g){
   fprintf(stderr, "%s %d factors\n", runmodesname[ps.algorithm], ac.D);
   double factor = 0.1/sqrt(ac.D);
#pragma omp parallel for
   for (int i=0; i<ps.M+ps.N; i++){
       vertex_data & vdata = _g->vertex_data(i);
       vdata.pvec = zeros(2*ac.D);
       vertex_data_svdpp data(_g->vertex_data(i));
       for (int j=0; j < ac.D; j++){
          data.weight[j] = (ac.debug ? 0.1 : (::randu()*factor));
          data.pvec[j] = (ac.debug ? 0.1 : (::randu()*factor));
       }
       vdata.bias = 0;
   } 
}

float bias_sgd_predict(const vertex_data_svdpp & user, const vertex_data_svdpp & movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction);

float bias_sgd_predict_new_user(const vertex_data_svdpp& user, const vertex_data_svdpp& movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction){
  prediction = ps.globalMean[0] + *movie.bias;
  prediction = std::min((double)prediction, ac.maxval);
  prediction = std::max((double)prediction, ac.minval);
  float err = rating - prediction;
  return err*err;
}


void bias_sgd_post_iter(){
  printf("Entering last iter with %d\n", ps.iiter);
  
  post_iter_stats<graph_type>();
 
  ac.sgd_gamma *= ac.sgd_step_dec;
}
float bias_sgd_predict(const vertex_data_svdpp& user, 
                const vertex_data_svdpp& movie, 
                const edge_data_mcmc * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){ assert(false); }

  
 
float bias_sgd_predict(const vertex_data_svdpp& user, 
                const vertex_data_svdpp& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
 
  if (user.num_edges == 0)
     return bias_sgd_predict_new_user(user, movie, edge, nothing, rating, prediction);

  prediction = ps.globalMean[0] + *user.bias + *movie.bias;
  for (int j=0; j< ac.D; j++)
    prediction += user.pvec[j]* movie.pvec[j];	

  //truncate prediction to allowed values
  prediction = std::min((double)prediction, ac.maxval);
  prediction = std::max((double)prediction, ac.minval);
  //return the squared error
  float err = rating - prediction;
  assert(!std::isnan(err));
  return err*err; 
 
}
float bias_sgd_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
   
  return bias_sgd_predict(vertex_data_svdpp((vertex_data&)user), vertex_data_svdpp((vertex_data&)movie), edge, nothing, rating, prediction);
}
 
 /***
 * UPDATE FUNCTION
 */
void bias_sgd_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    


  int id = scope.vertex();
  /* GET current vertex data */
  vertex_data_svdpp user(scope.vertex_data());
 
  /* print statistics */
  if (ps.to_print(id)){
    printf("biasSVD: entering user node  %u \n", id);   
    debug_print_vec("U", user.pvec, ac.D);
  }
  *user.rmse = 0;

  if (user.num_edges == 0){
		if (id == ps.M-1){
  	  bias_sgd_post_iter();
    }
		return; //if this user/movie have no ratings do nothing
  }

  gl_types::edge_list outs = scope.out_edge_ids();
  timer t;
  t.start(); 

   // for each rating
   //compute bias SVD Step 
   foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data_svdpp movie(scope.neighbor_vertex_data(scope.target(oedgeid)));
      float estScore;
      float sqErr = bias_sgd_predict(user, movie, &edge, NULL, edge.weight, estScore);
      *user.rmse += sqErr;
      assert(!std::isnan(*user.rmse));
      float err = edge.weight - estScore;
      *user.bias += ac.sgd_gamma*(err - ac.sgd_lambda* *user.bias);
      *movie.bias += ac.sgd_gamma*(err - ac.sgd_lambda* *movie.bias); 
      for (int j=0; j< ac.D; j++)
        movie.pvec [j] += ac.sgd_gamma*(err*user.pvec[j] - ac.sgd_lambda*movie.pvec[j]);
      for (int j=0; j< ac.D; j++)
				user.pvec[j] += ac.sgd_gamma*(err*movie.pvec[j] - ac.sgd_lambda*user.pvec[j]);
    }

    ps.counter[EDGE_TRAVERSAL] += t.current_time();

    if (id == ps.M-1){
  	  bias_sgd_post_iter();
    }


}

#include "graphlab/macros_undef.hpp"
#endif //__BIAS_SGD_HPP
