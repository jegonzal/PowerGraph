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


#ifndef __SGD_HPP
#define __SGD_HPP

#include "graphlab.hpp"
#include "als.hpp"
#include <graphlab/macros_def.hpp>


/**
 *
 *  Implementation of the SGD algorithm, as given in:
 *  Matrix Factorization Techniques for Recommender Systems
 *  by: Yehuda Koren, Robert Bell, Chris Volinsky
 *  In IEEE Computer, Vol. 42, No. 8. (07 August 2009), pp. 30-37. 
 *
 * */
extern advanced_config ac;
extern problem_setup ps;

using namespace graphlab;

void sgd_update_function(gl_types_mcmc::iscope &scope, 
			 gl_types_mcmc::icallback &scheduler) {
   assert(false);
} 
void sgd_update_function(gl_types_mult_edge::iscope &scope, 
			 gl_types_mult_edge::icallback &scheduler) {
   assert(false);
}

 /***
 * UPDATE FUNCTION
 */
void sgd_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  int id = scope.vertex();
  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  if (ps.to_print(id)){
    printf("SVDPP: entering user node  %u \n", id);   
    debug_print_vec("U" , user.pvec, ac.D);
  } 
  user.rmse = 0;
  if (user.num_edges == 0){
  	 if (id == ps.M - 1){
       last_iter<graph_type>();
       ac.sgd_gamma *= ac.sgd_step_dec;
     }
     return; //if this user/movie have no ratings do nothing
  }


  gl_types::edge_list outs = scope.out_edge_ids();
  timer t;
  t.start(); 

   // for each rating
   //compute SGD Step 
   foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      float estScore;
      float sqErr = predict(user, movie, NULL, edge.weight, estScore);
      user.rmse += sqErr;
      if (std::isnan(user.rmse))
        logstream(LOG_FATAL)<<"Numerical error occured. Aborting run." << std::endl;
      float err = edge.weight - estScore;
      vec oldmovie = movie.pvec;
      vec olduser = user.pvec;
      movie.pvec = movie.pvec + ac.sgd_gamma*(err*olduser - ac.sgd_lambda*oldmovie);
      user.pvec = user.pvec + ac.sgd_gamma*(err*oldmovie - ac.sgd_lambda*olduser);
   }


   ps.counter[EDGE_TRAVERSAL] += t.current_time();
 
   if (id == ps.M-1){
  	last_iter<graph_type>();
        ac.sgd_gamma *= ac.sgd_step_dec;
    }


}

#include "graphlab/macros_undef.hpp"
#endif //__SGD_HPP
