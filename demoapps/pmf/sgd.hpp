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
void sgd_update_function(gl_types_svdpp::iscope &scope, 
			 gl_types_svdpp::icallback &scheduler) {
   assert(false);
} 

 /***
 * UPDATE FUNCTION
 */
void sgd_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  //USER NODES    
  if ((int)scope.vertex() < ps.M){


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
 
  
  /* print statistics */
  if (ac.debug&& (scope.vertex() == 0 || ((int)scope.vertex() == ps.M-1) || ((int)scope.vertex() == ps.M) || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("SGD: entering %s node  %u \n", (((int)scope.vertex() < ps.M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < ps.M) ? "V " : "U") , user.pvec, ac.D);
  }

  assert((int)scope.vertex() < ps.M+ps.N);

  user.rmse = 0;

  if (user.num_edges == 0){
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
      assert(!std::isnan(user.rmse));
      float err = edge.weight - estScore;
      movie.pvec = movie.pvec + ac.sgd_gamma*(err*user.pvec - ac.sgd_lambda*movie.pvec);
      user.pvec = user.pvec + ac.sgd_gamma*(err*movie.pvec - ac.sgd_lambda*user.pvec);
   }


   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (scope.vertex() == (uint)ps.M-1){
  	last_iter<gl_types::graph,vertex_data,edge_data>();
        ac.sgd_gamma *= ac.sgd_step_dec;
    }

  }

}

#include "graphlab/macros_undef.hpp"
#endif //__SGD_HPP
