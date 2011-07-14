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
#include <graphlab/macros_def.hpp>


/**
 *
 *  Implementation of the SGD algorithm, as given in:
 *  Matrix Factorization Techniques for Recommender Systems
 *  by: Yehuda Koren, Robert Bell, Chris Volinsky
 *  In IEEE Computer, Vol. 42, No. 8. (07 August 2009), pp. 30-37. 
 *
 * */

extern string infile;
extern int iiter, L, Le;
extern bool ZERO;
extern timer gt;
extern graph_type validation_graph;
extern bool debug;
extern float sgd_lambda;
extern float sgd_gamma;
extern float sgd_step_dec;
using namespace graphlab;
using namespace itpp;

void last_iter();
double predict(const vertex_data& user, const vertex_data &movie, float rating, float & prediction);

/***
 * UPDATE FUNCTION
 */
void sgd_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  //USER NODES    
  if ((int)scope.vertex() < M){


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
 
  
  /* print statistics */
  if (debug&& (scope.vertex() == 0 || ((int)scope.vertex() == M-1) || ((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1))){
    printf("SGD: entering %s node  %u \n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , user.pvec, D);
  }

  assert((int)scope.vertex() < M+N);

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
#ifndef GL_NO_MULT_EDGES
      multiple_edges & edges = scope.edge_data(oedgeid);
      for (int j=0; j< (int)edges.medges.size(); j++){
         edge_data & edge = edges.medges[j];
#else
      edge_data & edge = scope.edge_data(oedgeid);
#endif
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      float estScore;
      float sqErr = predict(user, movie, NULL, edge.weight, estScore);
      user.rmse += sqErr;
      assert(!isnan(user.rmse));
      float err = edge.weight - estScore;
      movie.pvec = movie.pvec + sgd_gamma*(err*user.pvec - sgd_lambda*movie.pvec);
      user.pvec = user.pvec + sgd_gamma*(err*movie.pvec - sgd_lambda*user.pvec);
#ifndef GL_NO_MULT_EDGES
      }
#endif
   }


   counter[EDGE_TRAVERSAL] += t.current_time();

   if (scope.vertex() == (uint)M-1){
  	last_iter();
        sgd_gamma *= sgd_step_dec;
    }

  }

}

#include "graphlab/macros_undef.hpp"
#endif //__SGD_HPP
