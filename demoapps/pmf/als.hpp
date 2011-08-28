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


#ifndef _ALS_HPP
#define _ALS_HPP

#include "pmf.h"
#include <graphlab/macros_def.hpp>

extern advanced_config ac;
extern problem_setup ps;
mat eDT; 
vec vones; 
double pU = 10; //regularization for users
double pV = 10; //regularization for movies


vec CoSaMP(mat Phi, vec u, int K, int max_iter, double tol1, int D);

void init_pmf() {
  if (ps.BPTF)
    ps.pT=10;
  eDT = itpp::eye(ac.D)*ps.pT;
  vones = itpp::ones(ac.D);
  printf("pU=%g, pV=%g, pT=%g, D=%d\n", pU, pV, ps.pT,ac.D);  
}
 
 
// fill out the linear relation matrix (Q) between users/movies and movie/users
inline void parse_edge(const edge_data& edge, const vertex_data & pdata, mat & Q, vec & vals, int i, vec * weights){
      
  if (!ac.zero)
  	assert(edge.weight != 0);

  if (ps.tensor){
    dot2(pdata.pvec,  ps.times[(int)edge.time].pvec, Q, i, ac.D);  
  }
  else {
    for (int j=0; j<ac.D; j++)
      Q.set(j,i, pdata.pvec[j]); 
  }
 
  vals[i] = edge.weight;
  
  if (weights != NULL){
     if (!ac.zero) assert(edge.time!= 0);
     weights->set( i, edge.time);
     vals[i] *= edge.time;
  }
}
 /***
 * UPDATE FUNCTION
 */
void user_movie_nodes_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& vdata = scope.vertex_data();
 
  int id = scope.vertex();
  bool toprint = ac.debug && (id == 0 || (id == ps.M-1) || (id == ps.M) || (id == ps.M+ps.N-1)); 
  bool isuser = id < ps.M;
  /* print statistics */
  if (toprint){
    printf("entering %s node  %u \n", (!isuser ? "movie":"user"), id);   
    debug_print_vec((isuser ? "V " : "U") , vdata.pvec, ac.D);
  }

  vdata.rmse = 0;

  int numedges = vdata.num_edges;
  if (numedges == 0){
     if (!ps.tensor && (int)scope.vertex() == ps.M+ps.N-1)
       last_iter();

   return; //if this user/movie have no ratings do nothing
  }


  gl_types::edge_list outs = scope.out_edge_ids();
  gl_types::edge_list ins = scope.in_edge_ids();
  timer t;
  mat Q(ac.D,numedges); //linear relation matrix
  vec vals(numedges); //vector of ratings
  vec weight(numedges); // vector of weights (to be used in weighted ALS)
  int i=0;

  t.start(); 
  //USER NODES    
  if (isuser){

    foreach(gl_types::edge_id oedgeid, outs) {
      const vertex_data  & pdata = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
#ifndef GL_NO_MULT_EDGES
      multiple_edges &medges =scope.edge_data(oedgeid);
      for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
#else
        edge_data & edge = scope.edge_data(oedgeid);
#endif
        //go over each rating of a movie and put the movie vector into the matrix Q
        //and vector vals
        parse_edge(edge, pdata, Q, vals, i, ps.algorithm == WEIGHTED_ALS? &weight : NULL); 
        if (toprint && (i==0 || i == numedges-1))
          std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;
        i++;
#ifndef GL_NO_MULT_EDGES
      }   
#endif
    }
  }

  else {


    //MOVIE NODES
    foreach(gl_types::edge_id iedgeid, ins) {

      const vertex_data & pdata = scope.const_neighbor_vertex_data(scope.source(iedgeid)); 
#ifndef GL_NO_MULT_EDGES
       multiple_edges & medges =scope.edge_data(iedgeid);
        for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
#else
	edge_data & edge = scope.edge_data(iedgeid);
#endif   
        //go over each rating by user
        parse_edge(edge, pdata, Q, vals, i, ps.algorithm == WEIGHTED_ALS ? &weight: NULL); 
        if (toprint/* && (i==0 || i == numedges-1)*/)
          std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;

        i++;
        float prediction;     
        float trmse = predict(vdata, 
                              pdata, 
                              (ps.algorithm == WEIGHTED_ALS) ? &edge : NULL, 
                              ps.tensor?(&ps.times[(int)edge.time]):NULL, 
			      edge.weight, 
                              prediction);
#ifndef GL_NO_MCMC
        if (BPTF && iiter > BURN_IN){
          edge.avgprd += prediction;        
          trmse = pow((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
        }
#endif
      //weight rmse with edge weight 
      if (ps.algorithm == WEIGHTED_ALS)
          trmse *= edge.time;
      
       if (toprint)
          cout<<"trmse: " << trmse << endl;
      //aggregate RMSE
       vdata.rmse += trmse; 
 
#ifndef GL_NO_MULT_EDGES     
      }
#endif
      
    }

  }
  assert(i == numedges);
  ps.counter[EDGE_TRAVERSAL] += t.current_time();

  vec result;
      
  if (!ps.BPTF){
    //COMPUTE LEAST SQUARES (ALTERNATING LEAST SQUARES)
    t.start();
    double regularization = ac.als_lambda;
  
    //compute weighted regularization (see section 3.2 of Zhou paper)
   if (!ac.regnormal)
	   regularization*= Q.cols();

   //enforce sparsity priors on resulting factor vector, see algorithm 1, page 4 in Xi et. al paper
   if (ps.algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || (ps.algorithm == ALS_SPARSE_USR_FACTOR && isuser) || 
      (ps.algorithm == ALS_SPARSE_MOVIE_FACTOR && !isuser)){ 
       double sparsity_level = 1.0;
       if (isuser)
	  sparsity_level -= ac.user_sparsity;
       else sparsity_level -= ac.movie_sparsity;
       result = CoSaMP(Q*itpp::transpose(Q)+eDT*regularization, Q*vals, ceil(sparsity_level*(double)ac.D), ac.lasso_max_iter, 1e-4, ac.D); 
   }
    // compute regular least suqares
   else if (ps.algorithm != WEIGHTED_ALS){
       bool ret = itpp::ls_solve_chol(Q*itpp::transpose(Q)+eDT*regularization, Q*vals, result);
       assert(ret);
    } 
    //Weighted alternating least squares (see equations (6),(7) in paper 9)
    else {
       vec b = Q*vals;
       weight = sqrt(weight);
       //avoid explicit creation of W = diag(W) because it is wasteful in memory.
       //instead, compute directly the product Q*W*Q'
       for (int i=0; i<ac.D; i++)
	 for (int j=0; j<numedges; j++)
             Q._elem(i,j)*= weight[j];
       mat A = Q*transpose(Q)+(eDT*regularization);
       bool ret = itpp::ls_solve_chol(A, b, result);
       if (ac.debug)
          cout<<" eDT : " << eDT << "reg: " << regularization << " Q*vals " << b << "Q*W*Q'+eDT+Reg: " << A << endl;
       assert(ret);
    }
    ps.counter[ALS_LEAST_SQUARES] += t.current_time();
  }
  else 
  {
    //COMPUTE LEAST SQUARES (BPTF)
    //according to equation A.6 or A.7 in Xiong paper.
    assert(Q.rows() == ac.D);
    t.start();
    mat iAi_;
    bool ret =inv((isuser? A_U : A_V) + alpha *  Q*itpp::transpose(Q), iAi_);
    assert(ret);
    t.start();
    vec mui_ =  iAi_*((isuser? (A_U*mu_U) : (A_V*mu_V)) + alpha * Q * vals);
    ps.counter[BPTF_LEAST_SQUARES2]+= t.current_time();
       
    t.start();
    result = mvnrndex(mui_, iAi_, ac.D); 
    assert(result.size() == ac.D);
    ps.counter[BPTF_MVN_RNDEX] += t.current_time();
  }

  if (toprint){
    std::cout <<(ps.BPTF?"BPTF":"ALS")<<std::endl<<" result: " << result << " edges: " << numedges << std::endl;
  }
      
  //store new result
  vdata.pvec =  result;

  //calc post round tasks
  if (!ps.tensor && (int)scope.vertex() == ps.M+ps.N-1)
    last_iter();

}




#include <graphlab/macros_undef.hpp>
#endif //_ALS_HPP
