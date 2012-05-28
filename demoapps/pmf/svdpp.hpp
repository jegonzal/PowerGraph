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
 * Written By Danny Bickson, CMU
 * Based on Code by Yehuda Koren, Yahoo! Research
 * Send any question / comments to: danny.bickson@gmail.com
 *
 * This code implements the paper: Factorization Meets the Neighborhood: a Multifaceted 
 * Collaborative Filtering Model by Yehuda Koren, in KDD 2008.
 * Parallelization of the code is done by Danny Bickson, CMU
 */


#ifndef __SVDPP_HPP
#define __SVDPP_HPP

#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>
#include "pmf.h"

extern advanced_config ac;
extern problem_setup ps;

using namespace graphlab;


//constructor


vertex_data_svdpp::vertex_data_svdpp(vertex_data & vdata){
    rmse = &vdata.rmse;
    num_edges = vdata.num_edges;
    bias = &vdata.bias;
    assert(vdata.pvec.size() == 2*ac.D); //TO REMOVE
    pvec = &vdata.pvec[0];
    weight = &vdata.pvec[ac.D];
}


vertex_data_svdpp & vertex_data_svdpp::operator=(vertex_data & vdata){
    rmse = &vdata.rmse;
    num_edges = vdata.num_edges;
    bias = &vdata.bias;
    assert(vdata.pvec.size() == 2*ac.D); //TO REMOVE
    pvec = &vdata.pvec[0];
    weight = &vdata.pvec[ac.D];
    return *this;
}
/*
void vertex_data_svdpp::save(graphlab::oarchive& archive) const {  
    ////TODO archive << pvec;
    archive << rmse << num_edges << bias; 
    ///TODO archive << weight;
  }  
   
void vertex_data_svdpp::load(graphlab::iarchive& archive) {  
     //TODO archive >> pvec;
     archive >> rmse >> num_edges >> bias;  
     //TODO archive >> weight;
}*/

template<typename graph_type>
void init_svdpp(graph_type* _g){
  assert(false);
}

template<>
void init_svdpp/*<gc>*/(graph_type/*_svdpp*/ *_g){
   fprintf(stderr, "SVD++ %d factors\n", ac.D);
   double factor = 0.1/sqrt(ac.D);
#pragma omp parallel for
   for (int i=0; i<ps.M+ps.N; i++){
       vertex_data & vdata = _g->vertex_data(i);
       vdata.pvec = zeros(2*ac.D);
       vertex_data_svdpp data(_g->vertex_data(i));
       for (int j=0; j < ac.D; j++){
          data.weight[j] = (ac.debug ? 0.1 : (randu()*factor));
          data.pvec[j] = (ac.debug ? 0.1 : (randu()*factor));
       }
       data.bias = 0;
   } 
}

float svdpp_predict_new_user(const vertex_data_svdpp& user, const vertex_data_svdpp& movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction){
  prediction = ps.globalMean[0] + *movie.bias;
  prediction = std::min((double)prediction, ac.maxval);
  prediction = std::max((double)prediction, ac.minval);
  float err = rating - prediction;
  return err*err;
}




float svdpp_predict(const vertex_data_svdpp& user, const vertex_data_svdpp& movie, const edge_data_mcmc * edge, const vertex_data * nothing, const float rating, float & prediction){
  assert(false);
}
float svdpp_predict(const vertex_data_svdpp& user, const vertex_data_svdpp& movie, const edge_data * edge, const vertex_data * nothing, const float rating, float & prediction){
      assert(nothing == NULL);
      assert(ps.algorithm == SVD_PLUS_PLUS);

       if (user.num_edges == 0)
        return svdpp_predict_new_user(user, movie, edge, nothing, rating, prediction);

           //\hat(r_ui) = \mu + 
        prediction = ps.globalMean[0];
                 // + b_u  +    b_i +
        prediction += *user.bias + *movie.bias;
                 // + q_i^T   *(p_u      +sqrt(|N(u)|)\sum y_j)
        //prediction += dot_prod(movie.pvec,(user.pvec+user.weight));
        for (int j=0; j< ac.D; j++)
          prediction += movie.pvec[j] * (user.pvec[j] + user.weight[j]);

        prediction = std::min((double)prediction, ac.maxval);
        prediction = std::max((double)prediction, ac.minval);
        float err = rating - prediction;
        if (std::isnan(err))
          logstream(LOG_FATAL)<<"Got into numerical errors. Try to decrease step size using svdpp command line flags (see ./pmf --help for full list)" << std::endl;
        return err*err; 
      
}

float svdpp_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
  return svdpp_predict(vertex_data_svdpp((vertex_data&)user), vertex_data_svdpp((vertex_data&)movie), edge, nothing, rating, prediction);
}


/*void predict_missing_value(const vertex_data_svdpp&data, const vertex_data_svdpp& pdata, edge_data& edge, double & sq_err, int&e, int i){
    float prediction = 0;
    svdpp_predict(data, pdata, &edge, NULL, edge.weight, prediction);
    e++;
}*/
 

//calculate RMSE. This function is called only before and after grahplab is run.
//during run, agg_rmse_by_movie is called 0 which is much lighter function (only aggregate sums of squares)
double calc_svd_rmse(const graph_type * _g, bool test, double & res){

     graph_type * g = (graph_type*)ps.g<graph_type>(TRAINING);

     if (test && ps.Le == 0)
       return NAN;
      
     
     res = 0;
     double sqErr =0;
     int nCases = 0;

     for (int i=0; i< ps.M; i++){
       vertex_data_svdpp usr(g->vertex_data(i));
       int n = usr.num_edges; //+1.0 ? //regularization
       if (n == 0){}
       else {
         memset( usr.weight, 0 , ac.D * sizeof(double));
         foreach(edge_id_t oedgeid, g->out_edge_ids(i)) {
           vertex_data_svdpp movie(g->vertex_data(g->target(oedgeid))); 
           for (int j=0; j < ac.D; j++)
					   usr.weight[j] += movie.weight[j];
         }
         float usrnorm = double(1.0/sqrt(n));
         for (int j=0; j < ac.D; j++)
           usr.weight[j] *= usrnorm;
       }

       foreach(edge_id_t oedgeid, _g->out_edge_ids(i)){
         const edge_data & item = _g->edge_data(oedgeid);
         const vertex_data_svdpp movie(g->vertex_data(_g->target(oedgeid))); 
         float estScore;
        sqErr += svdpp_predict(usr, movie, (edge_data*)NULL, NULL, item.weight, estScore);
         nCases++;
       }
   }
   res = sqErr;
   assert(nCases == (test?ps.Le:ps.L));
   return sqrt(sqErr/(double)nCases);
}


void svd_post_iter(){
  printf("Entering last iter with %d\n", ps.iiter);

  double res,res2;
  double training_rmse = agg_rmse_by_user<graph_type, vertex_data>(res);
  double validation_rmse =  calc_svd_rmse(ps.g<graph_type>(VALIDATION), true, res2);
  printf("%g) Iter %s %d, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", ps.gt.current_time(), "SVD", ps.iiter,  training_rmse, validation_rmse);

  if (ac.calc_ap){
     logstream(LOG_INFO)<<"AP@3 for training: " << calc_ap<graph_type,vertex_data,edge_data>(ps.g<graph_type>(TRAINING)) << " AP@3 for validation: " << calc_ap<graph_type,vertex_data,edge_data>(ps.g<graph_type>(VALIDATION)) << std::endl;
  }

  //stop on divergence
  if (ac.halt_on_rmse_increase)
    if ((ps.validation_rmse && (ps.validation_rmse < validation_rmse)) ||
        (ps.training_rmse && (ps.training_rmse < training_rmse)))
          dynamic_cast<graphlab::core<vertex_data_svdpp,edge_data>*>(ps.glcore)->engine().stop();

  ps.validation_rmse = validation_rmse; 
  ps.training_rmse = training_rmse;

  ac.svdp.itmFctrStep *= ac.svdpp_step_dec;
  ac.svdp.itmFctr2Step *= ac.svdpp_step_dec;
  ac.svdp.usrFctrStep *= ac.svdpp_step_dec;
  ac.svdp.itmBiasStep *= ac.svdpp_step_dec;
  ac.svdp.usrBiasStep *= ac.svdpp_step_dec;

  ps.iiter++;
}


void svd_plus_plus_update_function(gl_types_mult_edge::iscope &scope, 
			 gl_types_mult_edge::icallback &scheduler) {
   assert(false); //mode not supported
} 
void svd_plus_plus_update_function(gl_types_mcmc::iscope &scope, 
			 gl_types_mcmc::icallback &scheduler) {
   assert(false); //mode not supported
} 
 
/***
 * UPDATE FUNCTION
 */
void svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    
  /* GET current vertex data */
  vertex_data_svdpp user(scope.vertex_data());
  int id = scope.vertex(); 
  
  /* print statistics */
  if (ps.to_print(id)){
    printf("SVDPP: entering user node  %u \n", id);   
    debug_print_vec("U" , user.pvec, ac.D);
  }

  *user.rmse = 0;
  if (user.num_edges == 0){
   if (id == ps.M-1)
    	svd_post_iter();
    return; //if this user/movie have no ratings do nothing
  }

  gl_types::edge_list outs = scope.out_edge_ids();
  timer t;

  t.start(); 
    memset(user.weight, 0, ac.D*sizeof(double));
    foreach(graphlab::edge_id_t oedgeid, outs) {
      vertex_data_svdpp movie(scope.neighbor_vertex_data(scope.target(oedgeid))); 
      //sum_{j \in N(u)} y_j 
      for (int j=0; j < ac.D; j++)
        user.weight[j] += movie.weight[j]; 
    }

   // sqrt(|N(u)|) 
   assert(user.num_edges != 0);  
   float usrNorm = double(1.0/sqrt(user.num_edges));
   //sqrt(|N(u)| * sum_j y_j
   for (int j=0; j< ac.D; j++)
     user.weight[j] *= usrNorm;

   vec step = zeros(ac.D);
 
   // main algorithm, see Koren's paper, just below below equation (16)
   foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data_svdpp movie(scope.neighbor_vertex_data(scope.target(oedgeid)));
      float estScore;
      *user.rmse += svdpp_predict(user, movie, (edge_data*)NULL, NULL, edge.weight, estScore); 
      // e_ui = r_ui - \hat{r_ui}
      float err = edge.weight - estScore;
      assert(!std::isnan(*user.rmse));
      vec itmFctr = _vec(movie.pvec, ac.D);
      vec usrFactor = _vec(user.pvec, ac.D) ;
   
      //q_i = q_i + gamma2     *(e_ui*(p_u      +  sqrt(N(U))\sum_j y_j) - gamma7    *q_i)
      for (int j=0; j< ac.D; j++)
        movie.pvec[j] += ac.svdp.itmFctrStep*(err*(usrFactor[j] +  user.weight[j])             - ac.svdp.itmFctrReg*itmFctr[j]);
      //p_u = p_u + gamma2    *(e_ui*q_i   -gamma7     *p_u)
      for (int j=0; j< ac.D; j++)
				user.pvec[j] += ac.svdp.usrFctrStep*(err *itmFctr[j] - ac.svdp.usrFctrReg*usrFactor[j]);
      step += err*itmFctr;

      //b_i = b_i + gamma1*(e_ui - gmma6 * b_i) 
      *movie.bias += ac.svdp.itmBiasStep*(err-ac.svdp.itmBiasReg* *movie.bias);
      //b_u = b_u + gamma1*(e_ui - gamma6 * b_u)
      *user.bias += ac.svdp.usrBiasStep*(err-ac.svdp.usrBiasReg* *user.bias);
   }

   step *= float(ac.svdp.itmFctr2Step*usrNorm);
   //gamma7 
   double mult = ac.svdp.itmFctr2Step*ac.svdp.itmFctr2Reg;
   foreach(graphlab::edge_id_t oedgeid, outs){
      vertex_data_svdpp  movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      //y_j = y_j  +   gamma2*sqrt|N(u)| * q_i - gamma7 * y_j
      for (int j=0; j< ac.D; j++)
        movie.weight[j] +=  step[j]                    -  mult  * movie.weight[j];
   }


   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (id == ps.M-1)
  	svd_post_iter();

}
void fill_factors_svdpp(){
		ps.U = zeros(ps.M,ac.D);
		ps.V = zeros(ps.N,ac.D);
		ps.svdpp_usr_bias = zeros(ps.M);
		ps.svdpp_movie_bias = zeros(ps.N);
      
		for (int i=0; i< ps.M+ps.N; i++){ 
      vertex_data_svdpp data = (vertex_data&)ps.g<graph_type>(TRAINING)->vertex_data(i);
      if (i < ps.M){
         for (int j=0; j< ac.D; j++)
            set_val(ps.U, i, j, data.pvec[j] + data.weight[j]);
         ps.svdpp_usr_bias[i] = *data.bias;
      }
      else {
        set_row(ps.V, i-ps.M, _vec(data.pvec, ac.D));
        ps.svdpp_movie_bias[i-ps.M] = *data.bias;
      }
		}
}
#include "graphlab/macros_undef.hpp"
#endif //__SVDPP_HPP
