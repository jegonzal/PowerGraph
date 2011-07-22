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


#ifndef _STATS_HPP
#define _STATS_HPP

#include <graphlab/macros_def.hpp>
extern int Lt;
extern mat dp;
extern double lasso_lambda;
extern runmodes algorithm;
//count the number of edges connecting a user/movie to its neighbors
//(when there are multiple edges in different times we count the total)
int count_edges(gl_types::edge_list es){
  
  if (algorithm != BPTF_TENSOR_MULT && algorithm != ALS_TENSOR_MULT)
      return es.size();

#ifndef GL_NO_MULT_EDGES
  int cnt = 0; 
  for (int j=0; j< (int)es.size(); j++){
    cnt += g->edge_data(es[j]).medges.size();
  }
  return cnt;
#else
  return es.size();
#endif
}


//go over all ratings and count how ratings for each node (user/movie)
void count_all_edges(graph_type * _g){
    for (int i=0; i<M+N; i++){
        vertex_data &vdata = _g->vertex_data(i);
        if (i < M)
          vdata.num_edges = count_edges(_g->out_edge_ids(i));
        else
          vdata.num_edges = count_edges(_g->in_edge_ids(i));
     }
}

int num_zeros(const vec & pvec){
   int ret = 0;
   for (int i=0; i< pvec.size(); i++)
      if (pvec[i] == 0)
         ret++;

   return ret;

}

// CALCULATE OBJECTIVE VALUE (Xiong paper)
double calc_obj(double res){
   
  double sumU = 0, sumV = 0, sumT = 0;
  double absSum = 0;
  int sparsity = 0;

  timer t;
  t.start(); 
  for (int i=0; i< M; i++){
    const vertex_data * data = &g->vertex_data(i);
    sumU += sum_sqr(data->pvec);
    if (algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || algorithm == ALS_SPARSE_USR_FACTOR){
	absSum += sum(abs(data->pvec));
        sparsity += num_zeros(data->pvec);
    }
     
  } 

  for (int i=M; i< M+N; i++){
    const vertex_data * data = &g->vertex_data(i);
    sumV += sum_sqr(data->pvec);
    if (algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || algorithm == ALS_SPARSE_MOVIE_FACTOR){
       absSum += sum(abs(data->pvec));
       sparsity += num_zeros(data->pvec);
    }
  } 


  mat T;
  if (tensor){
    T = zeros(D,K);
    for (int i=0; i<K; i++){
      vec tmp = times[i].pvec;
      sumT += pow(norm(tmp - vones, 2),2);
      T.set_col(i, tmp);
    }
    sumT *= pT;
  }
  counter[CALC_OBJ]+= t.current_time();
  
  double obj = (res + pU*sumU + pV*sumV + sumT + (tensor?trace(T*dp*T.transpose()):0)) / 2.0;

  if (algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || algorithm == ALS_SPARSE_USR_FACTOR || algorithm == ALS_SPARSE_MOVIE_FACTOR){ //add L1 penalty to objective
     int total;
     if (algorithm == ALS_SPARSE_MOVIE_FACTOR)
        total = N;
     else if (algorithm == ALS_SPARSE_USR_FACTOR)
        total = M;
     else if (algorithm == ALS_SPARSE_USR_MOVIE_FACTORS)
        total = M+N; 
     cout<<"Current sparsity : of " << runmodesname[algorithm] << " is: " << sparsity / ((double)D*total) << endl;
     obj += absSum;
  }

  if (debug)
     cout<<"OBJECTIVE: res: " << res << "sumU " << sumU << " sumV: " << sumV << " pu " << pU << " pV: " << pV << endl; 
  return obj;
}

// calc statistics about matrix/tensor and exit  
void calc_stats(testtype type){
   graph_type * gr = NULL;
   switch(type){ 
     case TRAINING: gr = g; break; 
     case VALIDATION: gr = &validation_graph; break;
     case TEST: gr = &test_graph; break;
   }

   if (gr->num_vertices() == 0){
     printf("%s is missing, skipping data\n", testtypename[type]);
     return;
   } 

   if (tensor && type == TRAINING){
     int firsttimeused=-1;
     int lasttimeused=-1;
     for (int i=0; i<K; i++){
       if (edges[i].size() > 0)
         firsttimeused = i;
     }
     for (int i=K-1; i>=0; i--){
       if (edges[i].size() > 0)
         lasttimeused = i;
     }
     printf("Out of total %d time components, first used is %d, last used is %d\n", K, firsttimeused, lasttimeused);
  }	

  double avgval=-1, minval=1e100, maxval=-1e100;
  double avgtime=-1, mintime=1e100, maxtime=-1e100;
  double minV=1e100, maxV=-1e100, minU=1e100, maxU=-1e100;
  int moviewithoutedges = 0;
  int userwithoutedges = 0;
  int numedges = 0;
  for (int i=M; i< M+N; i++){ 
    const vertex_data * data = &gr->vertex_data(i);
      if (itpp::min(data->pvec) < minU)
	 minU = itpp::min(data->pvec);
      if (itpp::max(data->pvec) > maxU)
	 maxU = itpp::max(data->pvec);
      if (gr->in_edge_ids(i).size() == 0)
	 moviewithoutedges++;
    foreach(edge_id_t iedgeid, gr->in_edge_ids(i)) {
#ifndef GL_NO_MULT_EDGES            
      multiple_edges & edges = gr->edge_data(iedgeid);
#endif
         //vertex_data * pdata = &gr->vertex_data(gr->source(iedgeid)); 
#ifndef GL_NO_MULT_EDGES        
      for (int j=0; j< (int)edges.medges.size(); j++){     
	edge_data & data = edges.medges[j];
#else
	edge_data & data = gr->edge_data(iedgeid);
#endif
	numedges++;
	avgval += data.weight;
	avgtime += data.time;
	if (data.weight<minval)
	   minval=data.weight;
	if (data.time <mintime)
	   mintime = data.time;
	if (data.weight>maxval)
	   maxval=data.weight;
	if (data.time > maxtime)
	   maxtime =data.time;
#ifndef GL_NO_MULT_EDGES	        
      }  
#endif
	
    }
 }
 for (int i=0; i< M; i++){ 
   const vertex_data * data = &gr->vertex_data(i);
   if (itpp::min(data->pvec) < minV)
     minV = itpp::min(data->pvec);
   if (itpp::max(data->pvec) > maxV)
     maxV = itpp::max(data->pvec);
	 
   if (gr->out_edge_ids(i).size() == 0)
     userwithoutedges++;
 }
 
 avgval /= numedges;
 avgtime /= numedges; 
 printf("%s Avg matrix value %g min val %g max value %g\n", testtypename[type],avgval, minval, maxval);
 printf("%s Avg time value %g min val %g max value %g\n", testtypename[type], avgtime, mintime, maxtime);
 printf("%s User without edges: %d movie without edges: %d\n", testtypename[type], userwithoutedges, moviewithoutedges);
 printf("%s Min V: %g Max V: %g Min U: %g, Max U: %g \n", testtypename[type], minV, maxV, minU, maxU);

 //verify we did not miss any ratings (nnz values)
 switch(type){
   case TRAINING: assert(numedges==L); break;
   case VALIDATION: assert(numedges==Le); break;
   case TEST: assert(numedges==Lt); break;
 }
}

//calculate RMSE. This function is called only before and after grahplab is run.
//during run, agg_rmse_by_movie is called 0 which is much lighter function (only aggregate sums of squares)
double calc_rmse(graph_type * _g, bool test, double & res){
     if (test && Le == 0)
       return NAN;
     
     if (algorithm == LANCZOS) //not implemented yet
       return NAN;
 
     res = 0;
     double RMSE = 0;
     int e = 0;
     for (int i=M; i< M+N; i++){
       vertex_data & data = g->vertex_data(i);
       foreach(edge_id_t iedgeid, _g->in_edge_ids(i)) {
         vertex_data & pdata = g->vertex_data(_g->source(iedgeid)); 
            
#ifndef GL_NO_MULT_EDGES
         multiple_edges & edges = _g->edge_data(iedgeid);
         for (int j=0; j< (int)edges.medges.size(); j++){       
           edge_data & edge = edges.medges[j];
#else
	        edge_data & edge = _g->edge_data(iedgeid);
#endif

           if (!ZERO)
           	assert(edge.weight != 0);

           float prediction = 0; 
           double sq_err = predict(data, 
                                   pdata, 
                                   ((algorithm == WEIGHTED_ALS) ? &edge: NULL), 
                                   tensor? (&times[(int)edge.time]):NULL, 
                                   edge.weight, 
                                   prediction);
           if (!ZERO)
	           assert(prediction != 0);         
           
           if (debug && (i== M || i == M+N-1) && ((e == 0) || ((e-1) == (test?Le:L))))
		cout<<"RMSE sq_err: " << sq_err << " prediction: " << prediction << endl; 

#ifndef GL_NO_MCMC
           if (BPTF && iiter > BURN_IN){
             edge.avgprd += prediction;
             sq_err = powf((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
           }
#endif
           if (algorithm == WEIGHTED_ALS)
              sq_err *= edge.time;
           RMSE+= sq_err;
           e++;
#ifndef GL_NO_MULT_EDGES        
         }
#endif
     }
   }
   res = RMSE;
   assert(e == (test?Le:L));
   return sqrt(RMSE/(double)e);

}
 

double calc_rmse_wrapper(graph_type* _g, bool test, double & res){
#ifdef SVD_PLUS_PLUS
   return calc_svd_rmse(_g, test, res);
#else
   if (algorithm == LANCZOS){
       res=-1; return -1; //not implemented yet
   }
   return calc_rmse(_g, test, res);
#endif
}

  
// go over all movie edges and aggregate RMSE by summing up the squares, computing the
// mean and then sqrt()
double agg_rmse_by_movie(double & res){

  timer t; t.start();
  res = 0;
  double RMSE = 0;
  for (int i=M; i< M+N; i++){ 
    const vertex_data * data = &g->vertex_data(i);
    RMSE+= data->rmse;
  }
  res = RMSE;
  counter[CALC_RMSE_Q] += t.current_time();
  return sqrt(RMSE/(double)L);

}
// go over all user edges and aggregate RMSE by summing up the squares, computing the
// mean and then sqrt()
double agg_rmse_by_user(double & res){

  timer t; t.start();
  res = 0;
  double RMSE = 0;
  for (int i=0; i< M; i++){ 
    const vertex_data * data = &g->vertex_data(i);
    RMSE+= data->rmse;
  }
  res = RMSE;
  counter[CALC_RMSE_Q] += t.current_time();
  return sqrt(RMSE/(double)L);

}

#include <graphlab/macros_undef.hpp>
#endif //_STATS_HPP
