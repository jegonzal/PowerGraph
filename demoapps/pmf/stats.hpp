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
#include "pmf.h"

extern advanced_config ac;
extern problem_setup ps;

using namespace graphlab;

extern const char * testtypename[];
extern const char * countername[];
extern std::vector<edge_id_t> * edges;


void print_runtime_counters(){
  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (ps.counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], ps.counter[i]); 
  }
}

//count the number of edges connecting a user/movie to its neighbors
//(when there are multiple edges in different times we count the total)
int count_edges(gl_types::edge_list es, const graph_type *_g){
   return es.size();
}
int count_edges(gl_types_mcmc::edge_list es, const graph_type_mcmc *_g){
   return es.size();
}
int count_edges(gl_types_svdpp::edge_list es, const graph_type_svdpp *_g){
   return es.size();
}

//count the number of edges connecting a user/movie to its neighbors
//(when there are multiple edges in different times we count the total)
int count_edges(graphlab::graph<vertex_data, multiple_edges>::edge_list es, const graphlab::graph<vertex_data, multiple_edges>*& _g){
  int cnt = 0; 
  for (int j=0; j< (int)es.size(); j++){
    cnt += _g->edge_data(es[j]).medges.size();
  }
  return cnt;
}




//go over all ratings and count how ratings for each node (user/movie)
template<typename graph_type>
void count_all_edges(const graph_type * _g){
    for (int i=0; i<ps.M+ps.N; i++){
        vertex_data &vdata = (vertex_data&)_g->vertex_data(i);
        if (i < ps.M)
          vdata.num_edges = count_edges(_g->out_edge_ids(i),_g);
        else
          vdata.num_edges = count_edges(_g->in_edge_ids(i),_g);
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
template<typename graph_type, typename vertex_data>
double calc_obj(double res){
   
  double sumU = 0, sumV = 0, sumT = 0;
  double absSum = 0;
  int user_sparsity = 0;
  int movie_sparsity = 0;
  int user_cnt = 0; 
  int movie_cnt = 0;
  int edges = 0;
  const graph_type * g = ps.g<graph_type>(TRAINING);

  graphlab::timer t;
  t.start();
  for (int i=0; i< ps.M; i++){
    const vertex_data * data = &g->vertex_data(i);
    if (data->num_edges > 0){
      sumU += sum_sqr(data->pvec);
      if (ps.algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || ps.algorithm == ALS_SPARSE_USR_FACTOR){
	absSum += abs_sum(data->pvec);
      } 
      user_sparsity += num_zeros(data->pvec);
      user_cnt++;
      edges+= data->num_edges;
    }
     
  } 

  for (int i=ps.M; i< ps.M+ps.N; i++){
    const vertex_data * data = &g->vertex_data(i);
    if (data->num_edges > 0 ){
      sumV += sum_sqr(data->pvec);
      if (ps.algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || ps.algorithm == ALS_SPARSE_MOVIE_FACTOR){
         absSum += abs_sum(data->pvec);
      }  
      movie_sparsity += num_zeros(data->pvec);
      movie_cnt++;
     edges+= data->num_edges;
    }
  } 


  mat T;
  if (ps.tensor){
    T = zeros(ac.D,ps.K);
    for (int i=0; i<ps.K; i++){
      vec tmp = ps.times[i].pvec;
      sumT += pow(norm(tmp - ps.vones, 2),2);
      set_col(T, i, tmp);
    }
    sumT *= ps.pT;
  }
  ps.counter[CALC_OBJ]+= t.current_time();
  
  double obj = (res +ps.pU*sumU + ps.pV*sumV + sumT + (ps.tensor?trace(T*ps.dp*T.transpose()):0)) / 2.0;

  if (ps.algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || ps.algorithm == ALS_SPARSE_USR_FACTOR || ps.algorithm == ALS_SPARSE_MOVIE_FACTOR){ //add L1 penalty to objective
     cout<<"Current user sparsity : " << ((double)user_sparsity / ((double)ac.D*user_cnt)) << " movie sparsity: "  << ((double)movie_sparsity / ((double)ac.D*movie_cnt)) << endl;
     obj += absSum;
  }

  if (ac.debug)
     cout<<"OBJECTIVE: res: " << res << "sumU " << sumU << " sumV: " << sumV << " pu " << ps.pU << " pV: " << ps.pV << endl; 

  assert(edges == 2*ps.L);
  return obj;
}
template<typename graph_type, typename vertex_data, typename edge_data>
void calc_stats2(const graph_type * gr, double&minU, double& maxU, int &moviewithoutedges, int& numedges, double &avgval, double& avgtime, double &minval, double &maxval, double&mintime, double&maxtime){
 for (int i=ps.M; i< ps.M+ps.N; i++){ 
    const vertex_data * data = &gr->vertex_data(i);
      if (min(data->pvec) < minU)
	 minU = min(data->pvec);
      if (max(data->pvec) > maxU)
	 maxU = max(data->pvec);
      if (gr->in_edge_ids(i).size() == 0)
	 moviewithoutedges++;
    foreach(edge_id_t iedgeid, gr->in_edge_ids(i)) {
	edge_data & data = (edge_data&)gr->edge_data(iedgeid);
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
    }
}
}

template<>
void calc_stats2<graph_type_mult_edge, vertex_data, multiple_edges>(const graph_type_mult_edge * gr, double&minU, double& maxU, int &moviewithoutedges, int& numedges, double &avgval, double& avgtime, double &minval, double &maxval,double&mintime, double&maxtime){
 for (int i=ps.M; i< ps.M+ps.N; i++){ 
    const vertex_data * data = &gr->vertex_data(i);
      if (min(data->pvec) < minU)
	 minU = min(data->pvec);
      if (max(data->pvec) > maxU)
	 maxU = max(data->pvec);
      if (gr->in_edge_ids(i).size() == 0)
	 moviewithoutedges++;
    foreach(edge_id_t iedgeid, gr->in_edge_ids(i)) {
      const multiple_edges & edges = gr->edge_data(iedgeid);
         //vertex_data * pdata = &gr->vertex_data(gr->source(iedgeid)); 
      for (int j=0; j< (int)edges.medges.size(); j++){     
	const edge_data_mcmc & data = edges.medges[j];
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
    }
 }
}
}
/// calc statistics about matrix/tensor and exit  
template<typename graph_type, typename vertex_data,typename edge_data>
void calc_stats(testtype type){
   const graph_type * gr = ps.g<graph_type>(type);

   if (gr->num_vertices() == 0){
     printf("%s is missing, skipping data\n", testtypename[type]);
     return;
   } 

   if (ps.tensor && type == TRAINING){
     int firsttimeused=-1;
     int lasttimeused=-1;
     for (int i=0; i<ps.K; i++){
       if (edges[i].size() > 0)
         firsttimeused = i;
     }
     for (int i=ps.K-1; i>=0; i--){
       if (edges[i].size() > 0)
         lasttimeused = i;
     }
     printf("Out of total %d time components, first used is %d, last used is %d\n", ps.K, firsttimeused, lasttimeused);
  }	
  double avgval=0, minval=1e100, maxval=-1e100;
  double avgtime=0, mintime=1e100, maxtime=-1e100;
  double minV=1e100, maxV=-1e100, minU=1e100, maxU=-1e100;
  int moviewithoutedges = 0;
  int userwithoutedges = 0;
  int numedges = 0;
  calc_stats2<graph_type, vertex_data, edge_data>(gr, minU, maxU, moviewithoutedges, numedges, avgval, avgtime, minval, maxval, mintime, maxtime);
 
  for (int i=0; i< ps.M; i++){ 
   const vertex_data * data = &gr->vertex_data(i);
   if (min(data->pvec) < minV)
     minV = min(data->pvec);
   if (max(data->pvec) > maxV)
     maxV = max(data->pvec);
	 
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
   case TRAINING: assert(numedges==ps.L); break;
   case VALIDATION: assert(numedges==ps.Le); break;
   case TEST: assert(numedges==ps.Lt); break;
 }
}

void predict_missing_value(const vertex_data& data, 
			   const vertex_data& pdata,
			   edge_data_mcmc& edge,
			   double &sq_err, int &e, int i);
void predict_missing_value(const vertex_data& data, 
			   const vertex_data& pdata,
			   const edge_data& edge,
			   double &sq_err, int &e, int i);
void predict_missing_value(const vertex_data_svdpp&data, 
			   const vertex_data_svdpp& pdata, 
		           edge_data& edge, 
			   double & sq_err, int&e, int i);
 
template<typename graph_type, typename vertex_data>
void calc_rmse_edge(edge_id_t iedgeid, const graph_type *_g, double & rmse, const vertex_data&data, const vertex_data&pdata, int&e, int i){
	   edge_data & edge = (edge_data&)_g->edge_data(iedgeid);
           double sq_err;
           predict_missing_value(data, pdata, edge, sq_err, e, i);
           rmse+= sq_err;
} 

template<>
void calc_rmse_edge<graph_type_mult_edge, vertex_data>(edge_id_t iedgeid, const graph_type_mult_edge *_g, double & rmse, const vertex_data & data, const vertex_data & pdata, int &e, int i){
         const multiple_edges & edges = _g->edge_data(iedgeid);
         for (int j=0; j< (int)edges.medges.size(); j++){       
           edge_data_mcmc & edge = (edge_data_mcmc&)edges.medges[j];
           double sq_err;
           predict_missing_value(data, pdata, (edge_data_mcmc&)edge, sq_err, e, i);
           rmse+= sq_err;
         }
} 




//calculate RMSE. This function is called only before and after grahplab is run.
//during run, agg_rmse_by_movie is called 0 which is much lighter function (only aggregate sums of squares)
template<typename graph_type, typename vertex_data>
double calc_rmse(const graph_type * _g, bool test, double & res){
     if (test && ps.Le == 0)
       return NAN;
     
     if (ps.algorithm == LANCZOS) //not implemented yet
       return NAN;
 
     res = 0;
     double RMSE = 0;
     int e = 0;
     for (int i=ps.M; i< ps.M+ps.N; i++){
       const vertex_data & data = ps.g<graph_type>(TRAINING)->vertex_data(i);
       foreach(edge_id_t iedgeid, _g->in_edge_ids(i)) {
         const vertex_data & pdata = ps.g<graph_type>(TRAINING)->vertex_data(_g->source(iedgeid)); 
         calc_rmse_edge<graph_type, vertex_data>(iedgeid, _g, RMSE, data, pdata, e, i);       
     }
   }
   res = RMSE;
   assert(e == (test?ps.Le:ps.L));
   return sqrt(RMSE/(double)e);

}
 

template<typename graph_type, typename vertex_data>
double calc_rmse_wrapper(const graph_type* _g, bool test, double & res){
   if (ps.algorithm == LANCZOS){
       res=-1; return -1; //not implemented yet
   }
   else return calc_rmse<graph_type, vertex_data>(_g, test, res);
}

  
// go over all movie edges and aggregate RMSE by summing up the squares, computing the
// mean and then sqrt()
template<typename graph_type, typename vertex_data>
double agg_rmse_by_movie(double & res){

  timer t; t.start();
  res = 0;
  double RMSE = 0;
  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    const vertex_data * data = &ps.g<graph_type>(TRAINING)->vertex_data(i);
    RMSE+= data->rmse;
  }
  res = RMSE;
  ps.counter[CALC_RMSE_Q] += t.current_time();
  return sqrt(RMSE/(double)ps.L);

}
// go over all user edges and aggregate RMSE by summing up the squares, computing the
// mean and then sqrt()
template<typename graph_type, typename vertex_data>
double agg_rmse_by_user(double & res){

  timer t; t.start();
  res = 0;
  double RMSE = 0;
  for (int i=0; i< ps.M; i++){ 
    const vertex_data * data = &ps.g<graph_type>(TRAINING)->vertex_data(i);
    RMSE+= data->rmse;
  }
  res = RMSE;
  ps.counter[CALC_RMSE_Q] += t.current_time();
  return sqrt(RMSE/(double)ps.L);

}

#include <graphlab/macros_undef.hpp>
#endif //_STATS_HPP
