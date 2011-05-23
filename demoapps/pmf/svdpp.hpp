#ifndef __SVD_HPP
#define __SVD_HPP

#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>


float itmBiasStep = 5e-3f;
float itmBiasReg = 1e-3f;
float usrBiasStep = 2e-4f;
float usrBiasReg = 5e-3f;
float usrFctrStep = 2e-2f;
float usrFctrReg = 2e-2f;
float itmFctrStep = 3e-3f;
float itmFctrReg = 1e-2f;
float itmFctr2Step = 1e-4f;
float itmFctr2Reg = 1e-2f;

extern string infile;
extern int iiter, L, Le;
extern bool ZERO;
extern timer gt;
extern graph_type validation_graph;
extern double globalMean[3];
extern bool debug;

double bestValidSqErr=DBL_MAX;
double stepSize=8e-3;
double regularization = 15e-3;

using namespace graphlab;
using namespace itpp;

void calc_user_moviebag(gl_types::iscope & scope, vertex_data & user, edge_list & outs){

    user.weight = zeros(D);
    foreach(graphlab::edge_id_t oedgeid, outs) {
       //edge_data & edge = scope.edge_data(oedgeid);
       const vertex_data  & movie = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
       user.weight += movie.weight;
    }
} 


void svd_init(){
   fprintf(stderr, "SVD++ %d factors (rate=%2.2e, reg=%2.2e)\n", D,stepSize,regularization);
   for (int i=0; i<M+N; i++){
       vertex_data & data = g->vertex_data(i);
       data.weight = debug ? itpp::ones(D) : itpp::randu(D);
   } 
}

// go over all edges and aggregate RMSE by summing up the squares, computing the
// mean and then sqrt()
double calc_rmse_q2(double & res){

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


//calculate RMSE. This function is called only before and after grahplab is run.
//during run, calc_rmse_q is called 0 which is much lighter function (only aggregate sums of squares)
double calc_svd_rmse(graph_type * _g, bool test, double & res){

     if (test && Le == 0)
       return NAN;
      
     
     res = 0;
     double sqErr =0;
     int nCases = 0;
     for (int i=0; i< M; i++){
       vertex_data & usr = g->vertex_data(i);
       int n = usr.num_edges; //+1.0 ? //regularization
       usr.weight = zeros(D);
       foreach(edge_id_t oedgeid, g->out_edge_ids(i)) {
         vertex_data & movie = g->vertex_data(g->target(oedgeid)); 
	 usr.weight += movie.weight;
        }
        float usrnorm = float(1.0/sqrt(n));
        usr.weight *= usrnorm;

       foreach(edge_id_t oedgeid, _g->out_edge_ids(i)){
         edge_data & item = _g->edge_data(oedgeid);
         vertex_data & movie = g->vertex_data(_g->target(oedgeid)); 
         float estScore = (float)globalMean[0];
         estScore += movie.pvec*(usr.weight+usr.pvec);
         estScore += movie.bias + usr.bias;
         estScore = min(estScore, maxval);
         estScore = max(estScore, minval);
         double err = item.weight - estScore;
         sqErr += err*err;
         nCases++;
       }
   }
   res = sqErr;
   assert(nCases == (test?Le:L));
   return sqrt(sqErr/(double)nCases);
}


void svd_post_iter(){
  printf("Entering last iter with %d\n", iiter);

  double res,res2;
  double rmse = calc_rmse_q2(res);
  printf("%g) Iter %s %d, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", gt.current_time(), "SVD", iiter,  rmse, calc_svd_rmse(&validation_graph, true, res2));

  itmFctrStep *= 0.9f;
  itmFctr2Step *= 0.9f;
  usrFctrStep *= 0.9f;
  itmBiasStep *= 0.9f;
  usrBiasStep *= 0.9f;

  iiter++;
}



/***
 * UPDATE FUNCTION
 */
void svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
 
  
  /* print statistics */
  if (debug&& (scope.vertex() == 0 || ((int)scope.vertex() == M-1) || ((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1) || ((int)scope.vertex() == 93712))){
    printf("SVDPP: entering %s node  %u \n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , user.pvec, D);
  }

  assert((int)scope.vertex() < M+N);

  user.rmse = 0;

  if (user.num_edges == 0){
    return; //if this user/movie have no ratings do nothing
  }


  edge_list outs = scope.out_edge_ids();
  edge_list ins = scope.in_edge_ids();
  timer t;

  int i=0;

  t.start(); 
  //USER NODES    
  if ((int)scope.vertex() < M){


    user.weight = zeros(D);
    
    foreach(graphlab::edge_id_t oedgeid, outs) {

      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid)); 
      user.weight += movie.weight; 
            
    }
   
   float usrNorm = float(1.0/sqrt(user.num_edges));
   user.weight *= usrNorm;

   vec step = zeros(D);
 

   foreach(graphlab::edge_id_t oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      float estScore = globalMean[0];
      estScore += movie.pvec*(user.pvec+user.weight);
      estScore += user.bias + movie.bias;
      estScore = min(estScore, maxval);
      estScore = max(estScore, minval);
      float err = edge.weight - estScore;
      user.rmse += err*err; 

     vec itmFctr = movie.pvec;
     vec usrFactor = user.pvec;
   
     movie.pvec += itmFctrStep*(err*(user.weight+usrFactor)-itmFctrReg*itmFctr);
     user.pvec += usrFctrStep*(err*itmFctr-usrFctrReg*usrFactor);
     step = sum(err*itmFctr);

     movie.bias += itmBiasStep*(err-itmBiasReg*movie.bias);
     user.bias += usrBiasStep*(err-usrBiasReg*user.bias);
   }

   step *= float(itmFctr2Step*usrNorm);

   double mult = itmFctr2Step*itmFctr2Reg;
   foreach(graphlab::edge_id_t oedgeid, outs){
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      movie.weight +=  step-mult*movie.weight;
   }


   counter[EDGE_TRAVERSAL] += t.current_time();

   if (scope.vertex() == M-1)
  	svd_post_iter();

  }

}

#include "graphlab/macros_undef.hpp"
#endif //__SVD_HPP
