#ifndef __SVD_HPP
#define __SVD_HPP

#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>

float LRATE1 = 0.003;               // Learning rate parameter
float LRATE2 = 0.003;               // Learning rate parameter for biases
float LRATE3 = 0.001;               // Learning rate parameter for weights
float LAMBDA1 = 0.015;               // Regularization parameter used to minimize over-fitting
float LAMBDA2 = 0.005;               // Biases regularization
float LAMBDA3 = 0.015;               // Regularization parameter for weights

extern string infile;
extern int iiter, L, Le;
extern bool ZERO;
extern timer gt;
extern graph_type validation_graph;

using namespace graphlab;
using namespace itpp;

void calc_users_moviebag(gl_types::iscope & scope, vertex_data & user, edge_list & outs){

}
void calc_user_moviebag(gl_types::iscope & scope, vertex_data & user, edge_list & outs){

    user.weight = zeros(D);
    foreach(graphlab::edge_id_t oedgeid, outs) {
       //edge_data & edge = scope.edge_data(oedgeid);
       const vertex_data  & movie = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
       user.weight += movie.weight;
    }
} 


float predict_svd_rating(const vertex_data & user, const vertex_data & movie, float ei){

  float dsum = itpp::dot(movie.pvec, user.pvec + ei*user.weight);
  dsum += user.bias + movie.bias;

  //truncate values
  if (infile == "kddcup" || infile == "kddcup2"){
    if (dsum > 100)
	dsum = 100;
    else if (dsum < 0)
	dsum = 0;
  }

  return dsum;
}
//calculate RMSE. This function is called only before and after grahplab is run.
//during run, calc_rmse_q is called 0 which is much lighter function (only aggregate sums of squares)
double calc_svd_rmse(graph_type * _g, bool test, double & res){

     if (test && Le == 0)
       return NAN;
      
     
     res = 0;
     double RMSE = 0;
     int e = 0;
     for (int i=0; i< M; i++){
       vertex_data & data = g->vertex_data(i);
       float ei = 1.0 / sqrtf(data.num_edges + 1.0); //regularization
       foreach(edge_id_t oedgeid, _g->out_edge_ids(i)) {
         vertex_data & pdata = g->vertex_data(_g->target(oedgeid)); 
            
#ifndef GL_NO_MULT_EDGES
         multiple_edges & edges = _g->edge_data(oedgeid);
         for (int j=0; j< (int)edges.medges.size(); j++){       
           edge_data & edge = edges.medges[j];
#else
	   edge_data & edge = _g->edge_data(oedgeid);
#endif
           if (!ZERO)
           	assert(edge.weight != 0);

           float p = predict_svd_rating(data, pdata, ei);
           float err = edge.weight -  p;
           if (debug && (i== M || i == M+N-1) && (e == 0 || e == (test?Le:L)))
             cout<<"RMSE:"<<i <<"u1"<< data.pvec << " v1 "<< pdata.pvec<<endl; 

           RMSE+= (err*err);
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

void svd_post_iter(){
  printf("Entering last iter with %d\n", iiter);

  double res,res2;
  double rmse = calc_rmse_q(res);
  printf("%g) Iter %s %d, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", gt.current_time(), "SVD", iiter,  rmse, calc_svd_rmse(&validation_graph, true, res2));
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

    calc_user_moviebag(scope, user, outs);

    foreach(graphlab::edge_id_t oedgeid, outs) {

#ifndef GL_NO_MULT_EDGES
      multiple_edges &medges =scope.edge_data(oedgeid);
#else
      edge_data & edge = scope.edge_data(oedgeid);
#endif
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid)); 
#ifndef GL_NO_MULT_EDGES	   
      for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
#endif
        //go over each rating of a movie and put the movie vector into the matrix Q
        //and vector vals
        float ei = 1.0 / sqrtf(outs.size() + 1.0); //regularization
        float p = predict_svd_rating(user,movie,ei);    
        float err = edge.weight -  p;
        user.rmse += err*err;

        float cf_bias = user.bias;
        float mf_bias = movie.bias;

        user.bias += (float)(LRATE2 * (err-LAMBDA2 * cf_bias));
        movie.bias += (float)(LRATE2 * (err-LAMBDA2 * mf_bias));

	//cache off old feature values
        vec cf = user.pvec;
        vec mf = movie.pvec;
        vec wf = movie.weight;


        user.pvec += (LRATE1 * (err*mf - LAMBDA1*cf));
        movie.pvec += (LRATE1 * (err*(itpp::dot(cf+ei,user.weight)) - LAMBDA1*mf));
        movie.weight += LRATE3*(err*ei*mf-LAMBDA3*wf);
        user.weight += movie.weight - wf;
	i++;
#ifndef GL_NO_MULT_EDGES
      }   
#endif
           
            
    }
   assert(i == user.num_edges);
   counter[EDGE_TRAVERSAL] += t.current_time();

   if (scope.vertex() == 0)
  	svd_post_iter();

  }

}

#include "graphlab/macros_undef.hpp"
#endif //__SVD_HPP
