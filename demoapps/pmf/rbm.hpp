#ifndef __RBM_HPP
#define __RBM_HPP
 #include "graphlab.hpp"
#include "pmf.h"
#include "stats.hpp"
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <boost/random.hpp>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <limits.h>
#include <time.h>
#include <graphlab/macros_def.hpp>

extern advanced_config ac;
extern problem_setup ps;

using namespace graphlab;
using namespace std;

void setRand2(double * a, int d, float c){
    for(int i = 0; i < d; ++i)
        a[i] = ((graphlab::random::rand01() - 0.5) * c);
}

float dot(double * a, double * b){
  float ret = 0;
  for(int i = 0; i < ac.D; ++i)
    ret += a[i] * b[i];
  return ret;
}

/*
 * h = pvec = D * DOUBLE
 * h0 = weight = D * DOUBLE
 * h1 = weight+ac.D = D * DOUBLE
 */
struct rbm_user{
   double * h;
   double * h0;
   double * h1;
   int num_edges; 

   rbm_user(const vertex_data & vdata){
     h = (double*)&vdata.pvec[0];
     h0 = h + ac.D;
     h1 = h0 + ac.D;
     num_edges = vdata.num_edges;
   }
   
   rbm_user & operator=(vertex_data & data){
     h = &data.pvec[0];
     h0 = h + ac.D;
     h1 = h0 + ac.D;
     num_edges = data.num_edges;
     return * this;
    }   
};


/**
 * ni = bias = DOUBLE
 * bi = pvec = rbm_bins * DOUBLE 
 * w = weight = rbm_bins * D * Double
 *
 */
struct rbm_movie{
   double * bi;
   float * ni;
   double * w;

   rbm_movie(const vertex_data& vdata){
     ni = (float*)&vdata.bias;
     bi = (double*)&vdata.pvec[0];
     w = bi + ac.D;
   }

   rbm_movie & operator=(vertex_data & data){
     ni = (float*)&data.bias;
     bi = (double*)&data.pvec[0];
     w = bi + ac.D;
     return * this;
   }
};

float predict(const rbm_user & usr, 
                const rbm_movie & mov, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
	  
    float ret = 0;
    double nn = 0;
    for(int r = 0; r < ac.rbm_bins; ++r){               
        double zz = exp(mov.bi[r] + dot(usr.h, &mov.w[r*ac.D]));
        ret += zz * (float)(r);
        assert(!std::isnan(ret));
        nn += zz;
    }
    assert (std::fabs(nn) > 1e-32);
    ret /= nn;
    if(ret < ac.minval) ret = ac.minval;
    else if(ret > ac.maxval) ret = ac.maxval;
    assert(!std::isnan(ret));
    prediction = ret * ac.rbm_scaling;
    assert(!std::isnan(prediction));
    return pow(prediction - rating,2);
}

float predict1(const rbm_user & usr, 
                const rbm_movie & mov, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
	  
    vec zz = zeros(ac.rbm_bins);
    float szz = 0;
    for(int r = 0; r < ac.rbm_bins; ++r){
        zz[r] = exp(mov.bi[r] + dot(usr.h0, &mov.w[r*ac.D]));
        szz += zz[r];
    }
    float rd = graphlab::random::rand01() * szz;
    szz = 0;
    int ret = 0;
    for(int r = 0; r < ac.rbm_bins; ++r){
        szz += zz[r];
        if(rd < szz){ 
            ret = r;
            break;
        }
    }
    prediction = ret * ac.rbm_scaling;
    assert(!std::isnan(prediction));
    return pow(prediction - rating, 2);
}

inline float sigmoid(float x){
    return 1 / (1 + exp(-1 * x));
}

double calc_rbm_rmse(const graph_type * _g, bool test, double & res){

  graph_type * g = (graph_type*)ps.g<graph_type>(TRAINING);
  if (test && ps.Le == 0)
       return NAN;
      
     
     res = 0;
     double sqErr =0;
     int nCases = 0;

     for (int i=0; i< ps.M; i++){
       rbm_user usr = g->vertex_data(i);
       int n = usr.num_edges; //+1.0 ? //regularization
       if (n == 0){
         nCases += _g->out_edge_ids(i).size();
       }
       else {
         foreach(graphlab::edge_id_t oedgeid, _g->out_edge_ids(i)){
           const edge_data & item = _g->edge_data(oedgeid);
           const rbm_movie movie = g->vertex_data(_g->target(oedgeid)); 
           float estScore;
           sqErr += predict(usr, movie, NULL, NULL, item.weight, estScore);
           nCases++;
         }
       }
   }
   assert(!std::isnan(sqErr));
   res = sqErr;
   ASSERT_EQ(nCases , (test?ps.Le:ps.L));
   return sqrt(sqErr/(double)nCases);
}



void rbm_post_iter(){
  //printf("Entering last iter with %d\n", ps.iiter);
  double res,res2;
  double training_rmse = agg_rmse_by_user<graph_type, vertex_data>(res);
  double validation_rmse = calc_rbm_rmse(ps.g<graph_type>(VALIDATION), true, res2);
  printf(ac.printhighprecision ? 
        "%g) Iter %s %d  TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n":
        "%g) Iter %s %d  TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n",
         ps.gt.current_time(), runmodesname[ps.algorithm], ps.iiter,  training_rmse, validation_rmse);

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

  ac.rbm_alpha *= ac.rbm_mult_step_dec;
  ps.iiter++;
}


/***
 * UPDATE FUNCTION
 */
void rbm_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
   
   /* GET current vertex data */
  vertex_data user = scope.vertex_data();
  user.pvec = zeros(ac.D*3);
  rbm_user usr(user); 
  int id = scope.vertex(); 

  /* print statistics */
  if (ps.to_print(id)){
    printf("RBM: entering user node  %u\n", id);   
  }

  user.rmse = 0;

  gl_types::edge_list outs = scope.out_edge_ids();
  if (outs.size() == 0){
     if (id == ps.M -1)
       rbm_post_iter();
     return;
  }

  timer t;
  t.start(); 
           
	unsigned int userRatings = outs.size(); 
  vec v1 = zeros(userRatings); 
  foreach(graphlab::edge_id_t oedgeid, outs){
    rbm_movie mov = scope.neighbor_vertex_data(scope.target(oedgeid));
    edge_data & edge = scope.edge_data(oedgeid);
    int r = edge.weight / ac.rbm_scaling;
            
    for(int k=0; k < ac.D; k++){
      usr.h[k] += mov.w[ac.D*r + k];
    }
  } 
    for(int k=0; k < ac.D; k++){
       usr.h[k] = sigmoid(usr.h[k]);
       if (graphlab::random::rand01() < usr.h[k]) 
          usr.h0[k] = 1;
       else usr.h0[k] = 0;
    }

    int i = 0;
    float prediction;
    foreach(graphlab::edge_id_t oedgeid, outs){
      rbm_movie movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      edge_data & edge  = scope.edge_data(oedgeid);
      predict1(user, movie, &edge, NULL, edge.weight, prediction);
      int vi = prediction / ac.rbm_scaling;
      v1[i] = vi;
      i++;
    }

    i = 0;
    foreach(graphlab::edge_id_t oedgeid, outs){
      rbm_movie mov = scope.neighbor_vertex_data(scope.target(oedgeid));
      int r = v1[i];
      for (int k=0; k< ac.D;k++){
        usr.h1[k] += mov.w[r*ac.D+k];
      }
      i++;
   }
  
   for (int k=0; k < ac.D; k++){
     usr.h1[k] = sigmoid(usr.h1[k]);
     if (graphlab::random::rand01() < usr.h1[k]) 
       usr.h1[k] = 1;
     else usr.h1[k] = 0;
  }

    i = 0;
    foreach(graphlab::edge_id_t oedgeid, outs){
      rbm_movie mov = scope.neighbor_vertex_data(scope.target(oedgeid));
      edge_data & edge = scope.edge_data(oedgeid);
      float prediction;
      predict(user, mov, &edge, NULL, edge.weight, prediction);
      double pui = prediction / ac.rbm_scaling;
      double rui = edge.weight / ac.rbm_scaling;
      user.rmse += (pui - rui) * (pui - rui);
      //nn += 1.0;
      int vi0 = (int)(rui);
      int vi1 = v1[i];
      for (int k = 0; k < ac.D; k++){
        mov.w[ac.D*vi0+k] += ac.rbm_alpha * (usr.h0[k] - ac.rbm_beta * mov.w[vi0*ac.D+k]);
        mov.w[ac.D*vi1+k] -= ac.rbm_alpha * (usr.h1[k] + ac.rbm_beta * mov.w[vi1*ac.D+k]);
     }
     i++; 
   }    
   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (id == ps.M-1)
  	rbm_post_iter();
 
}

 
void rbm_init(){
    graph_type * training = (graph_type*)ps.g<graph_type>(TRAINING);
    for(int i = 0; i < ps.N; ++i){
        vertex_data movie(training->vertex_data(ps.M+i));
        movie.pvec = zeros(ac.rbm_bins + ac.D * ac.rbm_bins);
        movie.bias = 0;
    }
    int currentRatingIdx = 0;
    for (int i=0; i< ps.M; i++){
        vertex_data user = training->vertex_data(i);
        user.pvec = zeros(ac.D*3);
         gl_types::edge_list outs = training->out_edge_ids(i);
         foreach(graphlab::edge_id_t oedgeid, outs){
            rbm_movie mov = training->vertex_data(training->target(oedgeid));
            edge_data & edge = training->edge_data(oedgeid);
            int r = edge.weight/ac.rbm_scaling;
            mov.bi[r]++;
            (*mov.ni)++;
        }
        currentRatingIdx++;
    }
    for(int i = 0; i < ps.N; ++i){
        rbm_movie mov = training->vertex_data(ps.M+i);
        setRand2(mov.w, ac.D*ac.rbm_bins, 0.001);
        if((*mov.ni) == 0) 
           continue;
        for(int r = 0; r < ac.rbm_bins; ++r){
            mov.bi[r] /= (*mov.ni);
            mov.bi[r] = log(1E-9 + mov.bi[r]);
        }
    }


    logstream(LOG_INFO) << "RBM initialization ok" << endl;

}

#include <graphlab/macros_undef.hpp>
#endif //__RBM_HPP
