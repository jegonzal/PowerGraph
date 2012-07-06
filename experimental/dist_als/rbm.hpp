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
     w = bi + ac.rbm_bins;
   }

   rbm_movie & operator=(vertex_data & data){
     ni = (float*)&data.bias;
     bi = (double*)&data.pvec[0];
     w = bi + ac.rbm_bins;
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

float rbm_predict(const vertex_data & usr, 
                const vertex_data & mov, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
  return predict(rbm_user((vertex_data&)usr), rbm_movie((vertex_data&)mov), edge, nothing, rating, prediction);
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

void rbm_post_iter(){

  post_iter_stats<graph_type>();
  ac.rbm_alpha *= ac.rbm_mult_step_dec;
}


/***
 * UPDATE FUNCTION
 */
void rbm_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
   
   /* GET current vertex data */
  vertex_data & user = scope.vertex_data();
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
    assert(r < ac.rbm_bins);      
    for(int k=0; k < ac.D; k++){
      usr.h[k] += mov.w[ac.D*r + k];
      assert(!std::isnan(usr.h[k]));
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
        assert(!std::isnan(mov.w[ac.D*vi0+k]));
        mov.w[ac.D*vi1+k] -= ac.rbm_alpha * (usr.h1[k] + ac.rbm_beta * mov.w[vi1*ac.D+k]);
        assert(!std::isnan(mov.w[ac.D*vi1+k]));
     }
     i++; 
   }    
   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (id == ps.M-1)
  	rbm_post_iter();
 
}

 
void rbm_init(){
    graph_type * training = (graph_type*)ps.g<graph_type>(TRAINING);
#pragma omp parallel for
    for(int i = 0; i < ps.N; ++i){
        vertex_data & movie = training->vertex_data(ps.M+i);
        movie.pvec = zeros(ac.rbm_bins + ac.D * ac.rbm_bins);
        movie.bias = 0;
    }
    int currentRatingIdx = 0;
    for (int i=0; i< ps.M; i++){
        vertex_data &user = training->vertex_data(i);
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
