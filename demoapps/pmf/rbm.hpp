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

float RandDouble()
{
    return rand() / static_cast<float>(RAND_MAX);
}

float RandInt(int i)
{
    return round(RandDouble() * (float)(i));
}
/*float rand01(){
     float r = (float)(rand() % 10000000);
     return r / 10000000; 
}*/

void setRand2(vec & a, int d, float c){
    for(int i = 0; i < d; ++i)
        a[i] = ((rand01() - 0.5) * c);
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

   rbm_user(const vertex_data_svdpp & vdata){
     h = (double*)&vdata.pvec[0];
     h0 = (double*)&vdata.weight[0];
     h1 = (double*)&vdata.weight[ac.D];
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

     rbm_movie(const vertex_data_svdpp& vdata){
     ni = (float*)&vdata.bias;
     bi = (double*)&vdata.pvec[0];
     w = (double*)&vdata.weight[0];
   }
};

/*inline float rand01(){
     float r = (float)(rand() % 10000000);
     return r / 10000000;
}*/
float rbm_predict(const vertex_data_svdpp& user, 
                const vertex_data_svdpp& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
//inline float SVDManager::rating(const ItemRating& itemRating, unsigned int user)
	  //unsigned int item = itemRating.item;
	  rbm_user usr(user);
    rbm_movie mov(movie);
    float ret = 0;
    double nn = 0;
    for(int r = 0; r < ac.rbm_bins; ++r){               
        //float zz = exp(bi[item][r] + dot(h[user], w[r][item]));
        double zz = exp(mov.bi[r] + dot(usr.h, &mov.w[r*ac.D]));
        ret += zz * (float)(r);
        assert(!std::isnan(ret));
        nn += zz;
    }
    if (std::fabs(nn) < 1e-32)
       ret = ac.minval;
    else {
    ret /= nn;
    if(ret < ac.minval) ret = ac.minval;
    else if(ret > ac.maxval) ret = ac.maxval;
    }
    assert(!std::isnan(ret));
    prediction = ret * ac.rbm_scaling;
    assert(!std::isnan(prediction));
    return pow(prediction - rating,2);
}

//inline int rating1(const ItemRating& itemRating, unsigned int user)
float rbm_predict1(const vertex_data_svdpp& user, 
                const vertex_data_svdpp& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
//inline float SVDManager::rating(const ItemRating& itemRating, unsigned int user)
	  //unsigned int item = itemRating.item;
	  rbm_user usr(user);
    rbm_movie mov(movie);
  	//unsigned int item = itemRating.item;
    vec zz = zeros(ac.rbm_bins);
    float szz = 0;
    for(int r = 0; r < ac.rbm_bins; ++r){
        //zz[r] = exp(bi[item][r] + dot(h0[user], w[r][item]));
        zz[r] = exp(mov.bi[r] + dot(usr.h0, &mov.w[r*ac.D]));
        szz += zz[r];
    }
    float rd = rand01() * szz;
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
/*
inline float SVDManager::RMSE(){
    float totalMse = 0;
	int currentRatingIdx = 0;
	float estScore = 0;

    for (unsigned int user=0; user<ValidationMetaData.nUsers; user++)
	{
		for (int i=0; i<RATINGS_PER_USER_VALIDATION; i++)
		{
			ItemRating currentItemRating = pItemRatings_validation[currentRatingIdx];
            estScore = rating(currentItemRating,user);
			float err = float(currentItemRating.score) - estScore * 10;
			totalMse += err*err;
			currentRatingIdx++;
		}
	}
    return sqrt(totalMse / ValidationMetaData.nRecords);
}
*/

double calc_rbm_rmse(const graph_type_svdpp * _g, bool test, double & res){

     graph_type_svdpp * g = (graph_type_svdpp*)ps.g<graph_type_svdpp>(TRAINING);
     if (test && ps.Le == 0)
       return NAN;
      
     
     res = 0;
     double sqErr =0;
     int nCases = 0;

     for (int i=0; i< ps.M; i++){
       vertex_data_svdpp & usr = (vertex_data_svdpp&)g->vertex_data(i);
       int n = usr.num_edges; //+1.0 ? //regularization
       if (n == 0){
         nCases += _g->out_edge_ids(i).size();
       }
       else {
         foreach(graphlab::edge_id_t oedgeid, _g->out_edge_ids(i)){
           const edge_data & item = _g->edge_data(oedgeid);
           const vertex_data_svdpp & movie = g->vertex_data(_g->target(oedgeid)); 
           float estScore;
           sqErr += rbm_predict(usr, movie, NULL, NULL, item.weight, estScore);
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
  printf("Entering last iter with %d\n", ps.iiter);

  double res,res2;
  double rmse = agg_rmse_by_user<graph_type_svdpp, vertex_data_svdpp>(res);
  printf("%g) Iter %s %d, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", ps.gt.current_time(), "SVD", ps.iiter,  rmse, calc_rbm_rmse(ps.g<graph_type_svdpp>(VALIDATION), true, res2));

  if (ac.calc_ap){
     logstream(LOG_INFO)<<"AP@3 for training: " << calc_ap<graph_type_svdpp,vertex_data_svdpp,edge_data>(ps.g<graph_type_svdpp>(TRAINING)) << " AP@3 for validation: " << calc_ap<graph_type_svdpp,vertex_data_svdpp,edge_data>(ps.g<graph_type_svdpp>(VALIDATION)) << std::endl;
  }

  ac.rbm_alpha *= ac.rbm_mult_step_dec;
  ps.iiter++;
}


/***
 * UPDATE FUNCTION
 */
void rbm_update_function(gl_types_svdpp::iscope &scope, 
			 gl_types_svdpp::icallback &scheduler) {
   
   /* GET current vertex data */
  vertex_data_svdpp& user = scope.vertex_data();
  //h[user_idx] = vector<float>(dim,0);
  //h0[user_idx] = vector<float>(dim,0);
   user.pvec = zeros(ac.D);
  user.weight = zeros(2*ac.D);
  rbm_user usr(user); 
  int id = scope.vertex(); 

  /* print statistics */
  if (ac.debug&& (id == 0 || id == ps.M-1)){
    printf("RBM: entering %s node  %u \n", "user", id);   
  }

  assert(id < ps.M);
  user.rmse = 0;

  gl_types_svdpp::edge_list outs = scope.out_edge_ids();
  if (outs.size() == 0)
     return;

  timer t;
  t.start(); 
           
  //unsigned int userRatings = pUsersData[user_idx].ratings;
	unsigned int userRatings = outs.size(); 
  //vector<int> v1(userRatings, 0); 
  vec v1 = zeros(userRatings); 
  //for (int i = 0; i < userRatings; i++){
  foreach(graphlab::edge_id_t oedgeid, outs){
    //ItemRating currentItemRating = pItemRatings_training[currentRatingIdx+i];
    //unsigned int item_idx = currentItemRating.item;
    //int r = currentItemRating.score / ac.rbm_scaling;
    vertex_data_svdpp & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
    rbm_movie mov(movie);
    edge_data & edge = scope.edge_data(oedgeid);
    int r = edge.weight / ac.rbm_scaling;
            
    //for (int k=0; k < dim; k++){
    for(int k=0; k < ac.D; k++){
      //h[user_idx][k] += w[r][item_idx][k]; 
      usr.h[k] += mov.w[ac.D*r + k];
    }
  } 
    //for (int k=0; k < dim; k++){
    //  h[user_idx][k] = sigmoid(h[user_idx][k]);
    //  if(rand01() < h[user_idx][k])  h0[user_idx][k] = 1;
    //    else h0[user_idx][k] = 0;
    //}
    for(int k=0; k < ac.D; k++){
       usr.h[k] = sigmoid(usr.h[k]);
       if (rand01() < usr.h[k]) 
          usr.h0[k] = 1;
       else usr.h0[k] = 0;
    }

   //for (int i = 0; i < userRatings; i++){
    int i = 0;
    float prediction;
    foreach(graphlab::edge_id_t oedgeid, outs){
      //ItemRating currentItemRating = pItemRatings_training[currentRatingIdx+i];
      //int vi = rating1(currentItemRating, user_idx);
      vertex_data_svdpp & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      edge_data & edge  = scope.edge_data(oedgeid);
      rbm_predict1(user, movie, &edge, NULL, edge.weight, prediction);
      int vi = prediction / ac.rbm_scaling;
      v1[i] = vi;
      i++;
    }

    //h1[user_idx] = vector<float>(dim,0);
    //usr.h1 = zeros(ac.D);
    i = 0;
    //for (int i = 0; i < userRatings; i++){
    foreach(graphlab::edge_id_t oedgeid, outs){
      //          ItemRating currentItemRating = pItemRatings_training[currentRatingIdx+i];
      //          unsigned int item_idx = currentItemRating.item;
      vertex_data_svdpp & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      rbm_movie mov(movie);
      //edge_data & edge = scope.edge_data(oedgeid);
      int r = v1[i];
      for (int k=0; k< ac.D;k++){
        //h1[user_idx][k] += w[r][item_idx][k];
        usr.h1[k] += mov.w[r*ac.D+k];
      }
      i++;
   }
  
   for (int k=0; k < ac.D; k++){
     //h1[user_idx][k] = sigmoid(h1[user_idx][k]);
     usr.h1[k] = sigmoid(usr.h1[k]);
     //if(rand01() < h1[user_idx][k])  h1[user_idx][k] = 1;
     //  else h1[user_idx][k] = 0;
     if (rand01() < usr.h1[k]) 
       usr.h1[k] = 1;
     else usr.h1[k] = 0;
  }

  //for (int i = 0; i < userRatings; i++){
    i = 0;
    foreach(graphlab::edge_id_t oedgeid, outs){
      //          ItemRating currentItemRating = pItemRatings_training[currentRatingIdx+i];
      //          unsigned int item_idx = currentItemRating.item;
      vertex_data_svdpp & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      rbm_movie mov(movie);
      edge_data & edge = scope.edge_data(oedgeid);
 
      //float pui = rating(currentItemRating, user_idx);
      float prediction;
      rbm_predict(user, movie, &edge, NULL, edge.weight, prediction);
      //float rui = currentItemRating.score / ac.rbm_scaling;
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

   if (scope.vertex() == (uint)(ps.M-1))
  	rbm_post_iter();
 
}

 
//SVDManager::SVDManager(int k){
void rbm_init(){
	  srand((unsigned)time(NULL)); 
    //int NUSER = TrainingMetaData.trainingTotalUsers;
    //int NITEM = TrainingMetaData.trainingTotalItems;
    //h = vector< vector<float> >(NUSER);
    //bi = vector< vector<float> >(NITEM);
    //h0 = vector< vector<float> >(NUSER);
    //h1 = vector< vector<float> >(NUSER);
    //w = vector< vector< vector<float> > >(ac.rbm_bins);
    graph_type_svdpp * training = (graph_type_svdpp*)ps.g<graph_type_svdpp>(TRAINING);
    for(int i = 0; i < ps.N; ++i){
        vertex_data_svdpp & movie = training->vertex_data(ps.M+i);
        movie.pvec = zeros(ac.rbm_bins);
        movie.weight = zeros(ac.D*ac.rbm_bins);
        movie.bias = 0;
    }
     
    //v1 = vector<uint8>(TrainingMetaData.trainingTotalRatings);
    //vector<float> ni(NITEM,0);
    int currentRatingIdx = 0;
    //for (unsigned int user_idx=0; user_idx<TrainingMetaData.nUsers; user_idx++) {
    for (int i=0; i< ps.M; i++){
        //unsigned int userRatings = pUsersData[user_idx].ratings;
        vertex_data_svdpp & user = training->vertex_data(i);
        user.pvec = zeros(ac.D);
        user.weight = zeros(2*ac.D); 
         gl_types_svdpp::edge_list outs = training->out_edge_ids(i);
        //for (unsigned int i = 0; i < userRatings; i++){
         foreach(graphlab::edge_id_t oedgeid, outs){
           // ItemRating currentItemRating = pItemRatings_training[currentRatingIdx];
           // unsigned int item_idx = currentItemRating.item;
            vertex_data_svdpp & movie = training->vertex_data(training->target(oedgeid));
            rbm_movie mov(movie);
            edge_data & edge = training->edge_data(oedgeid);
            int r = edge.weight/ac.rbm_scaling;
            mov.bi[r]++;
            (*mov.ni)++;
        }
        currentRatingIdx++;
    }
    for(int i = 0; i < ps.N; ++i){
        vertex_data_svdpp & movie = training->vertex_data(ps.M+i);
        rbm_movie mov(movie);
        setRand2(movie.weight, ac.D*ac.rbm_bins, 0.001);
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
