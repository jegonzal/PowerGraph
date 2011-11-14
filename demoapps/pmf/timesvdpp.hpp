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
 Koren time-SVD++ is described in the paper:
13) Yehuda Koren. 2009. Collaborative filtering with temporal dynamics. In Proceedings of the 15th ACM SIGKDD international conference on Knowledge discovery and data mining (KDD '09). ACM, New York, NY, USA, 447-456. DOI=10.1145/1557019.1557072 
* Based on code by Yao Wu, Chinese National Academy of Science, Beijing
* Parallelization of the code is done by Danny Bickson, CMU
 */


#ifndef __TIME_SVDPP_HPP
#define __TIME_SVDPP_HPP

#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>
#include "pmf.h"

extern advanced_config ac;
extern problem_setup ps;


extern advanced_config ac;
extern problem_setup ps;
double lrate = 0.0001;
double lrate2 = 0.00005;
double timesvdpp_beta = 0.00001; 
double garma = 0.0001; 
double garma2 = 0.001; 


using namespace graphlab;

void time_svd_post_iter();


double RandDouble1()
{
    return rand() / static_cast<double>(RAND_MAX);
}

double RandInt1(int i)
{
    return ::round(RandDouble1() * (double)(i));
}

double rand01(){
     double r = (double)(rand() % 10000000);
     return r / 10000000;
}


struct time_svdpp_usr{
   float * bu;
   double * p;
   double * pu;
   double * x;
   double * ptemp;
   int * tu;

   time_svdpp_usr(vertex_data_svdpp & vdata){
      bu = &vdata.bias;
      p = (double*)data(vdata.pvec);
      pu = p+ac.K;
      x = (double*)data(vdata.weight);
      ptemp = x+ac.K;
      tu = &vdata.num_edges;
   }

};

struct time_svdpp_movie{
   float * bi;
   double * q;
   double * y;

   time_svdpp_movie(vertex_data_svdpp& vdata){
     bi = &vdata.bias;
     q = (double*)data(vdata.pvec);
     y = (double*)data(vdata.weight);
   }
};

struct time_svdpp_time{
   float * bt;
   double * z;
   double * pt;

   time_svdpp_time(vertex_data_svdpp& vdata){
     bt = &vdata.bias;
     z = (double*)data(vdata.pvec);
     pt = (double*)data(vdata.weight);
   }
};

template<typename graph_type>
void init_time_svdpp(graph_type* _g){
  assert(false);
}



template<>
void init_time_svdpp<graph_type_svdpp>(graph_type_svdpp *_g){
   fprintf(stderr, "time-SVD++ %d factors\n", ac.D);

   int k = ac.D;

	srand((unsigned)time(NULL));
/*	mu = TrainingMetaData.totalMeanScore;
	dim = k;
	//ps.K = 7000;
	ps.M = TrainingMetaData.trainingTotalUsers;
	ps.N = TrainingMetaData.trainingTotalItems;
	bu = new double[ps.M];
	bi = new double[ps.N];
	bt = new double[ps.K];
	p = new double*[ps.M];
	pu = new double*[ps.M];
	q = new double*[ps.N];
	y = new double*[ps.N]; // impact feedback
	x = new double*[ps.M];
	z = new double*[ps.K];
	pt = new double*[ps.K];
	ptemp = new double*[ps.M];
*/ 
    //tu = new int[ps.M];
    
       // int currentRatingIdx_test = 0; 
	for (int u = 0; u < ps.M; u++) {
              vertex_data_svdpp & data = _g->vertex_data(u);
              data.pvec = zeros(2*k);
              data.weight = zeros(2*k);
              time_svdpp_usr usr(data);
/*		bu[u] = 0;
		p[u] = new double[k];
		pu[u] = new double[k];
		x[u] = new double[k];
		ptemp[u] = new double[k];
		for (int m = 0; m < k; m++){
			p[u][m] = 0.01 * rand01() / (double) (k);
			pu[u][m] = 0.001 * rand01() / (double) (k);
			x[u][m] = 0.001 * rand01() / (double) (k);
			ptemp[u][m] = p[u][m];
		}
*/
             *usr.bu = 0;
             //*usr.tu = 0;
             for (int m=0; m< k; m++){
                usr.p[m] = 0.01*rand01() / (double) (k);
	        usr.pu[m] = 0.001 * rand01() / (double) (k);
	        usr.x[m] = 0.001 * rand01() / (double) (k);
	        usr.ptemp[m] = usr.p[m];
	
             }
         
  
        //int userRatings = pUsersData[u].ratings;
        //currentRatingIdx_test += userRatings;
        //int day = pItemRatings_training[currentRatingIdx_test-1].day;
            //foreach(edge_id_t oedgeid, _g->out_edge_ids(u)){
	      // tu[u] = max<int>(tu[u], day);
	      // usr.tu = max<int>(usr.tu, _g->edge_data(oedgeid).time);
           // }
    }

	for (int i = ps.M; i < ps.N+ps.M; i++) {
                vertex_data_svdpp & data = _g->vertex_data(i);
                data.pvec = zeros(k);
                data.weight = zeros(k);
                time_svdpp_movie movie(data);
		/*bi[i] = 0;
		q[i] = new double[k];
		y[i] = new double[k];
		for (int m = 0; m < k; m++){
			q[i][m] = 0.01 * rand01() / (double) (k);
			y[i][m] = 0.001 * rand01() / (double) (k);
		}*/
		*movie.bi = 0;
		for (int m = 0; m < k; m++){
			movie.q[m] = 0.01 * rand01() / (double) (k);
			movie.y[m] = 0.001 * rand01() / (double) (k);
		}
		
	}
        for (int i = ps.M+ps.N; i < ps.M+ps.N+ps.K; i++) {
                vertex_data_svdpp & data = _g->vertex_data(i);
                data.pvec = zeros(k);
 		data.weight = zeros(k);
		time_svdpp_time timenode(data);

	/*	bt[i] = 0;
		z[i] = new double[k];
		pt[i] = new double[k];
		for (int m = 0; m < k; m++){
			z[i][m] = 0.001 * rand01() / (double) (k);
			pt[i][m] = 0.001 * rand01() / (double) (k);
        	}
        */
		*timenode.bt = 0;
              	for (int m = 0; m < k; m++){
			timenode.z[m] = 0.001 * rand01() / (double) (k);
			timenode.pt[m] = 0.001 * rand01() / (double) (k);
        	}
 
	}

}

double predict(const time_svdpp_movie & mov, const time_svdpp_usr & usr, float rating, int time, double & estScore);


/***
 * UPDATE FUNCTION
 */
void time_svd_plus_plus_update_function(gl_types_svdpp::iscope &scope, 
			 gl_types_svdpp::icallback &scheduler) {
    
  /* GET current vertex data */
  vertex_data_svdpp& user = scope.vertex_data();
  time_svdpp_usr usr(user); 
  int id = scope.vertex(); 
  /* print statistics */
  if (ac.debug&& (id == 0 || id == ps.M-1)){
    printf("SVD++: entering %s node  %u \n", "user", id);   
  }

  assert(id < ps.M);

  user.rmse = 0;

  gl_types_svdpp::edge_list outs = scope.out_edge_ids();
  gl_types_svdpp::edge_list ins = scope.in_edge_ids();
  
  if (outs.size() == 0)
     return;

  timer t;
  t.start(); 
 

/*		Start_t = time(NULL);
		double prmse = 0;
		int currentRatingIdx = 0;
        int currentRatingIdx_y = 0; //for y iteration
		int currentRatingIdx_val = 0;
		int currentRatingIdx_val1 = 0;
        int currentRatingIdx_test =0 ;
        double eui = 0;
*/  
//      for (unsigned int user_idx=0; user_idx<TrainingMetaData.nUsers; user_idx++) {



	unsigned int userRatings = outs.size(); //pUsersData[user_idx].ratings;
	double rRuNum = 1/sqrt(userRatings+10);
        int dim = ac.D;

        //for( int k=0; k<dim; ++k) {
            double sumY = 0.0;
           //for(unsigned int i=0; i < userRatings; ++i) {
           foreach(graphlab::edge_id_t oedgeid, outs){
               //unsigned int j = pItemRatings_training[currentRatingIdx+i].item;
               vertex_data_svdpp &movie = scope.neighbor_vertex_data(scope.target(oedgeid));
               sumY += sum(movie.weight); //y
           }
            graph_type_svdpp * validation = (graph_type_svdpp*)ps.g<graph_type_svdpp>(VALIDATION);
            if (validation != NULL){
            foreach(graphlab::edge_id_t oedgeid, validation->out_edge_ids(id)){
                //unsigned int item_idx = pItemRatings_validation[currentRatingIdx_val++].item;
		vertex_data_svdpp & movie = validation->vertex_data(validation->target(oedgeid));
		sumY += sum(movie.weight); //y
            }
            }
	    
            graph_type_svdpp * test = (graph_type_svdpp*)ps.g<graph_type_svdpp>(TEST);
            if (test != NULL){
            foreach(graphlab::edge_id_t oedgeid, test->out_edge_ids(id)){
                //unsigned int item_idx = pItemRatings_validation[currentRatingIdx_val++].item;
		vertex_data_svdpp & movie = test->vertex_data(test->target(oedgeid));
                sumY += sum(movie.weight); //y
            }
            }

         /*     TODO  for(unsigned int i=0; i < 4; ++i) {
                    unsigned int j = pItemRatings_validation[currentRatingIdx_val+i].item;
                    sumY += y[j][k];
                }
                for(unsigned int i=0; i < 6; ++i) {
                    unsigned int j = pItemRatings_test[currentRatingIdx_test+i].item;
                    sumY += y[j][k];
                }*/
            for( int k=0; k<dim; ++k) {
                usr.ptemp[k] = usr.pu[k] + rRuNum * sumY; // pTemp = pu + rRuNum*sumY
            }
            vector<double> sum(dim, 0);
            //currentRatingIdx_y = currentRatingIdx;
            //for (unsigned int i=0; i<userRatings; i++) {
           
            foreach(graphlab::edge_id_t oedgeid, outs) {
                //int rui = pItemRatings_training[currentRatingIdx].score;
                edge_data & edge = scope.edge_data(oedgeid);
                float rui = edge.weight;
                //unsigned int item_idx = pItemRatings_training[currentRatingIdx].item;
                vertex_data_svdpp& movie = scope.neighbor_vertex_data(scope.target(oedgeid));
                time_svdpp_movie mov(movie);
                //double pui = predict(pItemRatings_training[currentRatingIdx],user_idx);
                double pui; 
                predict(mov, usr, rui, edge.time, pui);

                //int t = pItemRatings_training[currentRatingIdx].day;
                int t = edge.time;

                double eui = rui - pui;

                //bu[user_idx] += lrate * (eui - timesvdpp_beta* bu[user_idx]);
	        *usr.bu += lrate*(eui - timesvdpp_beta * *usr.bu);

                //bi[item_idx] += lrate * (eui - timesvdpp_beta* bi[item_idx]);
                *mov.bi += lrate * (eui - timesvdpp_beta* *mov.bi);

                //bt[t] += lrate * (eui - timesvdpp_beta * bt[t]);
		//TODO: check this

               time_svdpp_time time(ps.times_svdpp[t]);

                for (int k = 0; k < dim; k++) {
                    //double oldValue = mov.q[item_idx][k];
                    double oldValue = mov.q[k];
                    //double userValue = ptemp[user_idx][k] + pu[user_idx][k] * pt[t][k];
                    double userValue = usr.ptemp[k] + usr.pu[k] * time.pt[k];
                    //sum[k] += eui * q[item_idx][k];
                    sum[k] += eui * mov.q[k];
                    //q[item_idx][k] += lrate * (eui * userValue - garma*q[item_idx][k]);
                    mov.q[k] += lrate * (eui * userValue - garma*mov.q[k]);
                    //ptemp[user_idx][k] += lrate * (eui * oldValue - garma * ptemp[user_idx][k]);
                    usr.ptemp[k] += lrate * ( eui * oldValue - garma * usr.ptemp[k]);
                    //p[user_idx][k] += lrate * (eui * oldValue - garma*p[user_idx][k]);
                    usr.p[k] += lrate * ( eui * oldValue - garma*usr.p[k] );
                    //pu[user_idx][k] += lrate * (eui * oldValue  * pt[t][k] - garma * pu[user_idx][k]);
                    usr.pu[k] += lrate * (eui * oldValue  * time.pt[k] - garma * usr.pu[k]);
                    //pt[t][k] += lrate * (eui * oldValue * pu[user_idx][k] - garma * pt[t][k]);
                    time.pt[k] += lrate * (eui * oldValue * usr.pu[k] - garma * time.pt[k]);
                    //double xOldValue = x[user_idx][k];
                    double xOldValue = usr.x[k];
                    //double zOldValue = z[t][k];
                    double zOldValue = time.z[k];

                    //x[user_idx][k] += lrate * (eui * zOldValue - garma * xOldValue);
                    usr.x[k] += lrate * (eui * zOldValue - garma * xOldValue);
                    //z[t][k] += lrate * (eui * xOldValue - garma * zOldValue);
                    time.z[k] += lrate * (eui * xOldValue - garma * zOldValue);
                }

                //prmse += eui * eui;
                user.rmse += eui*eui;
                //currentRatingIdx += 1;
            }

            //for(unsigned int j=0; j < userRatings; ++j) {
            foreach(graphlab::edge_id_t oedgeid, outs){
                //unsigned int item_idx = pItemRatings_training[currentRatingIdx_y++].item;
                vertex_data_svdpp & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
                time_svdpp_movie mov(movie);
                for(int k=0;k<dim;k++){
                    //y[item_idx][k] += lrate * (rRuNum * sum[k]- garma*y[item_idx][k]);
                    mov.y[k] += lrate * (rRuNum * sum[k]- garma*mov.y[k]);
                }
            }
            //for(unsigned int j=0; j < 4; ++j) {

	    //graph_type_svdpp * validation = ps.g<graph_type_svdpp>(VALIDATION);
            if (validation != NULL){
            foreach(graphlab::edge_id_t oedgeid, validation->out_edge_ids(id)){
                //unsigned int item_idx = pItemRatings_validation[currentRatingIdx_val++].item;
		vertex_data_svdpp & movie = validation->vertex_data(validation->target(oedgeid));
                time_svdpp_movie mov(movie);
                for(int k=0;k<dim;k++){
                    mov.y[k] += lrate * (rRuNum * sum[k]- garma*mov.y[k]);
                }
            }
            }
	    
            //graph_type_svdpp * test = ps.g<graph_type_svdpp>(TEST);
            /*for(unsigned int j=0; j < 6; ++j) {
                unsigned int item_idx = pItemRatings_test[currentRatingIdx_test++].item;
                for(int k=0;k<dim;k++){
                    y[item_idx][k] += lrate * (rRuNum * sum[k]- garma*y[item_idx][k]);
                }
            }*/
            if (test != NULL){
            foreach(graphlab::edge_id_t oedgeid, test->out_edge_ids(id)){
                //unsigned int item_idx = pItemRatings_validation[currentRatingIdx_val++].item;
		vertex_data_svdpp & movie = test->vertex_data(test->target(oedgeid));
                time_svdpp_movie mov(movie);
                for(int k=0;k<dim;k++){
                    mov.y[k] += lrate * (rRuNum * sum[k]- garma*mov.y[k]);
                }
            }
            }

   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (scope.vertex() == (uint)(ps.M-1))
  	time_svd_post_iter();
 
}



//inline double SVDPlusPlusManager::predict(ItemRating itemRating, unsigned int user)
double predict(const time_svdpp_movie & mov, const time_svdpp_usr & usr, float rating, int time, double & estScore){
	//unsigned int item = itemRating.item;
	unsigned int t = (unsigned int)time;//itemRating.day;
	    	//  pui = mu + bu[u] + bi[i];
	double pui  = ps.globalMean[0] + *usr.bu + *mov.bi;// + bt[t];
	time_svdpp_time ptime(ps.times_svdpp[t]);	
        int dim = ac.D;
	for(int k=0;k<dim;k++){
//pui += ptemp[user][k] * q[item][k];
		pui += (usr.ptemp[k] * mov.q[k]);
		pui += usr.x[k] * ptime.z[k];
		pui += usr.pu[k] * ptime.pt[k] * mov.q[k];
	}
	pui = min(pui,ac.maxval);
	pui = max(pui,ac.minval);
	estScore = pui;
	assert(!std::isnan(pui));
        return (rating-estScore)*(rating-estScore);
}
	    					                	        			





/*void predict_missing_value(const vertex_data_svdpp&data, const vertex_data_svdpp& pdata, edge_data& edge, double & sq_err, int&e, int i){
    float prediction = 0;
    predict(data, pdata, &edge, NULL, edge.weight, prediction);
    e++;
}*/
 

//calculate RMSE. This function is called only before and after grahplab is run.
//during run, agg_rmse_by_movie is called 0 which is much lighter function (only aggregate sums of squares)
double calc_time_svd_rmse(const graph_type_svdpp * _g, bool test, double & res){

     graph_type_svdpp * g = (graph_type_svdpp*)ps.g<graph_type_svdpp>(TRAINING);

     if (test && ps.Le == 0)
       return NAN;
      
     
     res = 0;
     double sqErr =0;
     int nCases = 0;

     for (int i=0; i< ps.M; i++){
       foreach(edge_id_t oedgeid, _g->out_edge_ids(i)){
         const edge_data & item = _g->edge_data(oedgeid);
         vertex_data_svdpp & movie = g->vertex_data(_g->target(oedgeid)); 
         vertex_data_svdpp & usr = g->vertex_data(i);
         double estScore;
         sqErr += predict(time_svdpp_movie(movie), time_svdpp_usr(usr), item.weight, item.time, estScore);
         nCases++;
       }
   }
   res = sqErr;
   assert(nCases == (test?ps.Le:ps.L));
   return sqrt(sqErr/(double)nCases);
}


void time_svd_post_iter(){
  printf("Entering last iter with %d\n", ps.iiter);

  double res,res2;
  double rmse = agg_rmse_by_user<graph_type_svdpp, vertex_data_svdpp>(res);
  printf("%g) Iter %s %d, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", ps.gt.current_time(), "SVD", ps.iiter,  rmse, calc_time_svd_rmse(ps.g<graph_type_svdpp>(VALIDATION), true, res2));

  lrate *= 0.9;
  ps.iiter++;
}


void time_svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
   assert(false); //mode not supported
}
void time_svd_plus_plus_update_function(gl_types_mult_edge::iscope &scope, 
			 gl_types_mult_edge::icallback &scheduler) {
   assert(false); //mode not supported
} 
void time_svd_plus_plus_update_function(gl_types_mcmc::iscope &scope, 
			 gl_types_mcmc::icallback &scheduler) {
   assert(false); //mode not supported
} 


#include "graphlab/macros_undef.hpp"
#endif //__SVDPP_HPP
