/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
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

		*timenode.bt = 0;
    for (int m = 0; m < k; m++){
			timenode.z[m] = 0.001 * rand01() / (double) (k);
			timenode.pt[m] = 0.001 * rand01() / (double) (k);
    }
 
	}

}


/*
 * predict missing rating for time-SVD++ algorithm
 */
float time_svdpp_predict(const vertex_data_svdpp& user, 
                const vertex_data_svdpp& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);


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
  if (outs.size() == 0){
    if (id == ps.M-1)
      time_svd_post_iter(); 
      return;
  }

  timer t;
  t.start(); 
 

	unsigned int userRatings = outs.size(); //pUsersData[user_idx].ratings;
	double rRuNum = 1/sqrt(userRatings+10);
  int dim = ac.D;
  double sumY = 0.0;
  foreach(graphlab::edge_id_t oedgeid, outs){
		vertex_data_svdpp &movie = scope.neighbor_vertex_data(scope.target(oedgeid));
		sumY += sum(movie.weight); //y
  }
            
  graph_type_svdpp * validation = (graph_type_svdpp*)ps.g<graph_type_svdpp>(VALIDATION);
  if (validation != NULL && validation->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, validation->out_edge_ids(id)){
			vertex_data_svdpp & movie = validation->vertex_data(validation->target(oedgeid));
			sumY += sum(movie.weight); //y
		}
  }
	    
  graph_type_svdpp * test = (graph_type_svdpp*)ps.g<graph_type_svdpp>(TEST);
  if (test != NULL && test->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, test->out_edge_ids(id)){
			vertex_data_svdpp & movie = test->vertex_data(test->target(oedgeid));
			sumY += sum(movie.weight); //y
		}
  }

   for( int k=0; k<dim; ++k) {
     usr.ptemp[k] = usr.pu[k] + rRuNum * sumY; // pTemp = pu + rRuNum*sumY
   }
   vector<double> sum(dim, 0);
   foreach(graphlab::edge_id_t oedgeid, outs) {
     edge_data & edge = scope.edge_data(oedgeid);
     float rui = edge.weight;
     vertex_data_svdpp& movie = scope.neighbor_vertex_data(scope.target(oedgeid));
     time_svdpp_movie mov(movie);
     float pui = 0; 
     time_svdpp_predict(user, movie, &edge, NULL, rui, pui);

     int t = edge.time;
     double eui = rui - pui;
	   *usr.bu += ac.tsp.lrate*(eui - beta* *usr.bu);
     *mov.bi += ac.tsp.lrate * (eui - beta* *mov.bi);
     time_svdpp_time time(ps.times_svdpp[t]);

     for (int k = 0; k < dim; k++) {
       double oldValue = mov.q[k];
       double userValue = usr.ptemp[k] + usr.pu[k] * time.pt[k];
       sum[k] += eui * mov.q[k];
       mov.q[k] += ac.tsp.lrate * (eui * userValue - ac.tsp.garma*mov.q[k]);
       usr.ptemp[k] += ac.tsp.lrate * ( eui * oldValue - ac.tsp.garma * usr.ptemp[k]);
       usr.p[k] += ac.tsp.lrate * ( eui * oldValue - ac.tsp.garma*usr.p[k] );
       usr.pu[k] += ac.tsp.lrate * (eui * oldValue  * time.pt[k] - ac.tsp.garma * usr.pu[k]);
       time.pt[k] += ac.tsp.lrate * (eui * oldValue * usr.pu[k] - ac.tsp.garma * time.pt[k]);
       double xOldValue = usr.x[k];
       double zOldValue = time.z[k];
       usr.x[k] += ac.tsp.lrate * (eui * zOldValue - ac.tsp.garma * xOldValue);
       time.z[k] += ac.tsp.lrate * (eui * xOldValue - ac.tsp.garma * zOldValue);
     }

     user.rmse += eui*eui;
   }

  foreach(graphlab::edge_id_t oedgeid, outs){
		vertex_data_svdpp & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
		time_svdpp_movie mov(movie);
		for(int k=0;k<dim;k++){
      mov.y[k] += ac.tsp.lrate * (rRuNum * sum[k]- ac.tsp.garma*mov.y[k]);
    }
  }
            
  if (validation != NULL && validation->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, validation->out_edge_ids(id)){
			vertex_data_svdpp & movie = validation->vertex_data(validation->target(oedgeid));
			time_svdpp_movie mov(movie);
			for(int k=0;k<dim;k++){
				mov.y[k] += ac.tsp.lrate * (rRuNum * sum[k]- ac.tsp.garma*mov.y[k]);
			}
		}
  }
	    
  if (test != NULL && test->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, test->out_edge_ids(id)){
			vertex_data_svdpp & movie = test->vertex_data(test->target(oedgeid));
     time_svdpp_movie mov(movie);
     for(int k=0;k<dim;k++){
			 mov.y[k] += ac.tsp.lrate * (rRuNum * sum[k]- ac.tsp.garma*mov.y[k]);
     }
    }
   }

   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (id == ps.M-1)
  	time_svd_post_iter();
 
}


  float time_svdpp_predict(const vertex_data_svdpp& user, 
                const vertex_data_svdpp& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){


  time_svdpp_movie mov((vertex_data_svdpp&)movie);
  time_svdpp_usr usr((vertex_data_svdpp&)user);
	unsigned int t = edge->time;//itemRating.day;
	double pui  = ps.globalMean[0] + *usr.bu + *mov.bi;// + bt[t];
	time_svdpp_time ptime(ps.times_svdpp[t]);	
  int dim = ac.D;
	for(int k=0;k<dim;k++){
		pui += (usr.ptemp[k] * mov.q[k]);
		pui += usr.x[k] * ptime.z[k];
		pui += usr.pu[k] * ptime.pt[k] * mov.q[k];
	}
	pui = min(pui,ac.maxval);
	pui = max(pui,ac.minval);
	prediction = pui;
	assert(!std::isnan(prediction));
  float err = rating - prediction;
  return err*err;
}
	    					                	        			

void time_svd_post_iter(){
  printf("Entering last iter with %d\n", ps.iiter);

  double res,res2;
  double rmse = agg_rmse_by_user<graph_type_svdpp, vertex_data_svdpp>(res);
  printf("%g) Iter %s %d, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", ps.gt.current_time(), "time-SVD++", ps.iiter,  rmse, calc_rmse<graph_type_svdpp,vertex_data_svdpp>(ps.g<graph_type_svdpp>(VALIDATION), true, res2));
  if (ac.calc_ap){
     logstream(LOG_INFO)<<"AP@3 for training: " << calc_ap<graph_type_svdpp,vertex_data_svdpp,edge_data>(ps.g<graph_type_svdpp>(TRAINING)) << " AP@3 for validation: " << calc_ap<graph_type_svdpp,vertex_data_svdpp,edge_data>(ps.g<graph_type_svdpp>(VALIDATION)) << std::endl;
  }
  ac.tsp.lrate *= ac.tsp.lrate_mult_dec;
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
