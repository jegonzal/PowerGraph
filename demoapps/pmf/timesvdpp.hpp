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

inline double sum(double * pvec){
  double tsum = 0;
  for (int j=0; j< ac.D; j++)
    tsum += pvec[j];
  return tsum;
}


struct time_svdpp_usr{
   float * bu;
   double * p;
   double * pu;
   double * x;
   double * ptemp;
   int * tu;
   float * rmse;

   time_svdpp_usr(vertex_data & vdata){
      bu = &vdata.bias;
      assert(vdata.pvec.size() == ac.D*4); //TO REMOVE
      p = &vdata.pvec[0];
      pu = p+ac.D;
      x = pu+ac.D;
      ptemp = x+ac.D;
      tu = &vdata.num_edges;
      rmse = &vdata.rmse;
   }

};

struct time_svdpp_movie{
   float * bi;
   double * q;
   double * y;

   time_svdpp_movie(vertex_data& vdata){
     bi = &vdata.bias;
     q = &vdata.pvec[0];
     y = q+ac.D;
   }
};

struct time_svdpp_time{
   float * bt;
   double * z;
   double * pt;

   time_svdpp_time(vertex_data& vdata){
     bt = &vdata.bias;
     z = &vdata.pvec[0];
     pt = z+ac.D;
   }
};

template<typename graph_type>
void init_time_svdpp(graph_type * _g){ assert(false); }

template<>
void init_time_svdpp<graph_type>(graph_type *_g){
   fprintf(stderr, "time-SVD++ %d factors\n", ac.D);

   int k = ac.D;

	for (int u = 0; u < ps.M; u++) {
    vertex_data data = _g->vertex_data(u);
    data.pvec = zeros(4*k);
    time_svdpp_usr usr(data);
    *usr.bu = 0;
    for (int m=0; m< k; m++){
          usr.p[m] = 0.01*graphlab::random::rand01() / (double) (k);
	        usr.pu[m] = 0.001 * graphlab::random::rand01() / (double) (k);
	        usr.x[m] = 0.001 * graphlab::random::rand01() / (double) (k);
	        usr.ptemp[m] = usr.p[m];
    }
  }

	for (int i = ps.M; i < ps.N+ps.M; i++) {
    vertex_data & data = _g->vertex_data(i);
    data.pvec = zeros(2*k);
    time_svdpp_movie movie(data);
		*movie.bi = 0;
		for (int m = 0; m < k; m++){
			movie.q[m] = 0.01 * graphlab::random::rand01() / (double) (k);
			movie.y[m] = 0.001 * graphlab::random::rand01() / (double) (k);
		}
	}
  for (int i = ps.M+ps.N; i < ps.M+ps.N+ps.K; i++) {
    vertex_data & data = _g->vertex_data(i);
    data.pvec = zeros(2*k);
		time_svdpp_time timenode(data);
		*timenode.bt = 0;
    for (int m = 0; m < k; m++){
			timenode.z[m] = 0.001 * graphlab::random::rand01() / (double) (k);
			timenode.pt[m] = 0.001 * graphlab::random::rand01() / (double) (k);
    }
	}

}


/*
 * predict missing rating for time-SVD++ algorithm
 */
float predict(const time_svdpp_usr& user, 
                const time_svdpp_movie& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);


float time_svdpp_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){
   return predict(time_svdpp_usr((vertex_data&)user), time_svdpp_movie((vertex_data&)movie), edge,nothing, rating, prediction);
}
float time_svdpp_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data_mcmc * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){ assert(false); }
 

                   
/***
 * UPDATE FUNCTION
 */
void time_svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    
  /* GET current vertex data */
  time_svdpp_usr usr(scope.vertex_data()); 
  int id = scope.vertex(); 
  /* print statistics */
  if (ac.debug&& (id == 0 || id == ps.M-1)){
    printf("SVD++: entering %s node  %u \n", "user", id);   
  }

  *usr.rmse = 0;
  gl_types::edge_list outs = scope.out_edge_ids();
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
		vertex_data_svdpp movie = scope.neighbor_vertex_data(scope.target(oedgeid));
		sumY += sum(movie.weight); //y
  }
            
  graph_type * validation = (graph_type*)ps.g<graph_type>(VALIDATION);
  if (validation != NULL && validation->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, validation->out_edge_ids(id)){
			vertex_data_svdpp movie = validation->vertex_data(validation->target(oedgeid));
			sumY += sum(movie.weight); //y
		}
  }
	    
  graph_type * test = (graph_type*)ps.g<graph_type>(TEST);
  if (test != NULL && test->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, test->out_edge_ids(id)){
			vertex_data_svdpp movie = test->vertex_data(test->target(oedgeid));
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
     time_svdpp_movie mov = scope.neighbor_vertex_data(scope.target(oedgeid));
     float pui = 0; 
     predict(usr, mov, &edge, NULL, rui, pui);

     int t = edge.time;
     double eui = rui - pui;
	   *usr.bu += ac.tsp.lrate*(eui - beta* *usr.bu);
     *mov.bi += ac.tsp.lrate * (eui - beta* *mov.bi);
     time_svdpp_time time(ps.times[t]);

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

     *usr.rmse += eui*eui;
   }

  foreach(graphlab::edge_id_t oedgeid, outs){
		time_svdpp_movie mov = scope.neighbor_vertex_data(scope.target(oedgeid));
		for(int k=0;k<dim;k++){
      mov.y[k] += ac.tsp.lrate * (rRuNum * sum[k]- ac.tsp.garma*mov.y[k]);
    }
  }
            
  if (validation != NULL && validation->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, validation->out_edge_ids(id)){
			time_svdpp_movie mov = validation->vertex_data(validation->target(oedgeid));
			for(int k=0;k<dim;k++){
				mov.y[k] += ac.tsp.lrate * (rRuNum * sum[k]- ac.tsp.garma*mov.y[k]);
			}
		}
  }
	    
  if (test != NULL && test->num_vertices() > 0){
		foreach(graphlab::edge_id_t oedgeid, test->out_edge_ids(id)){
			time_svdpp_movie mov = test->vertex_data(test->target(oedgeid));
      for(int k=0;k<dim;k++){
			  mov.y[k] += ac.tsp.lrate * (rRuNum * sum[k]- ac.tsp.garma*mov.y[k]);
      }
    }
  }

   ps.counter[EDGE_TRAVERSAL] += t.current_time();

   if (id == ps.M-1)
  	time_svd_post_iter();
 
}


  float predict(const time_svdpp_usr & usr, 
                const time_svdpp_movie & mov, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction){


	unsigned int t = edge->time;//itemRating.day;
	double pui  = ps.globalMean[0] + *usr.bu + *mov.bi;// + bt[t];
	time_svdpp_time ptime(ps.times[t]);	
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
  double training_rmse = agg_rmse_by_user<graph_type, vertex_data>(res);
  double validation_rmse =  calc_rmse<graph_type,vertex_data>(ps.g<graph_type>(VALIDATION),true, res2);
  printf(ac.printhighprecision ? 
        "%g) Iter %s %d  TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n":
        "%g) Iter %s %d  TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n",
  ps.gt.current_time(), runmodesname[ps.algorithm], ps.iiter, training_rmse, validation_rmse);
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
  ac.tsp.lrate *= ac.tsp.lrate_mult_dec;
  ps.iiter++;
}


void time_svd_plus_plus_update_function(gl_types_mult_edge::iscope &scope, 
			 gl_types_mult_edge::icallback &scheduler) {
   assert(false); //mode not supported
} 
void time_svd_plus_plus_update_function(gl_types_mcmc::iscope &scope, 
			 gl_types_mcmc::icallback &scheduler) {
   assert(false); //mode not supported
} 

void fill_factors_time_svd_plus_plus(){
      ps.timesvdpp_out.ptemp = zeros(ps.M,ac.D);
      ps.timesvdpp_out.x = zeros(ps.M,ac.D);
      ps.timesvdpp_out.pu = zeros(ps.M,ac.D);
      ps.timesvdpp_out.q = zeros(ps.N,ac.D);
      ps.timesvdpp_out.z = zeros(ac.K,ac.D);
      ps.timesvdpp_out.pt = zeros(ac.K,ac.D);
      for (int i=0; i< ps.M; i++){
        const time_svdpp_usr  data = (vertex_data&)ps.g<graph_type>(TRAINING)->vertex_data(i);
       
        set_row(ps.timesvdpp_out.ptemp, i, _vec(data.ptemp,ac.D));
        set_row(ps.timesvdpp_out.x, i, _vec(data.x, ac.D));
        set_row(ps.timesvdpp_out.pu, i, _vec(data.pu, ac.D));
      }
      for (int i=ps.M; i < ps.M+ps.N; i++){
        const time_svdpp_movie movie = (vertex_data&)ps.g<graph_type>(TRAINING)->vertex_data(i);
        set_row(ps.timesvdpp_out.q, i-ps.M, _vec(movie.q, ac.D));
      }
      for (int i=0; i< ac.K; i++){
        const time_svdpp_time  data = ps.times[i];
        set_row(ps.timesvdpp_out.z, i, _vec(data.z,ac.D));
        set_row(ps.timesvdpp_out.pt, i, _vec(data.pt, ac.D));
      }
}
#include "graphlab/macros_undef.hpp"
#endif //__SVDPP_HPP
