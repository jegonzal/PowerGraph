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
 */


#ifndef _IO_HPP
#define _IO_HPP

#include "stats.hpp"
#include "../../libs/matrixmarket/mmio.h" //matrix market format support
#include "implicit.hpp"
#include "svdpp.hpp"
#include <graphlab/macros_def.hpp>

extern advanced_config config;
extern problem_setup ps;


const static  int matlab_offset_user_movie = 1; //matlab array start from 1
static  int matlab_offset_time = 1; //matlab arrays start from 1
bool * flags = NULL;

FILE * open_file(const char * name, const char * mode){
  FILE * f = fopen(name, mode);
  if (f == NULL){
      perror("fopen failed");
      logstream(LOG_ERROR) <<" Failed to open file" << name << std::endl;
      exit(1);
   }
  return f;
}


template<typename graph_type, typename vertex_data>
void add_time_nodes(graph_type* _g){
    //init times
    
    int size = ac.D;
    if (ps.algorithm == TIME_SVD_PLUS_PLUS)
      size = 2*ac.D;
    ps.times = new vertex_data[ps.K];
    vec tones = ones(size)*(ps.K==1?1:0.1);
    //add T time node (ps.tensor dim 3)
    for (int i=0; i<ps.K; i++){
      ps.times[i].pvec =tones;
      _g->add_vertex(ps.times[i]);
      //if (ac.debug && (i <= 5 || i == ps.K-1))
      //  debug_print_vec("T: ", ps.times[i].pvec, ac.D);
    }
}; //nothing to be done here 

/*template<>
void add_time_nodes<graph_type_svdpp,vertex_data_svdpp>(graph_type_svdpp* _g){
    if (ps.K <= 1)
      logstream(LOG_FATAL)<<"When running a time based algorithm, don't forget to specify the total number of time bins using --K=XX command, where XX is the number of time bins. Note that time binds start from zero and should be integer." << std::endl;
    ps.times_svdpp = new vertex_data_svdpp[ps.K];
    //add T time node (ps.tensor dim 3)
    for (int i=0; i<ps.K; i++){
      _g->add_vertex(ps.times_svdpp[i]);
    }
}*/
/**
 * Add the graph nodes. We have nodes for each row (user), column (movies) and time bins.
 * 
 */
template<typename graph_type, typename vertex_data>
void add_vertices(graph_type * _g, testtype data_type){
  vertex_data vdata;
  // add M user nodes (ps.tensor dim 1)
  for (int i=0; i<ps.M; i++){
    //vdata.pvec = ac.debug? (ones(ac.D)*0.1) : (randu(ac.D)*(0.1/sqrt(ac.D)));
    _g->add_vertex(vdata);
    //if (ac.debug && (i<= 5 || i == ps.M-1))
    //  debug_print_vec("U: ", vdata.pvec, ac.D);
  }
  
  // add N movie node (ps.tensor dim 2) 
  for (int i=0; i<ps.N; i++){
    //vdata.pvec = ac.debug? (ones(ac.D)*0.1) : (randu(ac.D)*(0.1/sqrt(ac.D)));
    _g->add_vertex(vdata);
    //if (ac.debug && (i<=5 || i==ps.N-1))
    //  debug_print_vec("V: ", vdata.pvec, ac.D);
  }
  
  //add time nodes (if needed)
  if (data_type==TRAINING && ps.tensor){
    add_time_nodes<graph_type, vertex_data>(_g);
  }
}

template<typename graph_type, typename edge_data>
void verify_edges(graph_type * _g, testtype data_type){

  //verify edges
  for (int i=ps.M; i < ps.M+ps.N; i++){
    foreach(graphlab::edge_id_t eid, _g->in_edge_ids(i)){          
      int from = _g->source(eid);
      int to = _g->target(eid);
      assert(from < ps.M);
      assert(to >= ps.M && to < ps.M+ps.N);

      const edge_data & data = _g->edge_data(eid);
	    if (!ac.zero)
          assert(data.weight != 0);  
      if (ps.tensor && ps.algorithm != WEIGHTED_ALS)
          assert(data.time < ps.K);
  
      if (ps.K > 1 && data_type==TRAINING && ps.tensor)
          edges[(int)data.time].push_back(eid);
    }
  }
}
template<>
void verify_edges<graph_type_mult_edge,multiple_edges>(graph_type_mult_edge * _g, testtype data_type){

  //verify edges
  for (int i=ps.M; i < ps.M+ps.N; i++){
    foreach(graphlab::edge_id_t eid, _g->in_edge_ids(i)){          
     const  multiple_edges & tedges= _g->edge_data(eid);
      int from = _g->source(eid);
      int to = _g->target(eid);
      assert(from < ps.M);
      assert(to >= ps.M && to < ps.M+ps.N);

      for (int j=0; j< (int)tedges.medges.size(); j++){
        const edge_data_mcmc & data= tedges.medges[j];
	      if (!ac.zero)
          assert(data.weight != 0);  
        if (ps.algorithm != WEIGHTED_ALS)
          assert(data.time < ps.K);
  
        if (ps.K > 1 && data_type==TRAINING && ps.tensor)
          edges[(int)data.time].push_back(eid);
        }
    }
  }
}

void fill_factors_svdpp();
void fill_factors_time_svd_plus_plus();
void fill_factors_libfm();
#include "read_matrix_market.hpp"
/**
 * fill data structures used for writing output to file
 */
void fill_factors_uvt(const graph_type *g ){
  if (ac.bptf_additional_output != true)
       ((graph_type*)ps.g<graph_type>(TRAINING))->reduce_mem_consumption();
 
  if (ps.algorithm == SVD_PLUS_PLUS || ps.algorithm == BIAS_SGD){
     fill_factors_svdpp();
  }
  else if (ps.algorithm == TIME_SVD_PLUS_PLUS){
     fill_factors_time_svd_plus_plus();
   }
   else if (ps.algorithm == RBM){ //TODO

   }
   else if (ps.algorithm == LIBFM){
     fill_factors_libfm();
   }
   else if (ps.isals || ps.algorithm == STOCHASTIC_GRADIENT_DESCENT || ps.algorithm == NMF){
    ps.U = zeros(ps.M,ac.D);
     ps.V = zeros(ps.N,ac.D);

   for (int i=0; i< ps.M+ps.N; i++){ 
      const vertex_data & data = ps.g<graph_type>(TRAINING)->vertex_data(i);
      if (i < ps.M){
        set_row(ps.U, i, data.pvec);
      }
      else {
        set_row(ps.V, i-ps.M, data.pvec);
      }
   }

   if (ps.tensor){ 
     ps.T = zeros(ps.K,ac.D);
     for (int i=0; i<ps.K; i++){
        set_row(ps.T, i, ps.times[i].pvec);
     }
    } 
  }
  else assert(false);
} 


void fill_factors_uvt(const graph_type_mcmc * g){
   if (ps.algorithm != LANCZOS && ps.algorithm != SVD){
   //clear the graph only at the end of the run
   if (ac.bptf_additional_output != true)
     ((graph_type*)ps.g<graph_type_mcmc>(TRAINING))->reduce_mem_consumption();
   ps.U = zeros(ps.M,ac.D);
   ps.V = zeros(ps.N,ac.D);

   for (int i=0; i< ps.M+ps.N; i++){ 
      const vertex_data & data = ps.g<graph_type_mcmc>(TRAINING)->vertex_data(i);
      if (i < ps.M){
        set_row(ps.U, i, data.pvec);
      }
      else {
        set_row(ps.V, i-ps.M, data.pvec);
      }
   }

   if (ps.tensor){ 
     ps.T = zeros(ps.K,ac.D);
     for (int i=0; i<ps.K; i++){
        set_row(ps.T, i, ps.times[i].pvec);
     }
    } 
  }
}
void fill_factors_uvt(const graph_type_mult_edge * g){
   if (ps.algorithm != LANCZOS && ps.algorithm != SVD){
   //clear the graph only at the end of the run
   if (ac.bptf_additional_output != true)
     ((graph_type*)ps.g<graph_type_mult_edge>(TRAINING))->reduce_mem_consumption();
   ps.U = zeros(ps.M,ac.D);
   ps.V = zeros(ps.N,ac.D);

   for (int i=0; i< ps.M+ps.N; i++){ 
      const vertex_data & data = ps.g<graph_type_mult_edge>(TRAINING)->vertex_data(i);
      if (i < ps.M){
        set_row(ps.U, i, data.pvec);
      }
      else {
        set_row(ps.V, i-ps.M, data.pvec);
      }
   }

   if (ps.tensor){ 
     ps.T = zeros(ps.K,ac.D);
     for (int i=0; i<ps.K; i++){
        set_row(ps.T, i, ps.times[i].pvec);
     }
    } 
  }
}


template<typename edgedata, typename graph_type, typename edge_data>
int read_mult_edges(FILE * f, int nodes, testtype type, graph_type * g, graph_type * _g, bool symmetry = false);
 //write an output vector to file
template <typename T>
void write_vec(FILE * f, const int len, const T * array){
  assert(f != NULL && array != NULL);
  int total = 0;
  
  while(true){
    int rc = fwrite(array+total, sizeof(T), len-total, f);
    if (rc <= 0){
      if (errno == EINTR){
         logstream(LOG_WARNING) << "Interrupted system call, trying agin " << std::endl;
         continue;
      }
      perror("write failed");
      logstream(LOG_FATAL) << "Failed to write vector!" << std::endl;
    }
    total += rc;
    if (total >= len)
      break;
  }
}

void truncate_and_scale(float & prediction){
  if (prediction<ac.minval)
     prediction=ac.minval;
  else if (prediction>ac.maxval)
     prediction=ac.maxval; 

  if (ac.scalerating != 1)
     prediction *= ac.scalerating;
}
float rbm_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);
float bias_sgd_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);
float svdpp_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);
float time_svdpp_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);
float libfm_predict(const vertex_data& user, 
                const vertex_data& movie, 
                const edge_data * edge,
                const vertex_data* nothing,
                const float rating, 
                float & prediction);
    	
template<typename graph_type, typename vertex_data, typename edge_data>
void common_prediction(const graph_type &g, const graph_type & _g, const vertex_data& data,int i, int &lineNum, double& sumPreds, vec* test_predictions, bool dosave, double & RMSE, double &  MAE){
  
  foreach(edge_id_t iedgeid, _g.out_edge_ids(i)) {
    const vertex_data & pdata = g.vertex_data(_g.target(iedgeid)); 
	  const edge_data & edge = _g.edge_data(iedgeid);
          
    if (!ac.zero)
			assert(edge.weight != 0);

		float prediction = 0;
    if (ps.algorithm == BIAS_SGD)
      bias_sgd_predict(data, pdata, (edge_data*)NULL, NULL, edge.weight, prediction);
    else if (ps.algorithm == SVD_PLUS_PLUS)
      svdpp_predict(data, pdata, (edge_data*)NULL, NULL, edge.weight, prediction);
    else if (ps.algorithm == TIME_SVD_PLUS_PLUS)
      time_svdpp_predict(data, pdata, &edge, NULL, edge.weight, prediction);
    else if (ps.algorithm == RBM)
      rbm_predict(data, pdata, (edge_data*)NULL, NULL, edge.weight, prediction);
    else if (ps.algorithm == LIBFM)
      libfm_predict(data, pdata, &edge, &ps.times[(int)edge.time], edge.weight, prediction);
    else
     predict(data, pdata, ps.algorithm == WEIGHTED_ALS ? &edge : NULL, ps.tensor? (&ps.times[(int)edge.time]):NULL, edge.weight, prediction);
    truncate_and_scale(prediction);      
    RMSE += pow(prediction - edge.weight, 2);
    MAE += fabs(prediction - edge.weight);
    if (ac.debug && (i== 0 || i == ps.M))
			cout<<lineNum<<") prediction:"<<prediction<<endl; 
    if (dosave)
			test_predictions->operator[](lineNum) = prediction;
	  
    sumPreds += prediction;
 	  lineNum++; 
  }
}

//compute predictions 
void test_predict(vertex_data & usr, int i, int& lineNum, double & sumPreds, vec* test_predictions, bool dosave, const graph_type &g, const graph_type & _g, double & RMSE, double & MAE){
  if (ps.algorithm == SVD_PLUS_PLUS){
    vertex_data_svdpp user = usr;
		int n = user.num_edges; //+1.0 ? //regularization
		memset(&user.weight[0], 0, ac.D*sizeof(double));
		if (n > 0 ){
			foreach(edge_id_t oedgeid, g.out_edge_ids(i)) {
				vertex_data_svdpp & movie = (vertex_data_svdpp&)g.vertex_data(g.target(oedgeid)); 
				for (int j=0; j< ac.D; j++)
          user.weight[j] += movie.weight[j];
			}
			float usrnorm = float(1.0/sqrt(n));
			for (int j=0; j< ac.D; j++)
			  user.weight[j] *= usrnorm;
		}
  }
  else if (ps.algorithm == BIAS_SGD){}
   else if (ps.algorithm == TIME_SVD_PLUS_PLUS){ }
  else if (ps.algorithm == RBM){ }
  common_prediction<graph_type,vertex_data,edge_data>(g, _g,usr,i,lineNum, sumPreds, test_predictions, dosave, RMSE, MAE);


}

//compute predictions for BPTF/PMF
void test_predict(vertex_data & data, int i, int&lineNum, double& sumPreds, vec * test_predictions, bool dosave, const graph_type_mcmc& g, const graph_type_mcmc &_g, double & RMSE, double & MAE){

      foreach(edge_id_t iedgeid, _g.out_edge_ids(i)) {
        const vertex_data & pdata = g.vertex_data(_g.target(iedgeid)); 
	  edge_data_mcmc & edge = (edge_data_mcmc&)_g.edge_data(iedgeid);
          
          if (!ac.zero)
           	assert(edge.weight != 0);

          float prediction = 0;
          predict(data, 
                  pdata, 
                  NULL, 
                  ps.tensor? (&ps.times[(int)edge.time]):NULL, 
                  edge.weight, 
                  prediction);
          
          if (ps.BPTF && ps.iiter > ac.bptf_burn_in){
             edge.avgprd += prediction;
             //add = powf((edge.avgprd / (iiter - bptf_burn_in)) - edge.weight, 2);
              prediction = (edge.avgprd / (ps.iiter - ac.bptf_burn_in));
           }
          
          RMSE+= pow(prediction - edge.weight, 2);
          MAE += fabs(prediction - edge.weight);
	        truncate_and_scale(prediction);
          if (ac.debug && (i== 0 || i == ps.M))
            cout<<lineNum<<") prediction:"<<prediction<<endl; 
          if (dosave)
            test_predictions[lineNum] = prediction;
	        sumPreds += prediction;
 	        lineNum++; 
       }
}

//compute predictions of tensors with multiple edges between same pair of nodes
void test_predict(vertex_data & data, int i, int&lineNum, double & sumPreds, vec* test_predictions, bool dosave, const graph_type_mult_edge&g, const graph_type_mult_edge &_g, double & RMSE, double & MAE){
      foreach(edge_id_t iedgeid, _g.out_edge_ids(i)) {
        const multiple_edges & edges = _g.edge_data(iedgeid);
        const vertex_data & pdata = g.vertex_data(_g.target(iedgeid)); 
        for (int j=0; j< (int)edges.medges.size(); j++){  
          edge_data_mcmc & edge = (edge_data_mcmc&)edges.medges[j];
          
          if (!ac.zero)
           	assert(edge.weight != 0);

          float prediction = 0;
          predict(data, 
                  pdata, 
                  NULL, 
                  ps.tensor? (&ps.times[(int)edge.time]):NULL, 
                  edge.weight, 
                  prediction);
          if (ps.BPTF && ps.iiter > ac.bptf_burn_in){
             edge.avgprd += prediction;
              prediction = (edge.avgprd / (ps.iiter - ac.bptf_burn_in));
           }
	        truncate_and_scale(prediction);
          if (ac.debug && (i== 0 || i == ps.M))
            cout<<lineNum<<") prediction:"<<prediction<<endl; 
          if (dosave)
           test_predictions[lineNum] = prediction;
          RMSE += pow(prediction - edge.weight, 2);
          MAE += fabs(prediction - edge.weight);
	        sumPreds += prediction;
 	        lineNum++; 
         }
       }
 }


template<typename graph_type, typename vertex_data, typename edge_data>
void export_test_file(const graph_type & _g, testtype type, bool dosave) {
       
  const graph_type * g = ps.g<graph_type>(TRAINING);
  double sumPreds = 0;
  int lineNum = 0;
  if (!dosave)
    assert(ps.BPTF);	
  if (ps.Lt == 0 && type == TEST)
    logstream(LOG_FATAL) << "Empty or missing test data file, can not compute predictions!" << std::endl;
  if (ps.Lt2 == 0 && type == TEST2)
    logstream(LOG_FATAL) << "Empty or missing test data file, can not compute predictions!" << std::endl;

  int size = 0;
  string suffix = "";
  string comment = "";
  switch(type){
  
    case VALIDATION: 
       size = ps.Le; 
       suffix = ".validation.predictions"; 
       comment = "output predictions for validation data\n";
       break;
  
    case TEST: 
       size = ps.Lt; 
       suffix = ".test.predictions";
       comment = "output predictions for test data\n";
       break;

    case TEST2: 
       size = ps.Lt2; 
       suffix = ".test2.predictions"; 
       comment = "output predictions for test2 data\n";
       break;

    case TRAINING:
    default:
      assert(false);
  }
  vec out_predictions = zeros(size);
  double RMSE = 0, MAE = 0;
  for (int i=0; i< ps.M; i++){ 
      vertex_data & data = (vertex_data&)g->vertex_data(i);
      test_predict(data, i, lineNum, sumPreds, &out_predictions, dosave, *g, _g, RMSE, MAE);
  }

  ASSERT_EQ(lineNum, size);
  logstream(LOG_INFO)<< "**Completed successfully (mean prediction: " << (sumPreds/lineNum)-ac.shiftrating << std::endl;
  if (dosave)
    save_matrix_market_vector((ac.datafile + suffix).c_str(), 
     out_predictions, comment, false, false);
}



template<typename graph_type, typename vertex_data>
void export_uvt_to_matrixmarket(const graph_type * g){
  if (ps.algorithm != LANCZOS && ps.algorithm != SVD)
     fill_factors_uvt(g);
  char dfile[256] = {0};
  sprintf(dfile,"%s-%d-%d.out",ac.datafile.c_str(), ac.D,ps.iiter);
  remove(dfile);
  save_matrix_market_format(dfile, ps.U, ps.V);  
 
}



void set_num_edges(int val, testtype data_type){
  switch(data_type){
    case TRAINING: 
      ps.L = val; break;
    
    case VALIDATION: 
      ps.Le = val; 
      if (ac.aggregatevalidation)
	ps.L+=  ps.Le; //add edges of validation dataset into the training data set as well.
	break; 
     
    case TEST: 
      ps.Lt = val; break;

    case TEST2:
      ps.Lt2 = val; break;
  }  
}

/**
 * Verify that matrix size is consistent between training, validation and test files
 */
void verify_size(testtype data_type, int _M, int _N, int _K){
  if (data_type != TRAINING && ps.M != _M){
	  logstream(LOG_WARNING) << " wrong number of users: " << _M << " instead of " << ps.M << " in " << testtypename[data_type] << std::endl;
    if (ps.M < _M)
      logstream(LOG_FATAL)<<"Can not continue. Please fix your input file!" << std::endl;
  }
  if (data_type != TRAINING && ps.N != _N){
	  logstream(LOG_WARNING) << " wrong number of movies: " << _N << " instead of " << ps.N << " in " << testtypename[data_type] << std::endl;
    if (ps.N < _N)
      logstream(LOG_FATAL)<<"Can not continue. Please fix your input file!" << std::endl;
   }
  if (data_type != TRAINING && ps.K != _K && ac.K == 0){
	  logstream(LOG_WARNING) << " wrong number of time bins: " << _K << " instead of " << ps.K << " in " << testtypename[data_type] <<std::endl;
     if (ps.K < _K)
      logstream(LOG_FATAL)<<"Can not continue. Please fix your input file!" << std::endl;
  }
  printf("Matrix size is: USERS %d MOVIES %d TIME BINS %d\n", ps.M, ps.N, ps.K);
}




/* function that reads the ps.tensor from file */
/* Input format is:
 * M - number of users (int)
 * N - number of movies (int)
 * K - number of time bins (int), in case of weighted ALS this value is ignored
 * e - number of edges (int)
 * A list of edges in the format
 * [from] [to] [ time] [weight]  (4 floats)
 * where [from] is an integeter from 1 to M
 * [to] is an interger from 1 to N
 * [time] is an integer from 1 to K (for weighted ALS this is the weight, which is float)
 * [weight] - this is the rating, which is float. Rating is assumed non-zero unless the --zero=true flas is on 
 */
template<typename graph_type, typename gl_types, typename vertex_data, typename edge_data>
void load_pmf_graph(const char* filename, graph_type * g, graph_type * _g, testtype data_type) {


  if (ac.matrixmarket){
      //printf("Loading Matrix Market file %s %s\n", filename, testtypename[data_type]);
      load_matrix_market<graph_type, vertex_data, edge_data>(filename, _g, data_type);
      return;
  }

}

template<typename edgedata, typename edge_data>
void verify_edge(edgedata & ed, edge_data& edge, int i, testtype type){
    if (!ac.zero) //usually we do not allow zero ratings, unless --zero=true flag is set.
	 assert(ed.weight != 0); 
      //verify node ids are in allowed range
      if (i == 0 && ((int)ed.from < matlab_offset_user_movie || (int)ed.from > ps.last_node))
          logstream(LOG_FATAL) << " Wrong intput file format. Did you try to use --float=true " << endl;
      if ((int)ed.from < matlab_offset_user_movie || (int)ed.from > ps.last_node)
          logstream(LOG_FATAL) << " Edge from node " << ed.from << " where the allowed node range is [" <<
         matlab_offset_user_movie << "-" << ps.last_node << "]. In input line " << i << endl;
      if ((int)ed.to < matlab_offset_user_movie || (int)ed.to > ps.last_node)
          logstream(LOG_FATAL) << " Edge to node " << ed.to << " where the allowed node range is [" <<
         matlab_offset_user_movie << "-" << ps.last_node << "]. In input line " << i << endl;
       //no self edges
      if ((int)ed.to == (int)ed.from)
        logstream(LOG_FATAL) << " Self edge between node " << ed.to << " to itself is not allowed. In input line " << i << endl;

      edge.weight = (double)ed.weight;
    
      //if sacling of rating values is requested to it here.
      if (ac.scalerating != 1.0)
	     edge.weight /= ac.scalerating;
      ps.globalMean[type] += edge.weight;
     
      //if scaling of time bins request do it here
      double time  = ((ed.time - matlab_offset_time - ac.truncating)/(double)ac.scaling);
      edge.time = time;
}



template<typename edgedata, typename graph_type>
void add_edge(int i, edgedata &ed, graph_type *g, graph_type *_g, multiple_edges edges, testtype type){ 
      
      edge_data_mcmc edge;
     
      verify_edge<edgedata,edge_data_mcmc>(ed, edge, i, type); 
       std::pair<bool, edge_id_t> ret;
      if (flags[(int)ed.from-matlab_offset_user_movie] == true && flags[(int)ed.to-matlab_offset_user_movie] == true){
        ret = _g->find((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie);
      }
      else ret.first = false;

      if (ret.first == false){
        edges.medges.push_back(edge); 
        _g->add_edge((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie, edges); // Matlab export has ids starting from 1, ours start from 0
        if (type == VALIDATION && ac.aggregatevalidation)//add validation edges into training dataset as well
           g->add_edge((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie, edges); // Matlab export has ids starting from 1, ours start from 0
          flags[(int)ed.from-matlab_offset_user_movie] = true;
          flags[(int)ed.to-matlab_offset_user_movie] = true;
      }
      else {
        _g->edge_data(ret.second).medges.push_back(edge);
      }
 }





template<typename edgedata, typename graph_type>
void add_edge(int i, edgedata &ed, graph_type *g, graph_type *_g, edge_data edge, testtype type){ 
  verify_edge<edgedata, edge_data>(ed, edge, i,type);
  _g->add_edge((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie, edge);
  if (type == VALIDATION && ac.aggregatevalidation)//add validation edges into training dataset as well
  g->add_edge((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie, edge);
}

template<typename edgedata, typename graph_type>
void add_edge(int i, edgedata &ed, graph_type *g, graph_type *_g, edge_data_mcmc edge, testtype type){ 
  verify_edge<edgedata, edge_data_mcmc>(ed, edge, i,type);
  _g->add_edge((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie, edge);
  if (type == VALIDATION && ac.aggregatevalidation)//add validation edges into training dataset as well
  g->add_edge((int)ed.from-matlab_offset_user_movie, (int)ed.to-matlab_offset_user_movie, edge);
}


//LOAD FACTORS FROM FILE
template<typename graph_type>
void import_uvt_from_file(){
 
 const graph_type * g =  ps.g<graph_type>(TRAINING);
 mat U,V,T;
   load_matrix_market_matrix(ac.datafile + ".U", U);
   load_matrix_market_matrix(ac.datafile + ".V", V);
   if (ps.tensor)
     load_matrix_market_matrix(ac.datafile + ".T", T); 
   ASSERT_EQ(U.rows(), ps.M);
   ASSERT_EQ(V.rows(), ps.N);
   ASSERT_EQ(U.cols(), V.cols());
   ASSERT_EQ(U.cols(), ac.D);
 //initalizing feature vectors from file
#pragma omp parallel for
 for (int i=0; i< ps.M+ps.N; i++){ 
    vertex_data & data = (vertex_data&)g->vertex_data(i);
    if (i < ps.M)
        data.pvec = get_row(U, i); 
   else
        data.pvec = get_row(V, i-ps.M);
 }

 if (ps.tensor){ 
    for (int i=0; i<ps.K; i++){
        ps.times[i].pvec = get_row(T, i);
    }
 } 
}


template<typename graph_type, typename vertex_data>
void write_output(const graph_type * g){
  //write output matrices U,V,T to file
     export_uvt_to_matrixmarket<graph_type, vertex_data>(g);
}
#include <graphlab/macros_undef.hpp>
#endif
