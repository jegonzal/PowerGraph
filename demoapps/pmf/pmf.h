/* Copyright (c) 2009 Carnegie Mellon University. 
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

#ifndef PMF_H__	 
#define PMF_H__

//#define NDEBUG
#include "graphlab.hpp"
#include "graphlab/core_base.hpp"

/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU

For BPTF algorithm:
See algrithm description and explanation in: 
1) Liang Xiong, Xi Chen, Tzu-Kuo Huang, Jeff Schneider, Jaime G. Carbonell, Temporal Collaborative Filtering with Bayesian Probabilistic Tensor Factorization. In Proceedings of SIAM Data Mining, 2010.
2) Salakhutdinov and Mnih, Bayesian Probabilistic Matrix Factorization using Markov Chain Monte Carlo. in International Conference on Machine Learning, 2008.

For ALS (Alternating least squares) see algorithm and explanation:
3) Yunhong Zhou, Dennis Wilkinson, Robert Schreiber and Rong Pan. Large-Scale Parallel Collaborative Filtering for the Netflix Prize. Proceedings of the 4th international conference on Algorithmic Aspects in Information and Management. Shanghai, China pp. 337-348, 2008. 

For SVD++ see algorithm and explanations:
4) Koren, Yehuda. "Factorization meets the neighborhood: a multifaceted collaborative filtering model." In Proceeding of the 14th ACM SIGKDD 
international conference on Knowledge discovery and data mining, 426434. ACM, 2008. http://portal.acm.org/citation.cfm?id=1401890.1401944

For SGD, see algorithm:
5) Matrix Factorization Techniques for Recommender Systems
Yehuda Koren, Robert Bell, Chris Volinsky
In IEEE Computer, Vol. 42, No. 8. (07 August 2009), pp. 30-37. 
6) Tikk, D. (2009). Scalable Collaborative Filtering Approaches for Large Recommender Systems. Journal of Machine Learning Research, 10, 623-656.

For Lanczos algorithm (SVD) see:
7) http://en.wikipedia.org/wiki/Lanczos_algorithm 

For NMF (non-negative matrix factorization) see:
8) Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.

For Weighted alternating least squares see:
9) Rong Pan, Yunhong Zhou, Bin Cao, Nathan N. Liu, Rajan Lukose, Martin Scholz, and Qiang Yang. 2008. One-Class Collaborative Filtering. In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM '08). IEEE Computer Society, Washington, DC, USA, 502-511. 

For implicit collaborative filering see:
10) One-Class Collaborative Filtering
by: Rong Pan, Yunhong Zhou, Bin Cao, N. N. Liu, R. Lukose, M. Scholz, Qiang Yang. Data Mining, IEEE International Conference on In Data Mining, 2008. ICDM '08.

For sparsity enforcing priors see:
11) Xi Chen, Yanjun Qi, Bing Bai, Qihang Lin and Jaime Carbonell. Sparse Latent Semantic Analysis. In SIAM International Conference on Data Mining (SDM), 2011.

SVD is implemented using two sided Lanczos, see:
12) http://en.wikipedia.org/wiki/Singular_value_decomposition

Koren time-SVD++ is described in the paper:
13) Yehuda Koren. 2009. Collaborative filtering with temporal dynamics. In Proceedings of the 15th ACM SIGKDD international conference on Knowledge discovery and data mining (KDD '09). ACM, New York, NY, USA, 447-456. DOI=10.1145/1557019.1557072 
*/
#include <vector>
#define GL_NO_MULT_EDGES //comment this flag, if you want to have support for multiple edges in different times between the same user and movie

#include "mathlayer.hpp"
#include "../gabp/advanced_config.h"
#include <boost/algorithm/string/predicate.hpp>

//starts for holding edge data in file
struct edge_double{
  int from;
  int to;
  double time;
  double weight;
};
struct edge_float{
  float from;
  float to;
  float time;
  float weight;
};



/** Vertex and edge data types **/
struct vertex_data {
  vec pvec; //vector of learned values U,V,T
  float rmse; //root of mean square error
  int num_edges; //number of edges

  //constructor
  vertex_data();

  void save(graphlab::oarchive& archive) const; 
   
  void load(graphlab::iarchive& archive); 
};


struct vertex_data_svdpp {
  vec pvec; //vector of learned values U,V,T
  float rmse; //root of mean square error
  int num_edges; //number of edges
  float bias; //bias for this user/movie
  vec weight; //weight vector for this user/movie

  //constructor
  vertex_data_svdpp();

  void save(graphlab::oarchive& archive) const; 
   
  void load(graphlab::iarchive& archive); 
};


struct edge_data {
  float  weight;  //observation 
  float  time; //time of observation (for tensor algorithms)
  edge_data(){ 
	weight = 0; 
	time = 0; 
  }
 void save(graphlab::oarchive& archive) const {  
    archive << weight << time; 
  }  
   
  void load(graphlab::iarchive& archive) {  
    archive >> weight >> time;
  }
};

struct edge_data_mcmc {
  float  weight;  //observation 
  float  time; //time of observation (for tensor algorithms)
  float avgprd;

  edge_data_mcmc(){ 
	weight = 0; 
	time = 0; 
	avgprd = 0;
  }

 void save(graphlab::oarchive& archive) const {  
    archive << weight << time; 
    archive<< avgprd;; 
  }  
   
  void load(graphlab::iarchive& archive) {  
    archive >> weight >> time;
    archive >> avgprd; 
  }
};

//containiner for handling multiple edge
struct multiple_edges{
  std::vector<edge_data_mcmc> medges;
 void save(graphlab::oarchive& archive) const {  
       archive << medges; 
  }  
   
  void load(graphlab::iarchive& archive) {  
      archive >> medges;  
  }
};



 
 double get_rmse(const vertex_data & v);


//data file types
enum testtype{
    TRAINING = 0,
    VALIDATION = 1,
    TEST = 2,
    TEST2 = 3//second test file, for kdd cup 2012
};

enum linear_algebra_support{
  ITPP_SUPPORT = 1,
  EIGEN_SUPPORT = 2
};

//run modes
enum runmodes{
   ALS_MATRIX = 0,//alternating least squares, matrix (paper 3)
   BPTF_MATRIX = 1,//Bayesian prob. matrix factorization (paper 2)
   BPTF_TENSOR = 2,//Bayesian prob. tensor factorization (paper 1)
   BPTF_TENSOR_MULT = 3,
   ALS_TENSOR_MULT = 4,//alternating least squares, tensor
   SVD_PLUS_PLUS = 5, //SVD++, (paper 4)
   STOCHASTIC_GRADIENT_DESCENT = 6, //SGD (paper 5)
   LANCZOS = 7,// Lanczos algorithm (SVD) (reference 6)
   NMF = 8,
   WEIGHTED_ALS = 9,
   ALS_SPARSE_USR_FACTOR = 10,
   ALS_SPARSE_USR_MOVIE_FACTORS = 11,
   ALS_SPARSE_MOVIE_FACTOR = 12,
   SVD = 13, //simular value decompoistion via double sided Lanczos
   TIME_SVD_PLUS_PLUS = 14, //time-SVD++ (see reference 12)
   BIAS_SGD = 15,
   RBM = 16
};

#define MAX_RUNMODE 17

//counters for debugging running time of different modules
enum countervals{
   EDGE_TRAVERSAL=0,
   BPTF_SAMPLE_STEP=1,
   CALC_RMSE_Q=2,
   ALS_LEAST_SQUARES=3,
   BPTF_TIME_EDGES=4,
   BPTF_LEAST_SQUARES=5,
   CALC_OBJ =6,
   BPTF_MVN_RNDEX=7,
   BPTF_LEAST_SQUARES2=8, 
   SVD_MULT_A=9,
   SVD_MULT_A_TRANSPOSE=10,
};


/***
 * Different graph types for the different problems solved
 * 1) graph_type_vdpp - for running svd++
 * 2) graph_type_mcmc - for running mcmc methods
 * 3) graph_type_mult_edges - for having support for multiple edges in different times. model can support multiple ratings of user to the same movie in different times or a single rating. Single rating will run faster.
 * 4) graph_type - all other algorithms
 */
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;

typedef graphlab::graph<vertex_data, multiple_edges> graph_type_mult_edge;
typedef graphlab::types<graph_type_mult_edge> gl_types_mult_edge;

typedef graphlab::graph<vertex_data_svdpp, edge_data> graph_type_svdpp;
typedef graphlab::types<graph_type_svdpp> gl_types_svdpp;

typedef graphlab::graph<vertex_data, edge_data_mcmc> graph_type_mcmc;
typedef graphlab::types<graph_type_mcmc> gl_types_mcmc;



struct timesvdpp_output{
  mat ptemp;
  mat x;
  mat pu;
  mat q;
  mat z;
  mat pt;
};

/***
 * Class for  holding current information about the problem run
 */
class problem_setup{
public:

  bool BPTF; //is this a sampling MCMC algo?
  double pT; //regularization for tensor time nodes
  runmodes algorithm; //type of algorithm
 
  graphlab::timer gt;
  int iiter;//count number of time zero node run
  mat U,V,T; //for storing the output of ALS/BPTF

  vec svdpp_usr_bias;//for storing the output of SVD++
  vec svdpp_movie_bias; 
  
  timesvdpp_output timesvdpp_out; //for storing the ouptut of timesvd++

  mat dp;

/* Variables for PMF */
  int M,N,K,L;//training size: users, movies, times, number of edges
  int Le; //number of ratings in validation dataset 
  int Lt;//number of rating in test data set
  int Lt2;
  int last_node; //index of last node
  bool tensor; //is this tensor or a matrix
  bool isals; //is this algorithm an ALS variant
  double globalMean[3]; //store global mean of matrix/tensor entries

//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
 
  vertex_data * times;
  vertex_data_svdpp * times_svdpp;
  graphlab::core_base* glcore;
  graph_type * gg[4];
  graph_type_mcmc * g_mcmc[4];
  graph_type_mult_edge * g_mult_edge[4];
  graph_type_svdpp * g_svdpp[4];

  vec vones;
  mat eDT; 
  double pU; //regularization for users
  double pV; //regularization for movies

  double validation_rmse; //stores validation rmse 
  double training_rmse; //stored training rmse
  template<typename graph_type> const graph_type* g(testtype type);
    
  //template<typename graph_type> void set_graph(graph_type *g, testtype type);
  void set_graph(graph_type*_g, testtype type){gg[type]=_g;};
  void set_graph(graph_type_mcmc*_g, testtype type){g_mcmc[type]=_g;};
  void set_graph(graph_type_svdpp*_g, testtype type){g_svdpp[type]=_g;};
  void set_graph(graph_type_mult_edge*_g, testtype type){g_mult_edge[type]=_g;};

  
 problem_setup(){

  BPTF = false; //is this a sampling MCMC algo?
  pT = 1; //regularization for tensor time nodes
  algorithm = ALS_MATRIX; //type of algorithm
  iiter = 1;//count number of time zero node run

 /* Problem size */
  M=N=K=L=0;//training size: users, movies, times, number of edges
  Le = 0; //number of ratings in validation dataset 
  Lt = 0;//number of rating in test data set

  pU = pV = 0;
  memset(globalMean,0,3*sizeof(double));  //store global mean of matrix/tensor entries

  validation_rmse = 0; 
  training_rmse = 0;

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
  times = NULL;
  times_svdpp = NULL;
  glcore = NULL;
  //engine = NULL;
  memset(gg, 0, 3*sizeof(graph_type*));
  memset(g_mcmc, 0, 3*sizeof(graph_type_mcmc*));
  memset(g_mult_edge, 0, 3*sizeof(graph_type_mult_edge*));
  memset(g_svdpp, 0, 3*sizeof(graph_type_svdpp*));
}

  void verify_setup();
  bool to_print(int id);
};

extern advanced_config ac;

bool problem_setup::to_print(int id){
  return (ac.debug && (id == 0 || id == M-1 || id == M || M+N));
}

void problem_setup::verify_setup(){
  switch(algorithm){
  // iterative matrix factorization using alternating least squares
  // or SVD ++
  case ALS_MATRIX:
  case ALS_SPARSE_USR_FACTOR:
  case ALS_SPARSE_USR_MOVIE_FACTORS:
  case ALS_SPARSE_MOVIE_FACTOR:
  case WEIGHTED_ALS:
  case SVD_PLUS_PLUS:
  case STOCHASTIC_GRADIENT_DESCENT:
  case LANCZOS:
  case NMF:
  case SVD:
  case BIAS_SGD:
  case RBM:
    tensor = false; BPTF = false;
    break;

    // MCMC tensor factorization
  case BPTF_TENSOR:
    // tensor factorization , allow for multiple edges between user and movie in different times
  case BPTF_TENSOR_MULT:
    tensor = true; BPTF = true;
   break;
    //MCMC matrix factorization
  case BPTF_MATRIX:
    tensor = false; BPTF = true;
    break;
   // tensor factorization
  case ALS_TENSOR_MULT:
  case TIME_SVD_PLUS_PLUS:
    tensor = true; BPTF = false;
    break;
    
  default:
    assert(0);
  }

  if (algorithm == ALS_TENSOR_MULT || algorithm == ALS_MATRIX || algorithm == ALS_SPARSE_USR_FACTOR || algorithm == ALS_SPARSE_USR_MOVIE_FACTORS || algorithm == ALS_SPARSE_MOVIE_FACTOR || algorithm == WEIGHTED_ALS)
   isals = true;
  else isals = false;
 

  if (ac.bptf_delay_alpha != 0 && (algorithm != BPTF_TENSOR_MULT && algorithm != BPTF_TENSOR))
	logstream(LOG_WARNING) << "Delaying alpha (sampling of noise level) is ignored in non-MCMC methods" << std::endl;

  if (ac.bptf_burn_in != 10 && (algorithm != BPTF_TENSOR_MULT && algorithm != BPTF_TENSOR && algorithm != BPTF_MATRIX))
	logstream(LOG_WARNING) << "Markov chain burn in period is ignored in non-MCMC methods" << std::endl;

  if (ac.user_sparsity < 0.5 || ac.user_sparsity >= 1){
	logstream(LOG_FATAL) << "user_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
  }
  if (ac.movie_sparsity < 0.5 || ac.movie_sparsity >= 1){
	logstream(LOG_FATAL) << "movie_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
  }
    if (ac.minval == -DEF_MAX_VAL || ac.maxval == DEF_MAX_VAL){
      if (ac.algorithm == SVD_PLUS_PLUS || ac.algorithm == TIME_SVD_PLUS_PLUS)
         logstream(LOG_FATAL)<<"For SVD++/time-SVD++ you are required to specify min and max allowed matrix values using the flags --minval=XX, --maxval=XX. " << std::endl;
      else 
         logstream(LOG_WARNING)<<"It is recommended to set min and max allowed matrix values to improve prediction quality, using the flags --minval=XX, --maxval=XX" << std::endl;
    }

   if (ac.K > 0)
      if (tensor && ac.matrixmarket)
         K = ac.K; //set the number of time bins via command line for tensor read from matrix market file

   if (!ac.unittest){
		 if (boost::starts_with(ac.scheduler, "round_robin") && ac.algorithm == NMF)
       logstream(LOG_FATAL)<<"For NMF please do not specify a scheduler using the command --scheduler" << std::endl;
     else if (!ac.stats && !boost::starts_with(ac.scheduler, "round_robin") && (ac.algorithm == BIAS_SGD || isals || ac.algorithm == RBM ||
           ac.algorithm == STOCHASTIC_GRADIENT_DESCENT || ac.algorithm == TIME_SVD_PLUS_PLUS || ac.algorithm == SVD_PLUS_PLUS || tensor)){
       logstream(LOG_FATAL)<<"Please use round robin scheduler for this algorithm using the command line: --scheduler=\"round_robin(max_iteration=XX,block_size=1)\" XX is the number of required iterations" << std::endl;

     }
   }

   if (tensor && ac.matrixmarket && ac.matrixmarkettokensperrow == 3){
     logstream(LOG_WARNING)<<"When running tensor based factorization, with matrix market input format, number of matrix market tokens per row should be 4 (each row has [from] [to] [val] [time]\\n format). Use the command line argument --matrix_market_tokens_per_row=4 to avoid this warning" << std::endl;
     ac.matrixmarkettokensperrow = 4;
   }

}

template<> const graph_type *problem_setup::g(testtype type){ return gg[type]; }
template<> const graph_type_mcmc *problem_setup::g(testtype type){ return g_mcmc[type]; }
template<> const graph_type_mult_edge *problem_setup::g(testtype type){ return g_mult_edge[type]; }
template<> const graph_type_svdpp *problem_setup::g(testtype type){ return g_svdpp[type]; }


 

/**
 * functions forward declerations
 */
int do_main(int argc, const char * argv[]);
void set_num_edges(int val, testtype data_type);
void load_matrix_market(const char * filename, graph_type * _g, testtype data_type);
void verify_size(testtype data_type, int M, int N, int K);
void save_matrix_market_format(const char * filename, mat &U, mat &V);
float predict(const vertex_data& v1, const vertex_data & v2, const edge_data * edge, float rating, float & prediction);
float predict(const vertex_data_svdpp& user, const vertex_data_svdpp& movie, const edge_data * edge, const vertex_data * nothing, float rating, float & prediction);
float predict(const vertex_data& v1, const vertex_data& v2, const edge_data * edge, const vertex_data *v3, float rating, float &prediction);
void rbm_update_function(gl_types_svdpp::iscope &scope, gl_types_svdpp::icallback &scheduler);
void rbm_update_function(gl_types::iscope &scope, gl_types::icallback &scheduler){ assert(false); }
void rbm_update_function(gl_types_mcmc::iscope &scope, gl_types_mcmc::icallback &scheduler) { assert(false); }
void rbm_update_function(gl_types_mult_edge::iscope &scope, gl_types_mult_edge::icallback &scheduler) { assert(false); }
void rbm_init(); 

#endif

