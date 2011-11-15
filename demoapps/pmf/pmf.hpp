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



/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU


For BPTF algorithm:
See algrithm description and explanation in: 

1) Liang Xiong, Xi Chen, Tzu-Kuo Huang, Jeff Schneider, Jaime
G. Carbonell, Temporal Collaborative Filtering with Bayesian
Probabilistic Tensor Factorization. In Proceedings of SIAM Data
Mining, 2010.

2) Salakhutdinov and Mnih, Bayesian Probabilistic Matrix Factorization
using Markov Chain Monte Carlo. in International Conference on Machine
Learning, 2008.

For ALS (Alternating least squares) see algorithm and explanation:

3) Yunhong Zhou, Dennis Wilkinson, Robert Schreiber and Rong
Pan. Large-Scale Parallel Collaborative Filtering for the Netflix
Prize. Proceedings of the 4th international conference on Algorithmic
Aspects in Information and Management. Shanghai, China pp. 337-348,
2008.

For SVD++ see algorithm and explanations:

4) Koren, Yehuda. "Factorization meets the neighborhood: a
multifaceted collaborative filtering model." In Proceeding of the 14th
ACM SIGKDD international conference on Knowledge discovery and data
mining, 426434. ACM,
2008. http://portal.acm.org/citation.cfm?id=1401890.1401944

For SGD, see algorhtm: 

5) Matrix Factorization Techniques for Recommender Systems Yehuda
Koren, Robert Bell, Chris Volinsky In IEEE Computer, Vol. 42,
No. 8. (07 August 2009), pp. 30-37.

6) Tikk, D. (2009). Scalable Collaborative Filtering Approaches for
Large Recommender Systems. Journal of Machine Learning Research, 10,
623-656.

For Lanczos algorithm (SVD) see:
7) http://en.wikipedia.org/wiki/Lanczos_algorithm

For NMF (non-negative matrix factorization) see:

8) Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative
Matrix Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.

For Weigted alternating least squares see:

9) Rong Pan, Yunhong Zhou, Bin Cao, Nathan N. Liu, Rajan Lukose,
Martin Scholz, and Qiang Yang. 2008. One-Class Collaborative
Filtering. In Proceedings of the 2008 Eighth IEEE International
Conference on Data Mining (ICDM '08). IEEE Computer Society,
Washington, DC, USA, 502-511.

For implicit collaborative filering see:

10) One-Class Collaborative Filtering by: Rong Pan, Yunhong Zhou, Bin
Cao, N. N. Liu, R. Lukose, M. Scholz, Qiang Yang. Data Mining, IEEE
International Conference on In Data Mining, 2008. ICDM '08.

For sparsity enforcing priors see:

11) Xi Chen, Yanjun Qi, Bing Bai, Qihang Lin and Jaime
Carbonell. Sparse Latent Semantic Analysis. In SIAM International
Conference on Data Mining (SDM), 2011.

SVD is implemented using two sided Lanczos, see:
12) http://en.wikipedia.org/wiki/Singular_value_decomposition

*/



#ifndef PMF_HPP
#define PMF_HPP

#include <graphlab.hpp>

#include <Eigen/Dense.hpp>
#include <Eigen/Sparse.hpp>
#include <Eigen/Cholesky.hpp>
#include <Eigen/Eigenvalues.hpp>










/**
 * Define linear algbebra types
 */
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi ivec;
typedef EIgen::SparseVector<double> sparse_vec;


/** Vertex and edge data types **/
struct als_vertex_data {
  vec pvec; //! vector of learned values 
  float rmse; //!root of mean square error
  //constructor
  vertex_data() : rmse(0) { }
  void save(graphlab::oarchive& arc) const { 
    arc << pvec.size();
    for(size_t i = 0; i < pvec.size(); ++i) arc << pvec(i);
    arc << rmse;
  }   
  void load(graphlab::iarchive& arc) {
    size_t size; 
    arc >> size; pvec.resize(size);
    for(size_t i = 0; i < pvec.size(); ++i) arc >> pvec(i);
    arc >> rmse;
  }
};

struct als_edge_data {
  float  weight;  //observation 
  edge_data() : weight(0)  { } 
  void save(graphlab::oarchive& archive) const { archive << weight; }     
  void load(graphlab::iarchive& archive) { archive >> weight; }
}; // end of edge data

typedef graphlab::graph<als_vertex_data, als_edge_data> als_graph_type;













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
  mat U,V,T; //for storing the output
  vec svdpp_usr_bias;
  vec svdpp_movie_bias;
  mat dp;

/* Variables for PMF */
  int M,N,K,L;//training size: users, movies, times, number of edges
  int Le; //number of ratings in validation dataset 
  int Lt;//number of rating in test data set
  int last_node; //index of last node
  bool tensor; //is this tensor or a matrix
  bool isals; //is this algorithm an ALS variant
  double globalMean[3]; //store global mean of matrix/tensor entries

//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
  
  vertex_data * times;
  void* glcore;
  graph_type * gg[3];
  graph_type_mcmc * g_mcmc[3];
  graph_type_mult_edge * g_mult_edge[3];
  graph_type_svdpp*g_svdpp[3];

  vec vones;
  mat eDT; 
  double pU; //regularization for users
  double pV; //regularization for movies


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

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
  times = NULL;
  glcore = NULL;
  //engine = NULL;
  memset(gg, 0, 3*sizeof(graph_type*));
  memset(g_mcmc, 0, 3*sizeof(graph_type_mcmc*));
  memset(g_mult_edge, 0, 3*sizeof(graph_type_mult_edge*));
  memset(g_svdpp, 0, 3*sizeof(graph_type_svdpp*));
}

  void verify_setup();

};

extern advanced_config ac;

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
	logstream(LOG_ERROR) << "user_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
        exit(1);
  }
  if (ac.movie_sparsity < 0.5 || ac.movie_sparsity >= 1){
	logstream(LOG_ERROR) << "movie_sparsity of factor matrix has to be in the range [0.5 1)" << std::endl;
        exit(1);
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


#endif

