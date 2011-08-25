#ifndef PMF_H__	 
#define PMF_H__

//#define NDEBUG
#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "graphlab.hpp"
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

For SGD, see algorhtm:
5) Matrix Factorization Techniques for Recommender Systems
Yehuda Koren, Robert Bell, Chris Volinsky
In IEEE Computer, Vol. 42, No. 8. (07 August 2009), pp. 30-37. 
6) Tikk, D. (2009). Scalable Collaborative Filtering Approaches for Large Recommender Systems. Journal of Machine Learning Research, 10, 623-656.

For Lanczos algorithm (SVD) see:
7) http://en.wikipedia.org/wiki/Lanczos_algorithm 

For NMF (non-negative matrix factorization) see:
8) Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.

For Weigted alternating least squares see:
9) Rong Pan, Yunhong Zhou, Bin Cao, Nathan N. Liu, Rajan Lukose, Martin Scholz, and Qiang Yang. 2008. One-Class Collaborative Filtering. In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM '08). IEEE Computer Society, Washington, DC, USA, 502-511. 

For implicit collaborative filering see:
10) One-Class Collaborative Filtering
by: Rong Pan, Yunhong Zhou, Bin Cao, N. N. Liu, R. Lukose, M. Scholz, Qiang Yang. Data Mining, IEEE International Conference on In Data Mining, 2008. ICDM '08.

For sparsity enforcing priors see:
11) Xi Chen, Yanjun Qi, Bing Bai, Qihang Lin and Jaime Carbonell. Sparse Latent Semantic Analysis. In SIAM International Conference on Data Mining (SDM), 2011.

*/
#include <vector>
#define GL_NO_MULT_EDGES //comment this flag, if you want to have support for multiple edges in different times between the same user and movie
#define GL_NO_MCMC //comment this flag, if you want to have support for MCMC methods (BPTF)
//#define GL_SVD_PP //comment this flag, if you are not running svd++ algorithm

using namespace itpp;


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

#ifdef GL_SVD_PP //data structure for svd++/svd only
  float bias; //bias for this user/movie
  vec weight; //weight vector for this user/movie
#endif

  //constructor
  vertex_data();

  void save(graphlab::oarchive& archive) const; 
   
  void load(graphlab::iarchive& archive); 
};

struct edge_data {
  float  weight;  //observation 
  float  time; //time of observation (for tensor algorithms)
#ifndef GL_NO_MCMC  
  float avgprd;
#endif
  edge_data(){ 
	weight = 0; 
	time = 0; 
#ifndef GL_NO_MCMC
	avgprd = 0;
#endif
}

 void save(graphlab::oarchive& archive) const {  
    archive << weight << time; 
#ifndef GL_NO_MCMC	
	archive<< avgprd;; 
#endif
  }  
   
  void load(graphlab::iarchive& archive) {  
    archive >> weight >> time;
#ifndef GL_NO_MCMC
     archive >> avgprd; 
#endif 
  }
};

//containiner for handling multiple edge
struct multiple_edges{
  std::vector<edge_data> medges;
 void save(graphlab::oarchive& archive) const {  
       archive << medges; 
  }  
   
  void load(graphlab::iarchive& archive) {  
      archive >> medges;  
  }
};


float svd_predict(const vertex_data& user, const vertex_data& movie, const edge_data * edge, float rating, float & prediction);
float predict(const vec& x1, const vec& x2, const edge_data * edge, float rating, float & prediction);
double predict(const vertex_data& user, const vertex_data &movie, const edge_data * edge, float rating, float & prediction);
float predict(const vertex_data& v1, const vertex_data& v2, const edge_data * edge, const vertex_data *v3, float rating, float &prediction);


 
 double get_rmse(const vertex_data & v);


//data file types
enum testtype{
    TRAINING = 0,
    VALIDATION = 1,
    TEST = 2
};

static const char * testtypename[] = {"TRAINING", "VALIDATION", "TEST"};


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
   ALS_SPARSE_MOVIE_FACTOR = 12
};

#define MAX_RUNMODE 9

static const char * runmodesname[] = {"ALS_MATRIX (Alternating least squares)", "BPTF_MATRIX (Bayesian Prob. Matrix Factorization)", "BPTF_TENSOR (Bayesian Prob. Tensor Factorization)", "BPTF_TENSOR_MULT", "ALS_TENSOR_MULT", "SVD++", "SGD (Stochastic Gradient Descent)", "SVD (Singular Value Decomposition via LANCZOS)", "NMF (non-negative factorization)", "Weighted alternating least squares", "Alternating least squares with sparse user factor matrix", "Alternating least squares with doubly sparse (user/movie) factor matrices", "Alternating least squares with sparse movie factor matrix"};

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

static const char * countername[] = {"EDGE_TRAVERSAL", "BPTF_SAMPLE_STEP", "CALC_RMSE_Q", "ALS_LEAST_SQUARES", \
  "BPTF_TIME_EDGES", "BPTF_LEAST_SQUARES", "CALC_OBJ", "BPTF_MVN_RNDEX", "BPTF_LEAST_SQUARES2", "SVD_MULT_A", "SVD_MULT_A_TRANSPOSE"};



//types of graph nodes
enum colors{
   COLOR_USER=0,
   COLOR_MOVIE=1,
   COLOR_TIME=2,
   COLOR_LAST=3,
};


//model can support multiple ratings of user to the same movie in different times
//or a single rating. Single rating will run faster.
#ifndef GL_NO_MULT_EDGES
typedef graphlab::graph<vertex_data, multiple_edges> graph_type;
#else
typedef graphlab::graph<vertex_data, edge_data> graph_type;
#endif
typedef graphlab::types<graph_type> gl_types;





double agg_rmse_by_movie(double & res);
double agg_rmse_by_user(double & res);

void svd_init();
void svd_plus_plus_update_function(gl_types::iscope & scope, 
      gl_types::icallback & scheduler);


class problem_setup{
public:

  bool BPTF; //is this a sampling MCMC algo?
  double pT; //regularization for tensor time nodes
  runmodes algorithm; //type of algorithm
 
  graphlab::timer gt;
  int iiter;//count number of time zero node run
  mat U,V,T; //for storing the output
  mat dp;

/* Variables for PMF */
  int M,N,K,L;//training size: users, movies, times, number of edges
  int Le; //number of ratings in validation dataset 
  int Lt;//number of rating in test data set

  bool tensor; //is this tensor or a matrix

  double globalMean[3]; //store global mean of matrix/tensor entries

//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
  
  vertex_data * times;
  gl_types::iengine * engine;
  graph_type* g;
  graph_type validation_graph;
  graph_type test_graph;

 problem_setup(){

  BPTF = false; //is this a sampling MCMC algo?
  pT = 1; //regularization for tensor time nodes
  algorithm = ALS_MATRIX; //type of algorithm
  iiter = 1;//count number of time zero node run

 /* Problem size */
  M=N=K=L=0;//training size: users, movies, times, number of edges
  Le = 0; //number of ratings in validation dataset 
  Lt = 0;//number of rating in test data set

  memset(globalMean,0,3*sizeof(double));  //store global mean of matrix/tensor entries

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
  times = NULL;

  engine = NULL;
  g = NULL;
}

  void verify_setup();

};

void do_main(int argc, const char * argv[]);


#endif

