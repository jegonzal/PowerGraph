#ifndef PMF_H__	 
#define PMF_H__

//#define NDEBUG
#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>

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

*/
#include <vector>
#define GL_NO_MULT_EDGES //comment this flag, if you want to have support for multiple edges in different times between the same user and movie
#define GL_NO_MCMC //comment this flag, if you want to have support for MCMC methods (BPTF)
//#define GL_SVD_PP //comment this flag, if you are not running svd++ algorithm

using namespace itpp;

#define DEF_MAX_VAL 1e100
#define DEF_MIN_VAL -1e100


float maxval = DEF_MAX_VAL; //max allowed matrix/tensor entry
float minval = DEF_MIN_VAL;

float min(float a, float b){ return a<b?a:b; }
float max(float a, float b){ return a>b?a:b; }



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


int D=20;         //number of features

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
  vertex_data(){
    pvec = zeros(D);
    rmse = 0;
    num_edges = 0;
#ifdef GL_SVD_PP
    bias =0;
    weight = zeros(D);
#endif
  }

  void save(graphlab::oarchive& archive) const {  
    archive << pvec;
    archive << rmse << num_edges; 
#ifdef GL_SVD_PP
    archive << bias << weight;
#endif
  }  
   
  void load(graphlab::iarchive& archive) {  
     archive >> pvec;
     archive >> rmse >> num_edges;  
#ifdef GL_SVD_PP
     archive >> bias >> weight;
#endif
  }
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

//methods to compute the Root mean square error (RMSE)     
inline float predict(const vec& x1, const vec& x2, const edge_data * edge, float rating, float & prediction){
	prediction = dot(x1, x2);	
   //return the squared error
   prediction = min(prediction, maxval);
   prediction = max(prediction, minval);
	float sq_err = powf(prediction - rating, 2);
   return sq_err;
}

inline double predict(const vertex_data& user, const vertex_data &movie, const edge_data * edge, float rating, float & prediction){
#ifdef GL_SVD_PP	
   return svd_predict(user, movie, edge, rating, prediction);
#else
   return predict(user.pvec, movie.pvec, edge, rating, prediction);
#endif
} 


inline float predict(const vertex_data& v1, const vertex_data& v2, const edge_data * edge, const vertex_data *v3, float rating, float &prediction){
	if (v3 == NULL) //matrix	
		return predict(v1,v2,edge, rating,prediction);

	prediction = 0;
	for (int i=0; i< v1.pvec.size(); i++){
	   prediction += (v1.pvec[i] * v2.pvec[i] * v3->pvec.get(i));
	}
   prediction = min(prediction, maxval);
   prediction = max(prediction, minval);
   float sq_err = powf(prediction - rating, 2);
   return sq_err;
   
}


 
 double get_rmse(const vertex_data & v){
    return v.rmse;
 }



//data file types
enum testtype{
    TRAINING = 0,
    VALIDATION = 1,
    TEST = 2
};

const char * testtypename[] = {"TRAINING", "VALIDATION", "TEST"};

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
};

#define MAX_RUNMODE 9

const char * runmodesname[] = {"ALS_MATRIX (Alternating least squares)", "BPTF_MATRIX (Bayesian Prob. Matrix Factorization)", "BPTF_TENSOR (Bayesian Prob. Tensor Factorization)", "BPTF_TENSOR_MULT", "ALS_TENSOR_MULT", "SVD++", "SGD (Stochastic Gradient Descent)", "SVD (Singular Value Decomposition via LANCZOS)", "NMF (non-negative factorization)", "Weighted alternating least squares"};

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

const char * countername[] = {"EDGE_TRAVERSAL", "BPTF_SAMPLE_STEP", "CALC_RMSE_Q", "ALS_LEAST_SQUARES", \
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


#ifndef GL_SVD_PP
void svd_init(){};
void svd_plus_plus_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler){};
#else
void svd_init();
void svd_plus_plus_update_function(gl_types::iscope & scope, 
      gl_types::icallback & scheduler);
#endif

#endif

