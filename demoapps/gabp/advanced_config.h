#ifndef _ADVANCED_CONFIG
#define _ADVANCED_CONFIG

#include "graphlab.hpp"

class advanced_config{

public:
  bool debug;
  double threshold;
  int syncinterval;
  int iter;
  int ncpus;
  int algorithm;
  int init_mode; //initialization mode used for clustering
  bool cg_resid;
 
  bool omp_support;
 
  //marix proprties
  bool zero;//allow zero entries in matrix?
  bool isfloat;//input file in float format
  bool square;//is matrix square?
  int D; //width of factor matrices

  /* clustering related fields */
  int K; //number of clusters 
  bool clusterdump; //dump cluster locations into a text file
  bool tfidf; //deploy tf-idf transformation on matrix values

  int unittest;

  //input output related stuff
  std::string datafile;
  bool matrixmarket; //martix market input format
  bool loadfactors; //start from a previous solution instead of a random point
  bool savefactors; //save solution to file
  bool loadgraph; //load graph from a binary saved file
  bool savegraph; //save input file into graph binary file
  bool binaryoutput; //export the factors U,V,T to a binary file
   
  bool support_null_variance;
  //bool round_robin;
  bool supportgraphlabcf;
 bool stats; //print out statistics and exit
  bool regnormal; //regular normalization
  bool aggregatevalidation; //use validation dataset as training data
  bool outputvalidation; //experimental: output validation results of kdd format
  bool printhighprecision; //print RMSE output with high precision
  double scaling; //aggregate time values into bins? (default =1, no aggregation)
  double truncating ; // truncate unused time bins (optional, default = 0, no truncation)
  double scalerating; //scale the rating by dividing to the scalerating factor (optional)
#define DEF_MAX_VAL 1e100
  double minval; //minimal allowed value in matrix/tensor
  double maxval; //maximal allowed value in matrix/tensor
  std::string scheduler; 
  bool mainfunc;
  bool manualgraphsetup; //if true, graph is loaded in user code and not by us
  double regularization;

/* variables for ALS */
double als_lambda;

/* variables of BPTF */
int bptf_delay_alpha; //delay alpha sampling (optional, for BPTF)
int bptf_burn_in; //burn-in priod (for MCMC sampling - optional)
double bptf_alpha;

/* Variables for SVD++ */
float svdpp_step_dec;//step decrement size for SVD++

/* Variables for SGD */
float sgd_gamma; //step size
float sgd_lambda; //starting step size
float sgd_step_dec; //step decrement size

int svd_iter;

/* implicit ratings variables (see reference 10 in pmf.h) */
float implicitratingweight;
float implicitratingvalue;
std::string implicitratingtype;
float implicitratingpercentage;

/* sparsity enforcing priors (see reference 11 in pmf.h) */
int lasso_max_iter;
#define DEFAULT_SPARSITY 0.8
double user_sparsity;
double movie_sparsity;

bool FLOAT; //is data in float format

/* for LDA */
int em_max_inner_iter; //number of inner EM iterations

int distance_measure;

  advanced_config(){
    D = 20;
    K = 0;
    debug = true;
    threshold = 1e-10;
    syncinterval = 10000;
    iter = 20;
    ncpus = 1;
    cg_resid = true;
    zero=true;
    isfloat = true;
    unittest = 0;
    algorithm = 0;
    init_mode = 0;
    square = false;
    support_null_variance = false;
    //round_robin = true;
    supportgraphlabcf = false;
    debug = false; //debug mode
    zero = false; //support edges with zero weight
    loadfactors = false; //start from a previous solution instead of a random point
    savefactors = true; //save solution to file
    loadgraph = false; //load graph from a binary saved file
    savegraph = false; //save input file into graph binary file
    stats = false; //print out statistics and exit
    regnormal = false; //regular normalization
    aggregatevalidation = false; //use validation dataset as training data
    outputvalidation = false; //experimental: output validation results of kdd format
    binaryoutput = false; //export the factors U,V,T to a binary file
    matrixmarket = false; //is input/output in matrix market format
    printhighprecision = false; //print RMSE output with high precision
    scaling = 1.0; //aggregate time values into bins? (default =1, no aggregation)
    truncating = 0.0; // truncate unused time bins (optional, default = 0, no truncation)
    scalerating = 1.0; //scale the rating by dividing to the scalerating factor (optional)
    minval = -DEF_MAX_VAL;
    maxval = DEF_MAX_VAL;

    mainfunc = true;
    manualgraphsetup = false;

    omp_support = true;

    als_lambda = 1;
    bptf_delay_alpha = 0;
    bptf_burn_in = 10;
    bptf_alpha = 0;
    regularization = 0;
  /* Variables for SVD++ */
    svdpp_step_dec = 0.9;//step decrement size for SVD++

  /* Variables for SGD */
   sgd_gamma = 1e-2; //step size
   sgd_lambda = 0.3; //starting step size
   sgd_step_dec = 0.9; //step decrement size

  /* variables for SVD */
   svd_iter = 10; //number of iterations (which is the number of extracted eigenvectors)

  /* implicit ratings variables (see reference 10 in pmf.h) */
   implicitratingweight = 1;
   implicitratingvalue = 0;
   implicitratingtype = "none";
   implicitratingpercentage = 0;

  /* sparsity enforcing priors (see reference 11 in pmf.h) */
   lasso_max_iter = 10;
   user_sparsity = DEFAULT_SPARSITY;
   movie_sparsity= DEFAULT_SPARSITY;

   FLOAT = false;
   als_lambda = 1;

   /* for clustering */
   distance_measure = 0;
   clusterdump = false;
   tfidf = false;

   /* for LDA */
   em_max_inner_iter = 20;
 }


  void init_command_line_options(graphlab::command_line_options & clopts);
};


#endif //_ADVANCED_CONFIG
