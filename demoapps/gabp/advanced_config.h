#ifndef _ADVANCED_CONFIG
#define _ADVANCED_CONFIG



class advanced_config{

public:
  int D;
  bool debug;
  double threshold;
  int syncinterval;
  int iter;
  int ncpus;
  int algorithm;
  bool cg_resid;
  bool zero;
  bool isfloat;
  int unittest;
  bool square;
  std::string datafile;
  bool support_null_variance;
  bool round_robin;
  bool supportgraphlabcf;
  bool loadfactors; //start from a previous solution instead of a random point
  bool savefactors; //save solution to file
  bool loadgraph; //load graph from a binary saved file
  bool savegraph; //save input file into graph binary file
  bool stats; //print out statistics and exit
  bool regnormal; //regular normalization
  bool aggregatevalidation; //use validation dataset as training data
  bool outputvalidation; //experimental: output validation results of kdd format
  bool binaryoutput; //export the factors U,V,T to a binary file
  bool printhighprecision; //print RMSE output with high precision
  double scaling; //aggregate time values into bins? (default =1, no aggregation)
  double truncating ; // truncate unused time bins (optional, default = 0, no truncation)
  double scalerating; //scale the rating by dividing to the scalerating factor (optional)
#define DEF_MAX_VAL 1e100
  double minval; //minimal allowed value in matrix/tensor
  double maxval; //maximal allowed value in matrix/tensor
  
  bool mainfunc;
  bool manualgraphsetup; //if true, graph is loaded in user code and not by us

/* variables of BPTF */
double als_lambda;
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



  advanced_config(){
    D = 20;
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
    square = false;
    support_null_variance = false;
    round_robin = true;
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
    printhighprecision = false; //print RMSE output with high precision
    scaling = 1.0; //aggregate time values into bins? (default =1, no aggregation)
    truncating = 0.0; // truncate unused time bins (optional, default = 0, no truncation)
    scalerating = 1.0; //scale the rating by dividing to the scalerating factor (optional)
    minval = -DEF_MAX_VAL;
    maxval = DEF_MAX_VAL;

    mainfunc = true;
    manualgraphsetup = false;
    als_lambda = 1;
    bptf_delay_alpha = 0;
    bptf_burn_in = 10;
    bptf_alpha = 0;

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
 }


  graphlab::command_line_options init_command_line_options();
};


#endif //_ADVANCED_CONFIG
