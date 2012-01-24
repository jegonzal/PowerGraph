#ifndef CLUSTERING_H__	 
#define CLUSTERING_H__



//#define NDEBUG
#include "graphlab.hpp"
#include "../pmf/mathlayer.hpp"
#include "mathlayerf.hpp"

#define USE_DOUBLE
#ifdef USE_DOUBLE
typedef double flt_dbl;
typedef vec flt_dbl_vec;
typedef sparse_vec sparse_flt_dbl_vec;
typedef mat flt_dbl_mat;
#define FOR_ITERATOR_ FOR_ITERATOR
#define _randu randu
#define vec2fvec(a) a
#define fvec2vec(a) a
#define fmat2mat(a) a
#define mat2fmat(a) a
#else
typedef float flt_dbl;
typedef fvec flt_dbl_vec;
typedef sparse_fvec sparse_flt_dbl_vec;
typedef fmat flt_dbl_mat;
#ifdef HAS_EIGEN
#define FOR_ITERATOR_ FOR_ITERATOR2
#else
#define FOR_ITERATOR_ FOR_ITERATOR
#endif
#define zeros fzeros
#define init_vec init_fvec
#define _randu frandu
#endif

#include "kcores.h"

void knn_main();

//structs for holding edge data in file

struct edge_float{
  float from;
  float to;
  float weight;
};

struct edge_float_cf{
  float from;
  float to;
  float time;
  float weight;
};
struct edge_double_cf{
  int from;
  int to;
  double time;
  double weight;
};





/** Vertex and edge data types **/
struct vertex_data {
  sparse_flt_dbl_vec datapoint;
  int current_cluster;
  int prev_cluster;
  flt_dbl min_distance;
  bool reported;
  bool hot;
  flt_dbl_vec distances;
  bool clusterhead;

  //constructor
  vertex_data(){
    current_cluster = -1;
    prev_cluster = -1;
    min_distance = 0;
    reported = false;
    hot = false;
    clusterhead = false;
  }

  void save(graphlab::oarchive& archive) const; 
   
  void load(graphlab::iarchive& archive); 
};

struct edge_data {
  edge_data(){ 
  }


 void save(graphlab::oarchive& archive) const {  
//TODO  
}  
   
  void load(graphlab::iarchive& archive) {  
  //TODO
  }
};


struct cluster{
  flt_dbl_vec location;
  int num_assigned_points;
  flt_dbl_vec cur_sum_of_points;
  double sum_sqr;
  cluster(){
    num_assigned_points = 0;
    sum_sqr = 0;
  }
  
};

struct clusters{
  std::vector<cluster> cluster_vec;
};


//run modes
enum runmodes{
   K_MEANS = 0,//K-means algo
   K_MEANS_PLUS_PLUS = 1, //initalization for K_means
   K_MEANS_FUZZY = 2,
   LDA = 3,
   KSHELL_DECOMPOSITION = 4,
   ITEM_KNN = 5,
   USER_KNN = 6, 
   SVD_EXPERIMENTAL = 7
};

//#define MAX_RUNMODE 1

//data file types
enum testtype{
    TRAINING = 0,
    VALIDATION = 1,
    TEST = 2
};



enum initizliation_type{
   INIT_RANDOM = 0,
   INIT_ROUND_ROBIN = 1,
   INIT_KMEANS_PLUS_PLUS = 2,
   INIT_RANDOM_CLUSTER = 3,
   INIT_FROM_FILE = 4
};


//counters for debugging running time of different modules
enum countervals{
   DISTANCE_CALCULATION=0,
   LDA_NEWTON_METHOD=1,
   LDA_ACCUM_BETA=2,
   LDA_LIKELIHOOD=3,
   LDA_NORMALIZE=4,
   SVD_MULT_A=5,
   SVD_MULT_A_TRANSPOSE=6,
   CALC_RMSE_Q=7
};

enum svd_steps{
  SVD_FULL_PASS=0,
  SVD_ITERATIONS=1,
  SVD_EIGENDECOMPOTISION=2
};

typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;

typedef graphlab::graph<kcores_data, kcores_edge> graph_type_kcores;
typedef graphlab::types<graph_type_kcores> gl_types_kcores;

class problem_setup{
public:

  runmodes algorithm; //type of algorithm
  initizliation_type init_type; 
  graphlab::timer gt;
  int iiter;//count number of time zero node run

  int M,N; //data size MxN matrix
  int K;//number of clusters
  int L;//number of non zero elements in data
  int M_validation, M_test; //data points in validation and tests data
  int N_validation, N_test; //dimensions of validation /tets

  gl_types::core * glcore;
  gl_types_kcores::core * glcore_kcores;
  clusters clusts;

  double cost;
  double last_cost;
  double first_cost;
//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
  
  gl_types::iengine * engine;
  graph_type* gg;
  graph_type * test_graph;
  graph_type * validation_graph;
  graph_type_kcores* g_kcores;
  
  flt_dbl_mat output_clusters;
  flt_dbl_mat output_assignements;
  flt_dbl_mat V,U;
  bool output_assignements_integer;

  //explanatory comments to write into output files
  std::string output_assignements_comment;
  std::string output_clusters_comment;
  std::string V_comment, U_comment, output_comment3, output_comment4;

  flt_dbl_mat T; //for SVD_experimental
  int total_assigned;

  template<typename graph_type> graph_type* g();
  template<typename graph_type> graph_type* g(testtype type);
  void set_graph(graph_type*_g){gg=_g;};
  void set_graph(testtype type, graph_type*_g){
     switch(type){
       case TRAINING: gg=_g; break;
       case VALIDATION: validation_graph = _g; break;
       case TEST: test_graph = _g; break;
     }
  }
  void set_graph(testtype type, graph_type_kcores*_g){};
  void set_graph(graph_type_kcores*_g){g_kcores=_g;};
  void set_core(gl_types::core *_g){glcore=_g;};
  void set_core(gl_types_kcores::core *_g){glcore_kcores=_g;};
   
  problem_setup(){

  algorithm = K_MEANS; //type of algorithm
  iiter = 1;//count number of time zero node run

  output_assignements_integer = false;

 /* Problem size */
  M=N=K=L=0;
  M_validation = M_test = 0;

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
  glcore = NULL;
  engine = NULL;
  gg = NULL;
  g_kcores = NULL;
  glcore_kcores = NULL;
  total_assigned = 0;
  cost = last_cost = first_cost = 0;
}

  void verify_setup();

};
template<> inline graph_type *problem_setup::g(){ return gg; }
template<> inline graph_type_kcores *problem_setup::g(){ return g_kcores; }
template<> inline graph_type_kcores *problem_setup::g(testtype type){ return g_kcores; }
template<> inline graph_type *problem_setup::g(testtype type){ 
   switch(type){
      case TRAINING: return gg;
      case VALIDATION: return validation_graph;
      case TEST: return test_graph;
   }
   return NULL;
}


static graphlab::glshared<clusters> CLUSTERS;




int do_main(int argc, const char * argv[]);

void add_vertices(graph_type * _g, testtype type);
void add_vertices(graph_type_kcores * _g, testtype type);
void load_matrix_market(const char * filename, graph_type * _g, testtype type);
void load_matrix_market(const char * filename, graph_type_kcores * _g, testtype type);
void save_matrix_market_format(const char * filename);

void test_math();
void test_fmath();

void lda_main();
#endif

