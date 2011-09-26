#ifndef CLUSTERING_H__	 
#define CLUSTERING_H__

//#define NDEBUG
#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "graphlab.hpp"
#include "kcores.h"

using namespace itpp;


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
  sparse_vec datapoint;
  int current_cluster;
  int prev_cluster;
  float min_distance;
  bool reported;
  bool hot;
  vec distances;
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
  vec location;
  int num_assigned_points;
  vec cur_sum_of_points;
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
   KSHELL_DECOMPOSITION = 4
};

#define MAX_RUNMODE 1


enum initizliation_type{
   INIT_RANDOM = 0,
   INIT_ROUND_ROBIN = 1,
   INIT_KMEANS_PLUS_PLUS = 2,
   INIT_RANDOM_CLUSTER = 3
};


//counters for debugging running time of different modules
enum countervals{
   DISTANCE_CALCULATION=0,
   LDA_NEWTON_METHOD=1,
   LDA_ACCUM_BETA=2,
   LDA_LIKELIHOOD=3,
   LDA_NORMALIZE=4
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

  gl_types::core * glcore;
  gl_types_kcores::core * glcore_kcores;
  clusters clusts;

//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
  
  gl_types::iengine * engine;
  graph_type* gg;
  graph_type_kcores* g_kcores;

  mat output_clusters;
  mat output_assignements;
  int total_assigned;

  template<typename graph_type> graph_type* g();
  void set_graph(graph_type*_g){gg=_g;};
  void set_graph(graph_type_kcores*_g){g_kcores=_g;};
  void set_core(gl_types::core *_g){glcore=_g;};
  void set_core(gl_types_kcores::core *_g){glcore_kcores=_g;};
   
  problem_setup(){

  algorithm = K_MEANS; //type of algorithm
  iiter = 1;//count number of time zero node run

 /* Problem size */
  M=N=K=L=0;

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
  glcore = NULL;
  engine = NULL;
  gg = NULL;
  g_kcores = NULL;
  glcore_kcores = NULL;
  total_assigned = 0;
}

  void verify_setup();

};
template<> inline graph_type *problem_setup::g(){ return gg; }
template<> inline graph_type_kcores *problem_setup::g(){ return g_kcores; }

static graphlab::glshared<clusters> CLUSTERS;




void do_main(int argc, const char * argv[]);

void add_vertices(graph_type * _g);
void add_vertices(graph_type_kcores * _g);
void load_matrix_market(const char * filename, graph_type * _g);
void load_matrix_market(const char * filename, graph_type_kcores * _g);
void save_matrix_market_format(const char * filename);

void test_math();

void lda_main();
#endif

