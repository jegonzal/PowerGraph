#ifndef CLUSTERING_H__	 
#define CLUSTERING_H__

//#define NDEBUG
#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "graphlab.hpp"


using namespace itpp;


//structs for holding edge data in file

struct edge_float{
  float from;
  float to;
  float time;
  float weight;
};



/** Vertex and edge data types **/
struct vertex_data {
  sparse_vec datapoint;
  int current_cluster;
  float min_distance;
  //constructor
  vertex_data(){
    current_cluster = -1;
    min_distance = 0;
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




//run modes
enum runmodes{
   K_MEANS = 0//K-means algo
};

#define MAX_RUNMODE 1

static const char * runmodesname[] = {"K-means"};


//counters for debugging running time of different modules
enum countervals{
   DISTANCE_CALCULATION
};

static const char * countername[] = {"DISTANCE_CALCULTION"};



typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;



class problem_setup{
public:

  runmodes algorithm; //type of algorithm
 
  graphlab::timer gt;
  int iiter;//count number of time zero node run

  int M,N; //data size MxN matrix
  int K;//number of clusters
  int L;//number of non zero elements in data

  gl_types::core * glcore;

//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
  
  gl_types::iengine * engine;
  graph_type* g;

  mat output_clusters;
  vec output_assignements;

 problem_setup(){

  algorithm = K_MEANS; //type of algorithm
  iiter = 1;//count number of time zero node run

 /* Problem size */
  M=N=K=L=0;

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
  glcore = NULL;
  engine = NULL;
  g = NULL;
}

  void verify_setup();

};

static graphlab::glshared<mat> CLUSTER_LOCATIONS;




void do_main(int argc, const char * argv[]);

void add_vertices(graph_type * _g);
void load_matrix_market(const char * filename, graph_type * _g);
void save_matrix_market_format(const char * filename);

void test_math();
#endif

