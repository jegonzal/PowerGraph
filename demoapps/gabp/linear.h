#ifndef _LINEAR_H
#define _LINEAR_H

#include <graphlab.hpp>
#ifdef HAS_ITPP
#include <itpp/itbase.h>
#endif

//define sdouble as either float or double as needed
typedef double sdouble;


enum constant_offsets {GABP_PRIOR_MEAN_OFFSET = 0, //prior mean (b_i / A_ii)
                       GABP_PRIOR_PREC_OFFSET = 1, //prior precision P_ii = A_ii
                       GABP_REAL_OFFSET = 2, // the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
                       GABP_CUR_MEAN_OFFSET = 3, //intermediate value to store the mean \mu_i
                       GABP_CUR_PREC_OFFSET = 4, //current precision value P_i
                       GABP_PREV_MEAN_OFFSET = 5,// mean value from previous round (for convergence detection)
                       GABP_PREV_PREC_OFFSET = 6}; // precision value from previous round (for convergence detection)




/** Vertex data types for linear solvers**/
struct vertex_data {
  sdouble prior_mean;  //prior mean (b_i / A_ii)
  sdouble prior_prec;   //prior precision P_ii = A_ii
  sdouble real;   //the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
  sdouble cur_mean; //intermediate value to store the mean \mu_i
  sdouble cur_prec; // //current precision value P_i
  sdouble prev_mean; //  mean value from previous round (for convergence detection)
  sdouble prev_prec; //precision value from previous round (for convergence detection)

  vertex_data(){ 
     prev_mean = 1000;
     prior_mean = prior_prec = real = cur_mean = cur_prec = prev_prec = 0;
   };

};

/** Vertex data type for matrix inverse **/
struct vertex_data_inv {
  sdouble * prior_mean;  //prior mean (b_i / A_ii)
  sdouble prior_prec;   //prior precision P_ii = A_ii
  sdouble * real;   //the real solution (if known) x = inv(A)*b for square matrix or x = pinv(A) for skinny matrix
  sdouble * cur_mean; //intermediate value to store the mean \mu_i
  sdouble cur_prec; // //current precision value P_i
  sdouble * prev_mean; //  mean value from previous round (for convergence detection)
  sdouble prev_prec; //precision value from previous round (for convergence detection)

  vertex_data_inv(){ 
     prev_mean = NULL;
     cur_mean = NULL;
     prior_mean = NULL;
     real = NULL;
     prior_prec = cur_prec = prev_prec = 0;
   };

};

/** edge data type **/
//edge is a scalar  non zero entry A_{ij} in the matrix A (row i, col j)
struct edge_data {
  sdouble weight; //edge value
  sdouble mean; // message \mu_ij
  sdouble prec; // message P_ij

  edge_data(){
     weight = mean = prec = 0;
  }

};

/** edge data type **/
//edge is a scalar  non zero entry A_{ij} in the matrix A (row i, col j)
struct edge_data_inv {
  sdouble weight; //edge value
  sdouble *mean; // message \mu_ij
  sdouble prec; // message P_ij

  edge_data_inv(){
     weight = prec = 0;
     mean= NULL;
  }

};



struct vertex_data_shotgun{
  bool active;
  double x;
  double xjneg;
  double y;
  double Ax;
  double expAx;
#ifdef HAS_ITPP
  itpp::sparse_vec features;
#endif


  vertex_data_shotgun(){
    active = true;
    x = y = xjneg = Ax = 0;
    expAx = 1;
  }
};

struct edge_data_shotgun{
};

enum runmodes{
   GaBP = 0, 
   JACOBI = 1,
   CONJUGATE_GRADIENT = 2,
   GaBP_INV = 3,
   LEAST_SQUARES = 4,
   SHOTGUN_LASSO = 5,
   SHOTGUN_LOGREG = 6
};



//counters for debugging running time of different modules
enum countervals{
   EDGE_TRAVERSAL=0,
   NODE_TRAVERSAL=1,
   RECOMPUTE_EXP_AX_LOGREG=2
};

class problem_setup{
public:

  runmodes algorithm; //type of algorithm
  graphlab::timer gt;
  int iiter;//count number of time zero node run

  uint32_t m; // number of rows of A
  uint32_t n; // number of cols of A
  uint32_t e; // number of edges

//performance counters
#define MAX_COUNTER 20
  double counter[MAX_COUNTER];
  

  //for shotgun logreg
  bool cdn_all_zero; //we have reached an all zero solution
  int cdn_neg_y; //number of negative training instances
  int cdn_pos_y; //number of positive training instances
  unsigned long long int shotgun_numshoots;
  //vectors for storing the output
  std::vector<double> means;
  std::vector<double> prec;

  int last_node;

 problem_setup(){

  iiter = 1;//count number of time zero node run
 /* Problem size */

 cdn_all_zero = false;
  cdn_neg_y = cdn_pos_y = 0;
  shotgun_numshoots = 0;

//performance counters
  memset(counter, 0, MAX_COUNTER*sizeof(double));
}

  void verify_setup(graphlab::command_line_options & clopts);

};
typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;
typedef graphlab::graph<vertex_data_inv, edge_data_inv> graph_type_inv;
typedef graphlab::types<graph_type_inv> gl_types_inv;
typedef graphlab::graph<vertex_data_shotgun, edge_data_shotgun> graph_type_shotgun;
typedef graphlab::types<graph_type_shotgun> gl_types_shotgun;


#endif

