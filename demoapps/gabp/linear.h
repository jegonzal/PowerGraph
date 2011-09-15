#ifndef _LINEAR_H
#define _LINEAR_H

#include <graphlab.hpp>


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


enum algorithms{
   GaBP = 0, 
   JACOBI = 1,
   CONJUGATE_GRADIENT = 2,
   GaBP_INV = 3,
   LEAST_SQUARES = 4,
};

const char* algorithmnames[]= {"GaBP", "Jacobi", "Conjugate Gradient", "GaBP inverse", "Least Squares"};


//counters for debugging running time of different modules
enum countervals{
   EDGE_TRAVERSAL=0,
   NODE_TRAVERSAL=1
};

#define MAX_COUNTER 2
const char * countername[] = {"EDGE_TRAVERSAL", "NODE_TRAVERSAL"};


typedef graphlab::graph<vertex_data, edge_data> graph_type;
typedef graphlab::types<graph_type> gl_types;
graphlab::glshared<double> REAL_NORM_KEY;
graphlab::glshared<double> RELATIVE_NORM_KEY;
graphlab::glshared<size_t> ITERATION_KEY;
graphlab::glshared<double> THRESHOLD_KEY;
graphlab::glshared<bool> SUPPORT_NULL_VARIANCE_KEY;
graphlab::glshared<bool> ROUND_ROBIN_KEY;
graphlab::glshared<bool> DEBUG_KEY;
graphlab::glshared<size_t> MAX_ITER_KEY;


typedef graphlab::graph<vertex_data_inv, edge_data_inv> graph_type_inv;
typedef graphlab::types<graph_type_inv> gl_types_inv;
graphlab::glshared<int> MATRIX_WIDTH_KEY;



#endif

