/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
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
 */


/**
 * GRAPHLAB implementation of Gaussiabn Belief Propagation Code See
 * algrithm description and explanation in: Danny Bickson, Gaussian
 * Belief Propagation: Theory and Application. Ph.D. Thesis. The
 * Hebrew University of Jerusalem. Submitted October 2008.
 * http://arxiv.org/abs/0811.2518 By Danny Bickson, CMU. Send any bug
 * fixes/reports to bickson@cs.cmu.edu Code adapted to GraphLab by
 * Joey Gonzalez, CMU July 2010
 *
 * Functionality: The code solves the linear system Ax = b using
 * Gaussian Belief Propagation. (A is either square matrix or
 * skinny). A assumed to be full column rank.  Algorithm is described
 * in Algorithm 1, page 14 of the above Phd Thesis.
 *
 * If you are using this code, you should cite the above reference. Thanks!
 */

#define NO_CG_SUPPRT //comment this flag if you would like to compile cg code

#include <cmath>
#include <cstdio>
#include "linear.h"
#include "gabp.hpp"
#include "jacobi.hpp"
#include "cg.hpp"
#include "math.hpp"
#include "syncs.hpp"
#include "gabp_inv.hpp"
#include "advanced_config.h"

#ifdef HAS_ITPP
#include <itpp/itbase.h>
#endif

#include <graphlab/macros_def.hpp>

double counter[20];//timing counters

uint32_t m = 0; // number of rows of A
uint32_t n = 0; // number of cols of A
uint32_t e = 0; // number of edges
advanced_config config;

#define BUFSIZE 500000
void read_nodes(FILE * f, int len, int offset, int nodes,
                graph_type * g){

  assert(offset>=0 && offset < len);
  assert(nodes>0);

  int toread = nodes;
  if (nodes > BUFSIZE)
    toread = BUFSIZE;
  int remain = nodes;

  while(remain > 0){
    double * temp = new double[remain];
    toread = (remain < BUFSIZE)?remain:toread;
    int rc = (int)fread(temp, sizeof(double), toread, f);
    //assert(rc == toread);
    remain -= rc;
    vertex_data data;
    for (int i=0; i< rc; i++){
        if (offset == GABP_PRIOR_MEAN_OFFSET)
          data.prior_mean = temp[i];
        else if (offset == GABP_REAL_OFFSET)
          data.real = temp[i];
        else assert(false);
      g->add_vertex(data);
    }
    delete [] temp;
  }
}
void read_nodes(FILE * f, int len, int offset, int nodes,
                graph_type_inv * g){

  assert(offset>=0 && offset < len);
  assert(nodes>0);

  int toread = nodes;
  if (nodes > BUFSIZE)
    toread = BUFSIZE;
  int remain = nodes;

  while(remain > 0){
    double * temp = new double[remain];
    toread = (remain < BUFSIZE)?remain:toread;
    int rc = (int)fread(temp, sizeof(double), toread, f);
    //assert(rc == toread);
    remain -= rc;
    vertex_data_inv data;
    for (int i=0; i< rc; i++){
      data.prior_mean = new sdouble[nodes];
      memset(data.prior_mean, 0, nodes*sizeof(sdouble));
      data.prior_mean[i] = 1;
      data.cur_mean = new sdouble[nodes];
      memset(data.cur_mean, 0, nodes*sizeof(sdouble));
      data.prev_mean = new sdouble[nodes];
      memset(data.prev_mean, 1, nodes*sizeof(sdouble));
      g->add_vertex(data);
    }
    delete [] temp;
  }
}



//read a vector from file and return an array
double * read_vec(FILE * f, size_t len){
  double * vec = new double[len];
  assert(vec != NULL);
  fread(vec, len, sizeof(double), f);
  return vec;
}
//write an output vector to file
void write_vec(FILE * f, int len, double * array){
  assert(f != NULL && array != NULL);
  fwrite(array, len, sizeof(double), f);
}


template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset,
                  graph * g, double * vec, int len,
                  bool free_vec = false, int d = 1){
  assert(end_pos - start_pos > 1);
  assert(len == end_pos - start_pos);
  assert(vec != NULL);
  assert(offset>=0 && offset < len);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++){
      data[offset+d*j] = vec[d*(i-start_pos)+j];
    }
  }
  if (free_vec)
    delete [] vec;
}

template<typename graph>
void dispatch_vec(int start_pos, int end_pos, int offset,
                  graph * g,sdouble val,int d = 1){
  assert(end_pos - start_pos > 1);
  for (int i=start_pos; i < end_pos; i++){
    sdouble * data = (sdouble*)&g->vertex_data(i);
    for (int j=0; j<d; j++)
      data[offset+d*j] = val;
  }
}



//struct for holding edge data in file
struct edata{
  int from;
  int to;
  double weight;
};

struct edata2{
  float from;
  float to;
  float time;
  float weight;
};

struct edata3{
  int from;
  int to;
  double time;
  double weight;
};

//read edges from file into the graph
template<typename edatatype>
int read_edges(FILE * f, int len, int offset, int nodes,
               graph_type_inv * g, advanced_config & config, bool symmetry = false){
  assert(offset>=0 && offset < len);
  assert(nodes > 0);

  unsigned int e,g0;
  int rc = fread(&e,1,4,f);
  assert(rc == 4);
  if (!config.supportgraphlabcf){
    rc = fread(&g0,1,4,f); //zero pad
    assert(rc == 4);
    assert(g0 == 0);
  }
  printf("Creating %d edges...\n", e);
  assert(e>0);
  int total = 0;
  edatatype* ed = new edatatype[200000];
  printf("symmetry: %d\n", symmetry);
  int edgecount_in_file = e;
  if (symmetry) edgecount_in_file /= 2;
  while(true){
    memset(ed, 0, 200000*sizeof(edatatype));
    rc = (int)fread(ed, sizeof(edatatype),
                    std::min(200000, edgecount_in_file - total), f);
    total += rc;

    for (int i=0; i<rc; i++){
      //memset(tmp, 0, len/sizeof(sdouble));
      edge_data_inv tmp;
      tmp.mean = new sdouble[nodes];
      memset(tmp.mean, 0, sizeof(sdouble)*nodes);
      tmp.weight =  ed[i].weight;
      if (!config.zero)
        assert(ed[i].weight != 0);
      assert(ed[i].from >= 1 && ed[i].from <= nodes);
      assert(ed[i].to >= 1 && ed[i].to <= nodes);
      assert(ed[i].to != ed[i].from);
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(ed[i].from-1, ed[i].to-1, tmp);
      if (!config.square) { //add the reverse edge as well
        // Matlab export has ids starting from 1, ours start from 0
        g->add_edge(ed[i].to-1, ed[i].from-1, tmp);
      }
    }
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  delete [] ed; ed = NULL;
  return e;
}


//read edges from file into the graph
template<typename edatatype>
int read_edges(FILE * f, int len, int offset, int nodes,
               graph_type * g, advanced_config & config, bool symmetry = false){
  assert(offset>=0 && offset < len);


  unsigned int e,g0;
  int rc = fread(&e,1,4,f);
  assert(rc == 4);
  if (!config.supportgraphlabcf){
    rc = fread(&g0,1,4,f); //zero pad
    assert(rc == 4);
    assert(g0 == 0);
  }
  printf("Creating %d edges...\n", e);
  assert(e>0);
  int total = 0;
  edatatype* ed = new edatatype[200000];
  printf("symmetry: %d\n", symmetry);
  int edgecount_in_file = e;
  if (symmetry) edgecount_in_file /= 2;
  while(true){
    memset(ed, 0, 200000*sizeof(edatatype));
    rc = (int)fread(ed, sizeof(edatatype),
                    std::min(200000, edgecount_in_file - total), f);
    total += rc;

    edge_data tmp;
    for (int i=0; i<rc; i++){
      //memset(tmp, 0, len/sizeof(sdouble));
      tmp.weight =  ed[i].weight;
      if (!config.zero)
        assert(ed[i].weight != 0);
      assert(ed[i].from >= 1 && ed[i].from <= nodes);
      assert(ed[i].to >= 1 && ed[i].to <= nodes);
      assert(ed[i].to != ed[i].from);
      // Matlab export has ids starting from 1, ours start from 0
      g->add_edge(ed[i].from-1, ed[i].to-1, tmp);
      if (!config.square) { //add the reverse edge as well
        // Matlab export has ids starting from 1, ours start from 0
        g->add_edge(ed[i].to-1, ed[i].from-1, tmp);
      }
    }
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  delete [] ed; ed = NULL;
  return e;
}

FILE * load_matrix_metadata(const char * filename){
   printf("Loading %s\n", filename);
   FILE * f = fopen(filename, "r");
   assert(f!= NULL);

   fread(&m, 1, 4, f);
   fread(&n, 1, 4, f);
   if (n == 0) 
      n=m; //compatability with older file format, will be removed later
   return f;
}

/*
 *  READ A SQUARE INVERSE COV MATRIX A of size nxn
 *  Where the main digonal is the precision vector
 * */

void load_square_matrix(FILE * f, graph_type& graph, advanced_config & config) {

  assert(m == n);
  assert(n > 0);
  //read the observation vector y of size n
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble), GABP_PRIOR_MEAN_OFFSET,
             n, &graph);

  //read the real solution of size n (if given, otherwise it is zero)
  double * real = read_vec(f, n);
  dispatch_vec(0,n,GABP_REAL_OFFSET, &graph, real, n, true);

  //read the precition of size n (the main diagonal of the matrix A)
  double * prec = read_vec(f, n);
  dispatch_vec(0,n,GABP_PRIOR_PREC_OFFSET, &graph, prec, n, true);

  e = read_edges<edata>(f, sizeof(edge_data)/sizeof(sdouble), 0, n, &graph, config);
  fclose(f);
}

/*
 *  READ A SQUARE INVERSE COV MATRIX A of size nxn
 *  Where the main digonal is the precision vector
 * */

void load_square_matrix(FILE * f, graph_type_inv& graph, advanced_config & config) {

  assert(m == n);
  assert(n > 0);
  //read the observation vector y of size n
  read_nodes(f, sizeof(vertex_data)/sizeof(sdouble), GABP_PRIOR_MEAN_OFFSET,
             n, &graph);

  //unused, skip this vector to be in the right file offset
  double * real = read_vec(f, n);
  delete[] real;

  //read the precition of size n (the main diagonal of the matrix A)
  double * prec = read_vec(f, n);
  assert(n == graph.num_vertices());
  for (int i=0; i< (int)n; i++){
     vertex_data_inv &data= graph.vertex_data(i);
     data.prior_prec = prec[i];
  }
  delete[] prec;

  e = read_edges<edata>(f, sizeof(edge_data)/sizeof(sdouble), 0, n, &graph, config);
  fclose(f);
}


/**
 * READ A NON-SQUARE matrix of size m rows x n cols 
 * Where the observation y is a vector of size m, the solution vector x=A\y is 
 * a vector of size n.
 */
template<typename graph_type, typename vertex_data, typename edge_data>
void load_non_square_matrix(FILE * f, graph_type& graph, advanced_config & config) {
  
  assert( n > 0);
  assert( m > 0);
  assert(m!=n); 
  
  printf("Loading a non-square matrix A of size %d x %d\n", m,n);

  if (config.supportgraphlabcf){ //read matrix factorization file (GraphLab collabrative filtering format)
     int tmp;
     fread(&tmp, 1, 4, f); //skip over time bin number 
     vertex_data data;
     for (int i=0; i< (int)(m+n); i++){
       graph.add_vertex(data);
     }
     if (config.isfloat)
          e = read_edges<edata2>(f, sizeof(edge_data), 0, n+m, &graph, config);
     else  
          e = read_edges<edata3>(f, sizeof(edge_data), 0, n+m, &graph, config);
  }
  else { //read A, x, y, from file
  
    //read y (the observation) of size m
    read_nodes(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_PRIOR_MEAN_OFFSET,m,&graph);
    
    //read x (the real solution, if given) of size n
    read_nodes(f, sizeof(vertex_data)/sizeof(sdouble),
             GABP_REAL_OFFSET,n,&graph);
    //read the precision vector of size m+n (of both solution and observation)
    double * prec = read_vec(f, n+m);
    dispatch_vec(0,n+m,GABP_PRIOR_PREC_OFFSET, &graph, prec, n+m, true);
    dispatch_vec(0,n+m,GABP_PREV_MEAN_OFFSET, &graph, 1);
    dispatch_vec(0,n+m,GABP_PREV_PREC_OFFSET, &graph, 1);
    e = read_edges<edata>(f, sizeof(edge_data), 0, n+m, &graph,config);
  }
  fclose(f);
}

template <typename coretype>
double start_inv(graphlab::command_line_options &clopts, advanced_config &config){

  assert(config.algorithm == GaBP_INV);
  // Create a core
  coretype core;
  core.set_engine_options(clopts); // Set the engine options

  // Load the graph --------------------------------------------------
  FILE * f = load_matrix_metadata(config.datafile.c_str());
  if (m == n){ //square matrix
     config.square = true;
     load_square_matrix(f, core.graph(), config);
  }
  else {
     assert(false);
  }


  // CREATE INITIAL TASKS ******
  // TOOD is this the correct starting priority?
  double initial_priority = 1.0;
  core.add_task_to_all(gabp_update_inv_function, initial_priority); 
  
  // Create an atomic entry to track iterations (as necessary)
  ITERATION_KEY.set(0);
  // Set all cosntants
  THRESHOLD_KEY.set(config.threshold);
  SUPPORT_NULL_VARIANCE_KEY.set(config.support_null_variance);
  ROUND_ROBIN_KEY.set(config.round_robin);
  DEBUG_KEY.set(config.debug);
  MAX_ITER_KEY.set(config.iter);


  //create a vector for storing the output
  std::vector<double> means(n);


  // START GRAPHLAB *****
  double runtime;
  double diff = 0;

  runtime= core.start();
  // POST-PROCESSING *****
  std::cout << algorithmnames[config.algorithm] << " finished in " << runtime << std::endl;
#ifdef HAS_ITPP

  itpp::mat ret = itpp::zeros(core.graph().num_vertices(), core.graph().num_vertices());
  for (size_t i = 0; i < core.graph().num_vertices(); i++){
    const vertex_data_inv& vdata = core.graph().vertex_data(i);
    
    for (int j=0; j< (int)core.graph().num_vertices(); j++){
      ret.set(i,j,vdata.cur_mean[j]);
    }   

   //diff += pow(vdata.real - vdata.cur_mean,2);
       //x[i] = vdata.cur_mean;
       //prec[i] = vdata.cur_prec;
       //TODO
  }

    if (config.debug){
      std::cout<<ret<<std::endl;
      if (config.unittest == 6){
         itpp::mat A = ("1.7011078828 0.4972882085 1.01358835; 0.4972882085 2.0077549 1.09088855; 1.01358835 1.09088855 2.690904500");
         std::cout<<A*ret<<std::endl;
      }
         
   }
#endif
  //TODO, measure qulity of solutio
  return diff;
}



template <typename coretype>
double start(graphlab::command_line_options &clopts, advanced_config &config){


  assert(config.algorithm != GaBP_INV);
  // Create a core
  coretype core;
  core.set_engine_options(clopts); // Set the engine options


  // Load the graph --------------------------------------------------
  FILE * f = load_matrix_metadata(config.datafile.c_str());
  if (m == n){ //square matrix
     config.square = true;
     load_square_matrix(f, core.graph(), config);
  }
  else {
     config.square = false;
     load_non_square_matrix<graph_type, vertex_data, edge_data>(f, core.graph(), config);
     if (config.algorithm == JACOBI){
        logstream(LOG_ERROR)<<" Jacobi can not run with non-square mastrix. Run with --sqaure=true and provide a square mastrix in the input file!" << std::endl;
        return EXIT_FAILURE;
     }
  }



  // CREATE INITIAL TASKS ******
  // TOOD is this the correct starting priority?
  double initial_priority = 1.0;
  switch(config.algorithm){
    case GaBP: 
	  core.add_task_to_all(gabp_update_function, initial_priority); break;

    case JACOBI:
          core.add_task_to_all(jacobi_update_function, initial_priority); break;

    case CONJUGATE_GRADIENT:
          //deliberately empty, will be done later
          break;
   
    default:
         logstream(LOG_ERROR) << "Unknown algorithm" << std::endl;
         clopts.print_description(); 
         return EXIT_FAILURE;
  } 

 // Initialize the shared data --------------------------------------
  // Set syncs
  //
  switch(config.algorithm){
     case JACOBI:
     case GaBP:
      if (config.syncinterval > 0){
       core.set_sync(REAL_NORM_KEY,
                gl_types::glshared_sync_ops::sum<double, get_real_norm>,
                apply_func_real,
                double(0),  config.syncinterval,
                gl_types::glshared_merge_ops::sum<double>);

        core.set_sync(RELATIVE_NORM_KEY,
                gl_types::glshared_sync_ops::sum<double, get_relative_norm>,
                apply_func_relative,
                double(0),  config.syncinterval,
                gl_types::glshared_merge_ops::sum<double>);

  	core.engine().add_terminator(termination_condition);
     }
  }
  // Create an atomic entry to track iterations (as necessary)
  ITERATION_KEY.set(0);
  // Set all cosntants
  THRESHOLD_KEY.set(config.threshold);
  SUPPORT_NULL_VARIANCE_KEY.set(config.support_null_variance);
  ROUND_ROBIN_KEY.set(config.round_robin);
  DEBUG_KEY.set(config.debug);
  MAX_ITER_KEY.set(config.iter);

  //core.graph().compute_coloring();

  //create a vector for storing the output
  std::vector<double> means(n);


  // START GRAPHLAB *****
  double runtime;
  double diff = 0;

  switch(config.algorithm){
      case GaBP:
      case JACOBI:
        runtime= core.start();
        break;

      case CONJUGATE_GRADIENT:
        runtime = cg(&core,means,diff,config);
        break;
  }
  // POST-PROCESSING *****
  std::cout << algorithmnames[config.algorithm] << " finished in " << runtime << std::endl;

  std::vector<double> precs(n);
  for (size_t i = m; i < core.graph().num_vertices(); i++){
    const vertex_data& vdata = core.graph().vertex_data(i);
     if (config.algorithm == JACOBI || config.algorithm == GaBP){
       diff += pow(vdata.real - vdata.cur_mean,2);
       means[i-m] = vdata.cur_mean;
       precs[i-m] = vdata.cur_prec;
     }
     //TODO: else if (algorithm == CONJUGATE_GRADIENT)
     //   diff += pow(means[i-m] - vdata.real,2);
  }

  if (config.algorithm == JACOBI || config.algorithm == GaBP){
  std::cout << "Assuming the linear system is Ax=y, and the correct solution is x*," << algorithmnames[config.algorithm] << " converged to an accuracy norm(x-x*) of " << diff
            << " msg norm is: " << RELATIVE_NORM_KEY.get_val() << std::endl;
   }

   f = fopen((config.datafile+".out").c_str(), "w");
   assert(f!= NULL);

   std::cout<<"Writing result to file: "<<config.datafile<<".out"<<std::endl;
   std::cout<<"You can read the file in Matlab using the load_c_gl.m matlab script"<<std::endl;
   write_vec(f, means.size(), &means[0]);
   if (config.algorithm == GaBP)
     write_vec(f, means.size(), &precs[0]);

   fclose(f);

   return diff;
}



int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

#ifdef ITPP
  logstream(LOG_WARNING) << "it++ detected. adding support of conjugate gradient !" << std::endl;
#endif
  logstream(LOG_INFO) << "GraphLab Linear solver library code by Danny Bickson, CMU" << std::endl <<
                         "Send comments and bug reports to danny.bickson@gmail.com" << std::endl <<
                         "Currently implemented algorithms are: Gaussian Belief Propagation, Jacobi method, Conjugate Gradient" << std::endl;

  // Setup additional command line arguments for the GABP program
  config.threshold = 1e-5;
  config.support_null_variance = false;
  config.iter = 0;
  config.syncinterval = 10000;
  config.unittest = 0;

  clopts.attach_option("algorithm", &config.algorithm, "Algorithm 0=Gaussian BP, 1= Jacobi");
  clopts.add_positional("algorithm");
 
  clopts.attach_option("data", &config.datafile, "Binary input file (as created by the save_c_gl.m script)");
  clopts.add_positional("data");
  clopts.attach_option("threshold", &config.threshold, config.threshold, "termination threshold.");
  clopts.add_positional("threshold");
  clopts.attach_option("nullvar", &config.support_null_variance, config.support_null_variance,
                       "(optional) support invalid covariance matrices - with null variance.");
  clopts.attach_option("square", &config.square, config.square, "is the matrix square? ");
  clopts.attach_option("debug", &config.debug, config.debug, "Display debug output.");
  clopts.attach_option("syncinterval", &config.syncinterval, config.syncinterval, "sync interval (number of update functions before convergen detection");
  clopts.attach_option("supportgraphlabcf", &config.supportgraphlabcf, config.supportgraphlabcf, "input is given in GraphLab collaborative filtering format");
  clopts.attach_option("isfloat", &config.isfloat, config.isfloat, "input file is given in float format");
  clopts.attach_option("cg_resid", &config.cg_resid, config.cg_resid, "compute cg residual progress ");
  clopts.attach_option("zero", &config.zero, config.zero, "support sparse matrix entry containing zero val ");
  clopts.attach_option("unittest", &config.unittest, config.unittest, "unit testing ( allowed values: 1/2)");
 
  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  if (config.unittest == 1){
    logstream(LOG_WARNING)<< "Going to run GaBP unit testing using matrix of size 3x2" << std::endl;
    const char * args[] = {"gabp", "0", "mat3x2", "1e-10", "--unittest=1","--syncinterval=10000"};
    clopts.parse(6, (char**)args);
    clopts.set_scheduler_type("round_robin(max_iterations=10000,block_size=1)");
  }
  else if (config.unittest == 2){
    logstream(LOG_WARNING)<< "Going to run GaBP unit testing using matrix of size 3x3" << std::endl;
    const char * args[] = {"gabp", "0", "mat3x3", "1e-10", "--unittest=2", "--syncinterval=10000"};
    clopts.parse(6, (char**)args);
    clopts.set_scheduler_type("round_robin(max_iterations=10000,block_size=1)");
  }
  else if (config.unittest == 3){
    logstream(LOG_WARNING)<< "Going to run Jacobi unit testing using matrix of size 3x3" << std::endl;
    const char * args[] = {"gabp", "1", "mat3x3", "1e-10", "--unittest=3", "--syncinterval=10000"};
    clopts.parse(6, (char**)args);
    clopts.set_scheduler_type("round_robin(max_iterations=10000,block_size=1)");
  }
  else if (config.unittest == 4){
    logstream(LOG_WARNING)<< "Going to run CG unit testing using matrix of size 3x3" << std::endl;
    const char * args[] = {"gabp", "2", "mat3x3", "1e-10", "--unittest=4", "--syncinterval=10000"};
    clopts.parse(6, (char**)args);
    config.iter = 10;
    clopts.set_scheduler_type("fifo");
  }
 else if (config.unittest == 5){
    logstream(LOG_WARNING)<< "Going to run CG unit testing using matrix of size 3x2" << std::endl;
    const char * args[] = {"gabp", "2", "mat3x2", "1e-10", "--unittest=5", "--syncinterval=10000"};
    clopts.parse(6, (char**)args);
    config.iter = 10;
    clopts.set_scheduler_type("fifo");
  }
 else if (config.unittest == 6){
    logstream(LOG_WARNING)<< "Going to run GaBP inverse unit testing using matrix of size 3x3" << std::endl;
    const char * args[] = {"gabp", "3", "mat3x3", "1e-10", "--unittest=6", "--syncinterval=10000"};
    clopts.parse(6, (char**)args);
    config.iter = 10;
    MATRIX_WIDTH_KEY.set(3);
    clopts.set_scheduler_type("round_robin(max_iterations=20,block_size=1)");
 }

  // Ensure that a data file is provided
  if(!clopts.is_set("data")) {
    logstream(LOG_ERROR)<<"No data file provided!" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }
 // Ensure that an algorithm is provided
  if(!clopts.is_set("algorithm")) {
    logstream(LOG_ERROR)<< "No algorithm provided! Choose from: 0) GaBP 1) Jacobi 2) CG" << std::endl;
    clopts.print_description();
    return EXIT_FAILURE;
  }

  double diff = 0;

  switch(config.algorithm){
    case GaBP:
    case JACOBI:
    case CONJUGATE_GRADIENT:
       diff=start<gl_types::core>(clopts,config);
       break;
    case GaBP_INV:
       diff=start_inv<gl_types_inv::core>(clopts,config);
       break;
  }
 
  //print timing counters
  for (int i=0; i<MAX_COUNTER; i++){
    if (counter[i] > 0)
    	printf("Performance counters are: %d) %s, %g\n",i, countername[i], counter[i]); 
   }


   if (config.unittest == 1){
      assert(diff <= 1.7e-4);
   }
   else if (config.unittest == 2){
      assert(diff <= 1e-15);
   }
   else if (config.unittest == 3){
      assert(diff <= 1e-30);
   }
   else if (config.unittest == 4 || config.unittest == 5)
      assert(diff < 1e-14);
   
   return EXIT_SUCCESS;
}

#include <graphlab/macros_undef.hpp>

