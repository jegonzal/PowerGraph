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


#include "graphlab.hpp"
#include "../shared/io.hpp"
#include "../shared/types.hpp"
#include "../shared/mathlayer.hpp"

//#define USE_GRAPH2
#ifdef USE_GRAPH2
#include "graphlab/graph/graph2.hpp"
#else
#include "graphlab/graph/graph3.hpp"
#endif
using namespace graphlab;
using namespace std;

DECLARE_TRACER(matproduct)

/**
 *
 *  Implementation of the Lanczos algorithm, as given in:
 *  http://en.wikipedia.org/wiki/Lanczos_algorithm
 * 
 *  Code written by Danny Bickson, CMU, June 2011
 * */



//LANCZOS VARIABLES
int max_iter = 10;
bool debug;
int nv = 0;
int nsv = 0;
double tol = 1e-8;
bool finished = false;
double regularization = 0;
double ortho_repeats = 3;

struct vertex_data {
  vec pvec;
  double value;
  double A_ii;
  vertex_data(){}
  void add_self_edge(double value) { A_ii = value; }

  void set_val(double value, int field_type) { 
    pvec[field_type] = value;
  }
  //double get_output(int field_type){ return pred_x; }
}; // end of vertex_data

struct edge_data {
  real_type weight;
  edge_data(double weight = 0) : weight(weight) { }
};

int data_size = max_iter;
#ifdef USE_GRAPH2
typedef graphlab::graph2<vertex_data, edge_data> graph_type;
#else
typedef graphlab::graph3<vertex_data, edge_data> graph_type;
#endif
#include "../shared/math.hpp"
#include "../shared/printouts.hpp"


void init_lanczos(graph_type * g, bipartite_graph_descriptor & info){

  data_size = nsv + nv+1 + max_iter;
  for (int i=0; i< info.total(); i++)
      g->vertex_data(i).pvec = zeros(info.is_square() ? (2*data_size): data_size);
}



vec lanczos(graphlab::core<graph_type, Axb> & glcore, bipartite_graph_descriptor & info, timer & mytimer, vec & errest, 
            const std::string & vecfile){
   

   int nconv = 0;
   int its = 1;
   int mpd = 24;
   int N = std::min(info.rows, info.cols);
   DistMat A(info);
   DistSlicedMat U(info.is_square() ? data_size : 0, info.is_square() ? 2*data_size : data_size, true, info, "U");
   DistSlicedMat V(0, data_size, false, info, "V");
   vec alpha, beta, b;
   vec sigma = zeros(nv);
   errest = zeros(nv);
   DistVec v_0(info, 0, false, "v_0");
   if (vecfile.size() == 0)
     v_0 = randu(size(A,2));
   PRINT_VEC2("svd->V", v_0);
/* Example Usage:
  DECLARE_TRACER(classname_someevent)
  INITIALIZE_TRACER(classname_someevent, "hello world");
  Then later on...
  BEGIN_TRACEPOINT(classname_someevent)
  ...
  END_TRACEPOINT(classname_someevent)
 */
   DistDouble vnorm = norm(v_0);
   v_0=v_0/vnorm;
   PRINT_INT(nv);

   while(nconv < nv && its < max_iter){
     logstream(LOG_INFO)<<"Starting iteration: " << its << " at time: " << mytimer.current_time() << std::endl;
     int k = nconv;
     int n = nv;
     PRINT_INT(k);
     PRINT_INT(n);

     alpha = zeros(n);
     beta = zeros(n);

     U[k] = V[k]*A._transpose();
     orthogonalize_vs_all(U, k, alpha(0));
     //alpha(0)=norm(U[k]).toDouble(); 
     PRINT_VEC3("alpha", alpha, 0);
     //U[k] = U[k]/alpha(0);

     for (int i=k+1; i<n; i++){
       PRINT_INT(i);

       V[i]=U[i-1]*A;
       orthogonalize_vs_all(V, i, beta(i-k-1));
      
       //beta(i-k-1)=norm(V[i]).toDouble();
       //V[i] = V[i]/beta(i-k-1);
       PRINT_VEC3("beta", beta, i-k-1); 
      
       U[i] = V[i]*A._transpose();
       orthogonalize_vs_all(U, i, alpha(i-k));
       //alpha(i-k)=norm(U[i]).toDouble();

       //U[i] = U[i]/alpha(i-k);
       PRINT_VEC3("alpha", alpha, i-k);
     }

     V[n]= U[n-1]*A;
     orthogonalize_vs_all(V, n, beta(n-k-1));
     //beta(n-k-1)=norm(V[n]).toDouble();
     PRINT_VEC3("beta", beta, n-k-1);

  //compute svd of bidiagonal matrix
  PRINT_INT(nv);
  PRINT_NAMED_INT("svd->nconv", nconv);
  PRINT_NAMED_INT("svd->mpd", mpd);
  n = nv - nconv;
  PRINT_INT(n);
  alpha.conservativeResize(n);
  beta.conservativeResize(n);

  PRINT_MAT2("Q",eye(n));
  PRINT_MAT2("PT",eye(n));
  PRINT_VEC2("alpha",alpha);
  PRINT_VEC2("beta",beta);
 
  mat T=diag(alpha);
  for (int i=0; i<n-1; i++)
    set_val(T, i, i+1, beta(i));
  PRINT_MAT2("T", T);
  mat a,PT;
  svd(T, a, PT, b);
  PRINT_MAT2("Q", a);
  alpha=b.transpose();
  PRINT_MAT2("alpha", alpha);
  for (int t=0; t< n-1; t++)
     beta(t) = 0;
  PRINT_VEC2("beta",beta);
  PRINT_MAT2("PT", PT.transpose());

  //estiamte the error
  int kk = 0;
  for (int i=nconv; i < nv; i++){
    int j = i-nconv;
    PRINT_INT(j);
    sigma(i) = alpha(j);
    PRINT_NAMED_DBL("svd->sigma[i]", sigma(i));
    PRINT_NAMED_DBL("Q[j*n+n-1]",a(n-1,j));
    PRINT_NAMED_DBL("beta[n-1]",beta(n-1));
    errest(i) = abs(a(n-1,j)*beta(n-1));
    PRINT_NAMED_DBL("svd->errest[i]", errest(i));
    if (alpha(j) >  tol){
      errest(i) = errest(i) / alpha(j);
      PRINT_NAMED_DBL("svd->errest[i]", errest(i));
    }
    if (errest(i) < tol){
      kk = kk+1;
      PRINT_NAMED_INT("k",kk);
    }


    if (nconv +kk >= nsv){
      printf("set status to tol\n");
      finished = true;
    }
  }//end for
  PRINT_NAMED_INT("k",kk);


  vec v;
  if (!finished){
    vec swork=get_col(PT,kk); 
    PRINT_MAT2("swork", swork);
    v = zeros(size(A,1));
    for (int ttt=nconv; ttt < nconv+n; ttt++){
      v = v+swork(ttt-nconv)*(V[ttt].to_vec());
    }
    PRINT_VEC2("svd->V",V[nconv]);
    PRINT_VEC2("v[0]",v); 
  }


INITIALIZE_TRACER(matproduct, "computing ritz eigenvectors");
   //compute the ritz eigenvectors of the converged singular triplets
  if (kk > 0){
    PRINT_VEC2("svd->V", V[nconv]);
BEGIN_TRACEPOINT(matproduct);
    mat tmp= V.get_cols(nconv,nconv+n)*PT;
    V.set_cols(nconv, nconv+kk, get_cols(tmp, 0, kk));
    PRINT_VEC2("svd->V", V[nconv]);
    PRINT_VEC2("svd->U", U[nconv]);
    tmp= U.get_cols(nconv, nconv+n)*a;
END_TRACEPOINT(matproduct);
    U.set_cols(nconv, nconv+kk,get_cols(tmp,0,kk));
    PRINT_VEC2("svd->U", U[nconv]);
  }

  nconv=nconv+kk;
  if (finished)
    break;

  V[nconv]=v;
  PRINT_VEC2("svd->V", V[nconv]);
  PRINT_NAMED_INT("svd->nconv", nconv);

  its++;
  PRINT_NAMED_INT("svd->its", its);
  PRINT_NAMED_INT("svd->nconv", nconv);
  //nv = min(nconv+mpd, N);
  if (nsv < 10)
    nv = 10;
  PRINT_NAMED_INT("nv",nv);

} // end(while)

printf(" Number of computed signular values %d",nconv);
printf("\n");
  DistVec normret(info, nconv, false, "normret");
  DistVec normret_tranpose(info, nconv, true, "normret_tranpose");
  for (int i=0; i < nconv; i++){
    normret = V[i]*A._transpose() -U[i]*sigma(i);
    double n1 = norm(normret).toDouble();
    PRINT_DBL(n1);
    normret_tranpose = U[i]*A -V[i]*sigma(i);
    double n2 = norm(normret_tranpose).toDouble();
    PRINT_DBL(n2);
    double err=sqrt(n1*n1+n2*n2);
    PRINT_DBL(err);
    PRINT_DBL(tol);
    if (sigma(i)>tol){
      err = err/sigma(i);
    }
    PRINT_DBL(err);
    PRINT_DBL(sigma(i));
    printf("Singular value %d \t%13.6g\tError estimate: %13.6g\n", i, sigma(i),err);
  }

  return sigma;
}

int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string datafile, vecfile;
  std::string format = "matrixmarket";
  int unittest = 0;

  clopts.attach_option("data", &datafile, datafile,
                       "matrix input file");
  clopts.add_positional("data");
  clopts.attach_option("initial_vector", &vecfile, vecfile,"optional initial vector");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("max_iter", &max_iter, max_iter, "max iterations");
  clopts.attach_option("ortho_repeats", &ortho_repeats, ortho_repeats, "orthogonalization iterations. 1 = low accuracy but fast, 2 = medium accuracy, 3 = high accuracy but slow.");
  clopts.attach_option("nv", &nv, nv, "Number of vectors in each iteration");
  clopts.attach_option("nsv", &nsv, nsv, "Number of requested singular values to comptue"); 
  clopts.attach_option("regularization", &regularization, regularization, "regularization");
  clopts.attach_option("tol", &tol, tol, "convergence threshold");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab V2 matrix factorization library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Lanczos" << std::endl;


  // Create a core
  graphlab::core<graph_type, Axb> core;
  omp_set_num_threads(clopts.get_ncpus());
  core.set_options(clopts); // Set the engine options
  

#ifndef UPDATE_FUNC_IMPL
   Axb mathops;
   core.add_aggregator("Axb", mathops, 1000);
#endif

  //unit testing
  if (unittest == 1){
    datafile = "gklanczos_testA"; 
    vecfile = "gklanczos_testA_v0";
    nsv = 3; nv = 3;
    debug = true;
    core.set_scheduler_type("sweep(ordering=ascending,strict=true)");
    core.set_ncpus(1);
  }
  else if (unittest == 2){
    datafile = "gklanczos_testB";
    vecfile = "gklanczos_testB_v0";
    nsv = 10; nv = 10;
    debug = true;  max_iter = 100;
    core.set_scheduler_type("sweep(ordering=ascending,strict=true)");
    core.set_ncpus(1);
  }
  else if (unittest == 3){
    datafile = "gklanczos_testC";
    vecfile = "gklanczos_testC_v0";
    nsv = 4; nv = 10;
    debug = true;  max_iter = 100;
    core.set_scheduler_type("sweep(ordering=ascending,strict=true)");
    core.set_ncpus(1);
  }


  std::cout << "Load matrix " << datafile << std::endl;
#ifdef USE_GRAPH2
  load_graph(datafile, format, info, core.graph(), MATRIX_MARKET_3, false, false);
  core.graph().finalize();
#else  
  load_graph(datafile, format, info, core.graph(), MATRIX_MARKET_3, false, true);
  core.graph().load_directed(datafile, false, false);
#endif
  init_lanczos(&core.graph(), info);
  init_math(&core.graph(), &core, info, ortho_repeats, false);
  if (vecfile.size() > 0){
    std::cout << "Load inital vector from file" << vecfile << std::endl;
    load_vector(vecfile, format, info, core.graph(), 0, true, false);
  }  
 
  timer mytimer; mytimer.start(); 
  vec errest;
  vec eigenvalues = lanczos(core, info, mytimer, errest, vecfile);
 
  std::cout << "Lanczos finished in " << mytimer.current_time() << std::endl;
  std::cout << "\t Updates: " << core.last_update_count() << " per node: " 
     << core.last_update_count() / core.graph().num_vertices() << std::endl;

  //vec ret = fill_output(&core.graph(), bipartite_graph_descriptor, JACOBI_X);

  //write_output_vector(datafile + "x.out", format, ret);


  if (unittest == 1){
    assert(errest.size() == 3);
    for (int i=0; i< errest.size(); i++)
      assert(errest[i] < 1e-30);
  }
  else if (unittest == 2){
     assert(errest.size() == 10);
    for (int i=0; i< errest.size(); i++)
      assert(errest[i] < 1e-15);
  }


   return EXIT_SUCCESS;
}


