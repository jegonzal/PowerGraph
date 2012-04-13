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

//data file types
enum testtype{
  TRAINING = 0,
  VALIDATION = 1,
  TEST = 2
};
const char * testtypename[] = {"TRAINING", "VALIDATION", "TEST"};

struct problem_setup{
  int L; //number of edges in training file
  int Le; //number of edges in validation file
  int Lt; //number of edges in test file
  int its; //current iteration number
  timer gt; //timer of tracking program run time
  int M; //number of rows
  int N; //number of columns
  double pU; //regularization for users
  double pV; //regularization for movies

  problem_setup(){
    L = Le = Lt = 0;
    its = 0;
    gt.start();
  }
};

bool no_edge_data; //ugly, to be fixed later

struct advanced_config{

  int max_iter; //maximal number of iterations
  bool printhighprecision; //print RMSE using more decimal digits
  int D; //feature vector length
  mat eDT; //eye() matrix used for regularization
  //bool no_edge_data;
  double tol;
  bool finished;
  double maxval; //maximum allowed value
  double minval; //minimum allowed value
  bool update_function;
  bool save_vectors;
  std::string datafile; //input file name
  std::string format; //input file format
  int nodes;
  bool zero; //support zero ratings

  advanced_config(){ 
		max_iter = 10;
		debug = false;
		printhighprecision = false;
		D = 20; //feature vector length
		no_edge_data = false;
		tol = 1e-8;
		finished = false;
		maxval = 1e100;
		minval = -1e100;
		update_function = false;
		save_vectors = false;
		format = "matrixmarket";
		nodes = 0;
    zero = false;
  }
};

//

/**
 *
 *  Implementation of the alternating least squares via LAPACK / eigen
 * 
 *  Code written by Danny Bickson, CMU, June 2011
 * */


problem_setup ps;
advanced_config ac;


struct vertex_data {
  vec pvec;
  double value;
  double A_ii;
  uint num_edges;

  vertex_data(){ num_edges = 0; value = 0; A_ii = 0; }
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

#ifdef USE_GRAPH2
typedef graphlab::graph2<vertex_data, edge_data> graph_type;
#else
typedef graphlab::graph3<vertex_data, edge_data> graph_type;
#endif

int data_size;

#include "../shared/math.hpp"
#include "../shared/printouts.hpp"
#include "../shared/least_squares.hpp"

graph_type validation, test;
graph_type * training = NULL;


double calc_obj(double res){

  double sumU = 0, sumV = 0, sumT = 0;
  int user_cnt = 0;
  int movie_cnt = 0;
  int edges = 0;
  const graph_type * g = training;

#pragma omp parallel for
  for (int i=0; i< ps.M; i++){
    const vertex_data * data = &g->vertex_data(i);
    if (data->num_edges > 0){
      sumU += sum_sqr(data->pvec);
      user_cnt++;
      edges+= data->num_edges;
    }

  }
#pragma omp parallel for
  for (int i=ps.M; i< ps.M+ps.N; i++){
    const vertex_data * data = &g->vertex_data(i);
    if (data->num_edges > 0 ){
      sumV += sum_sqr(data->pvec);
      movie_cnt++;
      edges+= data->num_edges;
    }
  }
                                  
  double obj = (res +ps.pU*sumU + ps.pV*sumV + sumT) / 2.0;

  if (debug)
     cout<<"OBJECTIVE: res: " << res << "sumU " << sumU << " sumV: " << sumV << " pu " << ps.pU << " pV: " << ps.pV << endl;

  //assert(edges == 2*ps.L);
  return obj;
}



float calc_rmse_edge(const edge_data & edge, const graph_type *_g, double & rmse, const vertex_data&data, const vertex_data&pdata, int& e, int i){
   float prediction;
   e++;
   return predict(data, pdata, edge.weight, prediction);
} 

double calc_rmse(const graph_type * _g, bool validation, double & res, const bipartite_graph_descriptor & info){
     if (validation && ps.Le == 0)
       return NAN;
     
     res = 0;
     double RMSE = 0;
     int e = 0;
     uint i;
//#pragma omp parallel for private(i) reduction(+: RMSE)
//#pragma omp parallel for
     for (i=info.get_start_node(false); i< (uint)info.get_end_node(false); i++){
       const vertex_data & data = training->vertex_data(i);
       for (uint j=0; j< _g->num_in_edges(i); j++) {
         const edge_type & edget = _g->in_edges(i)[j];
         assert(edget.source() != i);
         const vertex_data & pdata = training->vertex_data(edget.source()); 
         RMSE = RMSE + calc_rmse_edge(_g->edge_data(edget), _g, RMSE, data, pdata, e, i);       
     }
   }
   res = RMSE;
   if (e != (validation?ps.Le:ps.L))
      logstream(LOG_FATAL)<<"Missing ratings in " << testtypename[validation] << " file. Expected to have "
      << (validation?ps.Le:ps.L) << " while encountered only " << e << std::endl;
   return sqrt(RMSE/(double)e);

}

double agg_rmse(){
  double res = 0;
  for (int i=info.get_start_node(false); i< info.get_end_node(false); i++){
     res += training->vertex_data(i).value;
  }
  return res;
}


void last_iter(){
  printf("Entering last iter with %d\n", ps.its ); 
  double res,res2;
  //double rmse = calc_rmse(training, false, res, info); //TODO (ps.algorithm != STOCHASTIC_GRADIENT_DESCENT && ps.algorithm != NMF) ? agg_rmse_by_movie<graph_type,vertex_data>(res) : agg_rmse_by_user<graph_type,vertex_data>(res);
  //rmse=0;
  res = agg_rmse();
  printf(ac.printhighprecision ? 
  "%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.12f VALIDATION RMSE=%0.12f.\n":
  "%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n"
  , ps.gt.current_time(), "ALS", ps.its, calc_obj(res),  sqrt(res/ps.L), calc_rmse(&validation, true, res2, info));
 ps.its++;
}

 
void init_als(graph_type * g, const bipartite_graph_descriptor & info){

  if (g->num_vertices() == 0)
     logstream(LOG_FATAL)<<"Failed to load graph. Aborting" << std::endl;

#pragma omp parallel for
  for (int i=0; i< info.total(); i++){
      g->vertex_data(i).pvec = (debug ?  ones(ac.D)*0.1 : randu(ac.D)*0.1);
   }
   logstream(LOG_INFO)<<"Allocated a total of: " << 
     ((double)g->num_vertices() * ac.D * sizeof(double)/ 1e6) << " MB for storing factor matrices." << std::endl;

   mi.eDT = eye(ac.D);
   data_size = ac.D;
   mi.maxval = ac.maxval; 
   mi.minval = ac.minval;
  ps.pU = ps.pV = regularization;
   printf("pU=%g, pV=%g, D=%d\n", ps.pU, ps.pV, ac.D);
}


void als(graphlab::core<graph_type, als_lapack> & glcore,
	  const bipartite_graph_descriptor & info, vec & errest, 
            const std::string & vecfile){
   

   ps.its = 1;
   DistMat A(info);
   assert(data_size > 0);
   DistSlicedMat U(0, data_size, true, info, "U");
   DistSlicedMat V(0, data_size, false, info, "V");

   while(ps.its < ac.max_iter){
     logstream(LOG_INFO)<<"Starting ALS iteration: " << ps.its << " at time: " << ps.gt.current_time() << std::endl;
     V = A.backslash(U);
     U = A.backslash(V);
     last_iter();
   } // end(while)

}

void common_prediction(const graph_type &g, const graph_type & _g, const vertex_data& data,int i, int &lineNum, double& sumPreds, vec& test_prediction){
      edge_list edges = _g.in_edges(i);
      for (uint j=0; j< edges.size(); j++) {
          const vertex_data & pdata = g.vertex_data(edges[j].target()); 
	        const edge_data & edge = _g.edge_data(edges[j]);
          
          if (!ac.zero)
           	assert(edge.weight != 0);

          float prediction = 0;
          predict(data, pdata, edge.weight, prediction);
    
         if (debug && (i== 0 || i == ps.M))
            cout<<lineNum<<") prediction:"<<prediction<<endl; 
            
         test_prediction[lineNum] = prediction;
	       sumPreds += prediction;
 	       lineNum++; 
       }
}


int main(int argc,  char *argv[]) {
  
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);

  graphlab::command_line_options clopts("GraphLab Linear Solver Library");

  std::string vecfile;
  int unittest = 0;

  clopts.attach_option("data", &ac.datafile, ac.datafile,
                       "matrix input file");
  clopts.add_positional("data");
  //clopts.attach_option("initial_vector", &vecfile, vecfile,"optional initial vector");
  clopts.attach_option("debug", &debug, debug, "Display debug output.");
  clopts.attach_option("unittest", &unittest, unittest, 
		       "unit testing 0=None, 1=3x3 matrix");
  clopts.attach_option("max_iter", &ac.max_iter, ac.max_iter, "max iterations");
  clopts.attach_option("regularization", &regularization, regularization, "regularization");
  clopts.attach_option("update_function", &ac.update_function, ac.update_function, "true = use update function. false = user aggregator");
  clopts.attach_option("save_vectors", &ac.save_vectors, ac.save_vectors, "save output matrices U and V.");
  clopts.attach_option("nodes", &ac.nodes, ac.nodes, "number of rows/cols in square matrix (optional)");
  clopts.attach_option("no_edge_data", &no_edge_data, no_edge_data, "matrix is binary (optional)");
  clopts.attach_option("maxval", &ac.maxval, ac.maxval, "maximum allowed value in matrix");
  clopts.attach_option("minval", &ac.minval, ac.minval, "minimum allowed value in matrix");
  clopts.attach_option("regnormal", &regnormal, regnormal, "Weight regularization with the number of neighbors for each graph node");
  clopts.attach_option("D", &ac.D, ac.D, "Feature vector width. Larger D is more accurate but slower."); 
  clopts.attach_option("printhighprecision", &ac.printhighprecision, ac.printhighprecision, "print high precision RMSE values");

  // Parse the command line arguments
  if(!clopts.parse(argc, argv)) {
    std::cout << "Invalid arguments!" << std::endl;
    return EXIT_FAILURE;
  }


  if (ac.update_function)
    logstream(LOG_INFO)<<"Going to use update function mechanism. " << std::endl;

  logstream(LOG_WARNING)
    << "Eigen detected. (This is actually good news!)" << std::endl;
  logstream(LOG_INFO) 
    << "GraphLab V2 matrix factorization library code by Danny Bickson, CMU" 
    << std::endl 
    << "Send comments and bug reports to danny.bickson@gmail.com" 
    << std::endl 
    << "Currently implemented algorithms are: Lanczos, One Sides Lanzaos and Alternating least squares." << std::endl;


  // Create a core
  graphlab::core<graph_type, als_lapack> core;
  omp_set_num_threads(clopts.get_ncpus());
  core.set_options(clopts); // Set the engine options
  

#ifndef UPDATE_FUNC_IMPL
   als_lapack mathops;
   core.add_aggregator("als_lapack", mathops, 1000);
#endif

  //unit testing
  if (unittest == 1){
    ac.datafile = "ALSA"; 
    ac.D = 20; ac.max_iter = 50; regularization = 0.001;
    //debug = true;
    core.set_scheduler_type("sweep(ordering=ascending,strict=true)");
    core.set_ncpus(1);
  }
  else if (unittest == 2){
    ac.datafile = "gklanczos_testB";
    debug = true;  ac.max_iter = 100;
    core.set_scheduler_type("sweep(ordering=ascending,strict=true)");
    core.set_ncpus(1);
  }
  else if (unittest == 3){
    ac.datafile = "gklanczos_testC";
    debug = true;  ac.max_iter = 100;
    core.set_scheduler_type("sweep(ordering=ascending,strict=true)");
    core.set_ncpus(1);
  }

  std::cout << "Load matrix " << ac.datafile << std::endl;
  info.force_non_square = true;
#ifdef USE_GRAPH2
  load_graph(ac.datafile, ac.format, info, core.graph(), MATRIX_MARKET_3, false, false);
  core.graph().finalize();
#else
  load_graph(ac.datafile, ac.format, info, core.graph(), MATRIX_MARKET_3, false, true);
   core.graph().load_directed(ac.datafile, false, no_edge_data);
   info.nonzeros = core.graph().num_edges();
#endif
  ps.L = core.graph().num_edges();
  ps.M = info.rows;
  ps.N = info.cols;
  training = &core.graph();


#pragma omp parallel for
  for (int i=0; i< info.total(); i++){
    training->vertex_data(i).num_edges = training->num_in_edges(i) + training->num_out_edges(i);
     //cout << "node: " << i << " num edges: " << training->vertex_data(i).num_edges << std::endl;
   }

  bipartite_graph_descriptor infoe;
  infoe.force_non_square = true;
#ifdef USE_GRAPH2  
  bool found = load_graph(ac.datafile + "e", ac.format, infoe, validation, MATRIX_MARKET_3, false, false, true);
  validation.finalize();
#else
  bool found = load_graph(ac.datafile + "e", ac.format, infoe, validation, MATRIX_MARKET_3, false, true, true);
  if (found){
   validation.load_directed(ac.datafile + "e", false, no_edge_data);
   infoe.nonzeros = core.graph().num_edges();
  }
#endif
  ps.Le = validation.num_edges();
  if (ps.Le > 0 && (info.rows != infoe.rows || info.cols != infoe.cols))
    logstream(LOG_FATAL) << "Validation file has dimension " << infoe.rows << "x"<< infoe.cols << " while training file has dimension "
                         << info.rows << "x" << info.cols << ". Please fix your input to have the same dimensions. " << std::endl;
  bipartite_graph_descriptor infot;
  infot.force_non_square = true;
#ifdef USE_GRAPH2
  bool found = load_graph(ac.datafile + "t", ac.format, infot, test, MATRIX_MARKET_3, false, false, true);
  test.finalize();
#else
   found = load_graph(ac.datafile + "t", ac.format, infot, validation, MATRIX_MARKET_3, false, true, true);
   if (found){
     test.load_directed(ac.datafile + "t", false, no_edge_data);
     infot.nonzeros = core.graph().num_edges(); 
   }
#endif
  ps.Lt = test.num_edges();
  if (ps.Lt > 0 && (info.rows != infot.rows || info.cols != infot.cols))
    logstream(LOG_FATAL) << "Test file has dimension " << infot.rows << "x"<< infot.cols << " while training file has dimension "
                         << info.rows << "x" << info.cols << ". Please fix your input to have the same dimensions. " << std::endl;

  init_als(&core.graph(), info);
  init_math(&core.graph(), &core, info, ac.update_function);
  vec errest;
  als(core, info, errest, vecfile);
 
  std::cout << "Lanczos finished in " << ps.gt.current_time() << std::endl;
  std::cout << "\t Updates: " << core.last_update_count() << " per node: " 
     << core.last_update_count() / core.graph().num_vertices() << std::endl;


  if (ps.Lt > 0){
    double sumPreds = 0;
    int lineNum = 0;
    vec out_predictions = zeros(ps.Lt);

    for (int i=0; i< ps.M; i++){
      vertex_data & data = (vertex_data&)training->vertex_data(i);
      common_prediction(*training, test, data, i, lineNum, sumPreds, out_predictions);
    }

    assert(lineNum==ps.Lt); 
  logstream(LOG_INFO)<< "**Completed successfully (mean prediction: " << sumPreds/lineNum << std::endl;
    save_matrix_market_format_vector(ac.datafile + ".test.predictions",
     out_predictions, false, "output predictions for test data\n");
   }


  //write_output_vector(ac.datafile + ".singular_values", ac.format, singular_values,false, "%GraphLab SVD Solver library. This file contains the singular values.");

  if (unittest == 1){
  }
  else if (unittest == 2){
     assert(errest.size() == 10);
    for (int i=0; i< errest.size(); i++)
      assert(errest[i] < 1e-15);
  }


   return EXIT_SUCCESS;
}


