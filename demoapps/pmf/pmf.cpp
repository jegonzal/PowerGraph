// #define NDEBUG
#include <fstream>
#include <cmath>
#include <cstdio>

#include "graphlab.hpp"
#include "itppvecutils.hpp"
#include "pmf.h"
#include "prob.hpp"
#include "bptf.hpp"
//#include "svdpp.hpp"

#include <graphlab/macros_def.hpp>
/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU
 See documentation in header file pmf.h

*/

using namespace graphlab;
bool BPTF = true; //is this a sampling MCMC algo?
bool tensor = true; //is this tensor or a matrix
bool debug = false; //debug mode
bool ZERO = false; //support edges with zero weight
runmodes options; //type of algorithm
timer gt;
using namespace itpp;
using namespace std;
double scaling = 1; //aggregate time values into bins? (default =1, no aggregation)
double truncating =0; // truncate unused time bins (optional, default = 0, no truncation)
std::vector<edge_id_t> * edges;  //edge list for pointing from time nodes into their connected users/movies (for tensor)
std::string infile;
bool loadfactors = false; //start from a previous solution instead of a random point
bool savefactors = true; //save solution to file
bool loadgraph = false; //load graph from a binary saved file
bool savegraph = false; //save input file into graph binary file
bool stats = false; //print out statistics and exit
bool regnormal = false; //regular normalization
bool aggregatevalidation = false; //use validation dataset as training data
extern bool finish; //defined in convergence.hpp
int iiter = 1;//count number of time zero node run
double scalerating = 0; //scale the rating by dividing to the scalerating factor (optional)
int delayalpha = 0; //delay alpha sampling (optional, for BPTF)

/* Variables for PMF */
int M,N,K,L;//training size: users, movies, times, number of edges
int Le = 0; //number of ratings in validation dataset 
int Lt = 0;//number of rating in test data set
double pU = 10; //regularization for users
double pT = 1; //regularization for tensor time nodes
double pV = 10; //regularization for movies
double muT = 1; //mean of time nodes
vec vones; 
mat eDT; 
mat dp;

/* Variables for SVD++ */
float *svd_m_bias, * svd_c_bias, * svd_movie_weight, * svd_sum_user_weight;

double counter[20];

vertex_data * times = NULL;

gl_types::iengine * engine;
graph_type* g;
graph_type validation_graph;
graph_type test_graph;

#ifndef _min
#define _min(a,b) (a>b)?b:a
#endif


const size_t RMSE = 0;


 
/**
 * pmf APPLICATION CLASS
 */
    
/** CONSTRUCTOR **/
void init_pmf() {
  if (BPTF)
	   pT=10;
  eDT = itpp::eye(D)*pT;
  vones = itpp::ones(D);
}
    
void load_pmf_graph(const char* filename, graph_type * g, testtype flag,gl_types::core & glcore);    
void calc_T(int id);    
double calc_obj(double res);
void last_iter();
void export_kdd_format(graph_type * _g, bool dosave);
void calc_stats(testtype type);

// update function for time nodes
// this function is called only in tensor mode
void time_node_update_function(gl_types::iscope &scope, gl_types::icallback &scheduler) {

  assert(tensor);

  int id = scope.vertex() - M-N;
  assert(id >=0 && id < K);	
  if (debug && (id == 0 || id == K-1)){
    printf("entering time node  %d \n", id);   
  } 
  if (K > 1)
    calc_T(id); 
  else last_iter();
}


//calculate RMSE. This function is called only before and after grahplab is run.
//during run, calc_rmse_q is called 0 which is much lighter function (only aggregate sums of squares)
double calc_rmse(graph_type * _g, bool test, double & res){

     if (test && Le == 0)
       return NAN;
      
     res = 0;
     double RMSE = 0;
     int e = 0;
     for (int i=M; i< M+N; i++){
       vertex_data * data = &g->vertex_data(i);
       foreach(edge_id_t iedgeid, _g->in_edge_ids(i)) {
         vertex_data * pdata = &g->vertex_data(_g->source(iedgeid)); 
            
#ifndef GL_NO_MULT_EDGES
         multiple_edges & edges = _g->edge_data(iedgeid);
         for (int j=0; j< (int)edges.medges.size(); j++){       
           edge_data & edge = edges.medges[j];
#else
	   edge_data & edge = _g->edge_data(iedgeid);
#endif

           if (!ZERO)
           	assert(edge.weight != 0);

           double sum = 0; 
           double add = rmse(data->pvec, pdata->pvec, tensor? (&times[(int)edge.time].pvec):NULL, D, edge.weight, sum);
           if (!ZERO)
	      assert(sum != 0);         
           if (debug && (i== M || i == M+N-1) && (e == 0 || e == (test?Le:L)))
             cout<<"RMSE:"<<i <<"u1"<< data->pvec << " v1 "<< pdata->pvec<<endl; 

#ifndef GL_NO_MCMC
           if (BPTF && iiter > BURN_IN){
             edge.avgprd += sum;
             add = powf((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
           }
#endif
            
           RMSE+= add;
           e++;
#ifndef GL_NO_MULT_EDGES        
         }
#endif
     }
   }
   res = RMSE;
   assert(e == (test?Le:L));
   return sqrt(RMSE/(double)e);

}
   
// go over all edges and aggregate RMSE by summing up the squares, computing the
// mean and then sqrt()
double calc_rmse_q(double & res){

  timer t; t.start();
  res = 0;
  double RMSE = 0;
  for (int i=M; i< M+N; i++){ 
    const vertex_data * data = &g->vertex_data(i);
    RMSE+= data->rmse;
  }
  res = RMSE;
  counter[CALC_RMSE_Q] += t.current_time();
  return sqrt(RMSE/(double)L);

}

// fill out the linear relation matrix (Q) between users/movies and movie/users
inline void parse_edge(const edge_data& edge, const vertex_data & pdata, mat & Q, vec & vals, int i){
      
  if (!ZERO)  
  	assert(edge.weight != 0);

  if (tensor){
    dot2(pdata.pvec,  times[(int)edge.time].pvec, Q, i, D);  
  }
  else {
    for (int j=0; j<D; j++)
      Q.set(j,i, pdata.pvec[j]); 
  }
 
  vals[i] = edge.weight;
}
 

//count the number of edges connecting a user/movie to its neighbors
//(when there are multiple edges in different times we count the total)
int count_edges(edge_list es){
  
  if (options != BPTF_TENSOR_MULT && options != ALS_TENSOR_MULT)
      return es.size();

#ifndef GL_NO_MULT_EDGES
  int cnt = 0; 
  for (int j=0; j< (int)es.size(); j++){
    cnt += g->edge_data(es[j]).medges.size();
  }
  return cnt;
#else
  return es.size();
#endif
}


/***
 * UPDATE FUNCTION
 */
void user_movie_nodes_update_function(gl_types::iscope &scope, 
			 gl_types::icallback &scheduler) {
    

  /* GET current vertex data */
  vertex_data& vdata = scope.vertex_data();
 
  
  /* print statistics */
  if (debug&& (scope.vertex() == 0 || ((int)scope.vertex() == M-1) || ((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1) || ((int)scope.vertex() == 93712))){
    printf("entering %s node  %u \n", (((int)scope.vertex() < M) ? "movie":"user"), (int)scope.vertex());   
    debug_print_vec((((int)scope.vertex() < M) ? "V " : "U") , vdata.pvec, D);
  }

  assert((int)scope.vertex() < M+N);

  vdata.rmse = 0;

  int numedges = vdata.num_edges;
  if (numedges == 0){
    return; //if this user/movie have no ratings do nothing
  }


  edge_list outs = scope.out_edge_ids();
  edge_list ins = scope.in_edge_ids();
  timer t;
  mat Q(D,numedges); //linear relation matrix
  vec vals(numedges); //vector of ratings

  int i=0;

  t.start(); 
  //USER NODES    
  if ((int)scope.vertex() < M){

    foreach(graphlab::edge_id_t oedgeid, outs) {
#ifndef GL_NO_MULT_EDGES
      multiple_edges &medges =scope.edge_data(oedgeid);
#else
      edge_data & edge = scope.edge_data(oedgeid);
#endif
      const vertex_data  & pdata = scope.const_neighbor_vertex_data(scope.target(oedgeid)); 
#ifndef GL_NO_MULT_EDGES	   
      for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
#endif
        //go over each rating of a movie and put the movie vector into the matrix Q
        //and vector vals
        parse_edge(edge, pdata, Q, vals, i); 
        if (debug && ((int)scope.vertex() == 0 || (int)scope.vertex() == M-1) && (i==0 || i == numedges-1))
          std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;
        i++;
#ifndef GL_NO_MULT_EDGES
      }   
#endif
    }
  }

  else {


    //MOVIE NODES
    foreach(edge_id_t iedgeid, ins) {

      const vertex_data & pdata = scope.const_neighbor_vertex_data(scope.source(iedgeid)); 
#ifndef GL_NO_MULT_EDGES
       multiple_edges & medges =scope.edge_data(iedgeid);
        for (int j=0; j< (int)medges.medges.size(); j++){
        edge_data& edge = medges.medges[j];
#else
	edge_data & edge = scope.edge_data(iedgeid);
#endif   
        //go over each rating by user
        parse_edge(edge, pdata, Q, vals, i); 
        if (debug && (((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1)) && (i==0 || i == numedges-1))
          std::cout<<"set col: "<<i<<" " <<Q.get_col(i)<<" " <<std::endl;

        i++;
        double sum;     
        double trmse = rmse(vdata.pvec, pdata.pvec, tensor?(&times[(int)edge.time].pvec):NULL, D, edge.weight, sum);
        //assert(sum != 0);
#ifndef GL_NO_MCMC
        if (BPTF && iiter > BURN_IN){
          edge.avgprd += sum;        
          trmse = pow((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
        }
#endif
        vdata.rmse += trmse; 
 
#ifndef GL_NO_MULT_EDGES     
      }
#endif
      
    }

  }
  assert(i == numedges);
  counter[EDGE_TRAVERSAL] += t.current_time();

  vec result;
      
  if (!BPTF){
    //COMPUTE LEAST SQUARES (ALTERNATING LEAST SQUARES)
    t.start();
    double regularization = LAMBDA;
    //compute weighted regularization (see section 3.2 of Zhou paper)
    if (!regnormal)
	regularization*= Q.cols();

    bool ret = itpp::ls_solve(Q*itpp::transpose(Q)+eDT*regularization, Q*vals, result);
    assert(ret);
    counter[ALS_LEAST_SQUARES] += t.current_time();
  }
  else {
    //COMPUTE LEAST SQUARES (BPTF)
    //according to equation A.6 or A.7 in Xiong paper.
    assert(Q.rows() == D);
    t.start();
    mat iAi_;
    bool ret =inv(((int)scope.vertex() < M? A_U : A_V) + alpha *  Q*itpp::transpose(Q), iAi_);
    assert(ret);
    t.start();
    vec mui_ =  iAi_*((((int)scope.vertex() < M)? (A_U*mu_U) : (A_V*mu_V)) + alpha * Q * vals);
    counter[BPTF_LEAST_SQUARES2]+= t.current_time();
       
    t.start();
    result = mvnrndex(mui_, iAi_, D); 
    assert(result.size() == D);
    counter[BPTF_MVN_RNDEX] += t.current_time();
  }

  if (debug && (((int)scope.vertex()  == 0) || ((int)scope.vertex() == M-1) || ((int)scope.vertex() == M) || ((int)scope.vertex() == M+N-1))){
    std::cout <<(BPTF?"BPTF":"ALS")<<" Q: " << Q << std::endl<<" result: " << result << " edges: " << numedges << std::endl;
  }
      
  //store new result
  vdata.pvec =  result;

  //calc post round tasks
  if (!tensor && (int)scope.vertex() == M+N-1)
    last_iter();

}

void last_iter(){
  printf("Entering last iter with %d\n", iiter);

  double res,res2;
  double rmse = calc_rmse_q(res);
  //rmse=0;
  printf("%g) Iter %s %d  Obj=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", gt.current_time(), BPTF?"BPTF":"ALS", iiter,calc_obj(res),  rmse, calc_rmse(&validation_graph, true, res2));
  iiter++;
  if (iiter == BURN_IN && BPTF){
    printf("Finished burn-in period. starting to aggregate samples\n");
  }
         
  if (BPTF){
    timer t;
    t.start();
    if (iiter > delayalpha)
    	sample_alpha(res);
    sample_U();
    sample_V();
    if (tensor) 
      sample_T();
    counter[BPTF_SAMPLE_STEP] += t.current_time();
    if (infile == "kddcup" || infile == "kddcup2")
	export_kdd_format(&test_graph, !(iiter %15));
   }
}


// Calculate time nodes (for tensor)
void calc_T(int i){

  assert(tensor);
 
  assert(i >=0 && i < K);
  if (ZERO && edges[i].size() == 0)
	  return;

  assert(edges[i].size() > 0);

  if (debug && (i==0 || i == K-1))
    printf("node %d with Q size: %d %d\n", i, (int)edges[i].size(), (int)D);

  int k=0;
  int batchSize = 100;
  mat Q(batchSize, D); 
  vec vals(batchSize);
  mat QQ = zeros(D,D);
  vec RQ = zeros(D);
  int cnt =0;

  timer t; t.start();
  foreach (edge_id_t edge, edges[i]){
    if (k % batchSize == 0){
      Q.zeros();
      vals.zeros();
      cnt = 1;
    }

    //find the right edge which matches the current time
#ifndef GL_NO_MULT_EDGES    
     multiple_edges * medges= &g->edge_data(edge);
      edge_data data;
      for (int j=0; j< (int)medges->medges.size(); j++){
        data = medges->medges[j];
        if (data.time == i)
          break;
      }
#else
     edge_data data= g->edge_data(edge);
#endif
    assert(data.time == i);

    assert((int)g->target(edge)>= M);
    assert((int)g->source(edge)< M);
    vertex_data * v1 = &g->vertex_data(g->target(edge));
    vertex_data * v2 = &g->vertex_data(g->source(edge));
    vec ret = elem_mult(v1->pvec, v2->pvec);
     
    for (int s=0; s<D; s++)
      Q(k%batchSize,s)=ret(s);
    if (debug && (i==0 || i == K-1) && (k == 0 || k == (int)edges[i].size() - 1))
      std::cout<<" clmn "<<k<< " vec: " << ret<<std::endl;

    vals[k%batchSize] = data.weight;
    k++;

    if ((cnt  == batchSize) || (cnt < batchSize && k == (int)edges[i].size()-1)){
      QQ += transpose(Q)*Q;
      RQ += transpose(Q)*vals;
      assert(QQ.rows() == D && QQ.cols() == D);
    }
    cnt++;
  }
  counter[BPTF_TIME_EDGES] += t.current_time();


  if (debug && (i == 0 ||i== K-1 )){
    std::cout<<"QQ:"<<QQ<<std::endl;
    std::cout<<"RQ:"<<RQ<<std::endl;
  }

  assert(RQ.size() == D);
  assert((unsigned int)k == edges[i].size());
  vec sol;    
  vec out;

  t.start();
  if (i == 0){
    vec t1 = times[1].pvec; 
    //calc least squares estimation of time nodes
    if (!BPTF){
      QQ = QQ+ 2*eDT;
      bool ret = itpp::ls_solve(QQ, pT*(t1 + vones*muT)+ RQ, out);
      assert(ret);
    }
    //calc equation A.8 in Xiong paper
    else {
      mat A_k = 2*A_T+alpha*QQ;
      mat iAk_;
      bool ret = inv(A_k, iAk_);
      assert(ret);
      vec muk_ = iAk_*(A_T*(t1 + mu_T)+alpha*RQ); //TODO
      out = mvnrndex(muk_, iAk_, D);
    }
  }
  else if (i == K-1){
    vec tk_2 = times[K-2].pvec; 
    //calc least squares estimation of time nodes
    if (!BPTF){
      QQ = QQ + eDT;
      bool ret = itpp::ls_solve(QQ, pT*tk_2 + RQ, out);
      assert(ret); 
    }
    //calc equation A.8 in Xiong paper
    else {
      mat A_k = A_T+alpha*QQ;
      mat iAk_; 
      bool ret = inv(A_k, iAk_);
      assert(ret);
      vec muk_ = iAk_*(A_T*tk_2+alpha*RQ); //TODO
      out = mvnrndex(muk_, iAk_, D);
    }
  }
  else {
    vec tsum; 
    vec t1 = times[i-1].pvec; 
    vec t2 = times[i+1].pvec;
    tsum = t1+t2;
    //calc least squares estimation of time nodes
    if (!BPTF){
      QQ = QQ + 2*eDT;
      bool ret = itpp::ls_solve(QQ, pT*tsum + RQ, out);
      assert(ret);
    }
    //calc equation A.8 in Xiong paper
    else {
      mat A_k = 2*A_T+alpha*QQ;
      mat iAk_;
      bool ret = inv(A_k, iAk_);
      assert(ret);
      vec muk_ = iAk_*(A_T* tsum +alpha*RQ); //TODO
      out = mvnrndex(muk_, iAk_, D);
    }
  }

  times[i].pvec = out;
  if (debug && (i == 0|| i == K-1)) 
    std::cout<<times[i].pvec<<std::endl;
  assert(QQ.rows() == D && QQ.cols() == D);
  counter[BPTF_LEAST_SQUARES] += t.current_time();
  //}

  if (i == K-1){
    last_iter();
  }
          


}

// CALCULATE OBJECTIVE VALUE
double calc_obj(double res){
   
  double sumU = 0, sumV = 0, sumT = 0;
  timer t;
  t.start(); 
  for (int i=0; i< M; i++){
    const vertex_data * data = &g->vertex_data(i);
    sumU += sum_sqr(data->pvec);
  } 

  for (int i=M; i< M+N; i++){
    const vertex_data * data = &g->vertex_data(i);
    sumV += sum_sqr(data->pvec);
  } 


  mat T;
  if (tensor){
    T = zeros(D,K);
    for (int i=0; i<K; i++){
      vec tmp = times[i].pvec;
      sumT += pow(norm(tmp - vones, 2),2);
      T.set_col(i, tmp);
    }
    sumT *= pT;
  }
  counter[CALC_OBJ]+= t.current_time();
  
  double obj = (res + pU*sumU + pV*sumV + sumT + (tensor?trace(T*dp*T.transpose()):0)) / 2.0;
  return obj;
}


//SAVE FACTORS TO FILE
void export_uvt_to_file(){
 //saving output to file 
 mat U = zeros(M,D);
 mat V = zeros(N,D);
 mat T = zeros(K,D);
 for (int i=0; i< M+N; i++){ 
    vertex_data & data = g->vertex_data(i);
    if (i < M)
	U.set_row(i, data.pvec);
    else
	V.set_row(i-M, data.pvec);
 }

 if (tensor){ 
    for (int i=0; i<K; i++){
     	T.set_row(i, times[i].pvec);
    }
 } 

 char dfile[256] = {0};
 sprintf(dfile,"%s%d.out",infile.c_str(), D);
 it_file output(dfile);
 output << Name("User") << U;
 output << Name("Movie") << V;
  if (tensor){
    output << Name("Time") << T;
 }
 output.close();
}

//LOAD FACTORS FROM FILE
void import_uvt_from_file(){

 mat U,V,T;
 char dfile[256] = {0};
 sprintf(dfile,"%s%d.out",infile.c_str(), D);
 printf("Loading factors U,V,T from file\n");
 it_file input(dfile);
 input >> Name("User") >> U;
 input >> Name("Movie") >> V;
  if (tensor){
    input >> Name("Time") >> T;
 }
 input.close();
 //saving output to file 
 for (int i=0; i< M+N; i++){ 
    vertex_data & data = g->vertex_data(i);
    if (i < M)
        data.pvec = U.get_row(i); 
   else
        data.pvec =V.get_row(i-M);
 }

 if (tensor){ 
    for (int i=0; i<K; i++){
        times[i].pvec = T.get_row(i);
    }
 } 
 
}


 
/** 
 * ==== SETUP AND START
 */
void start(int argc, char ** argv) {
      
  command_line_options clopts;
  //clopts.scheduler_type = "round_robin";
  clopts.attach_option("debug", &debug, debug, "Display debug output. (optional)");
  clopts.attach_option("burn_in", &BURN_IN, BURN_IN, "burn-in period");
  clopts.attach_option("float", &FLOAT, FLOAT, "is data in float format?");
  clopts.attach_option("D", &D, D, "dmension of weight vector");
  clopts.attach_option("lambda", &LAMBDA, LAMBDA, "regularization weight");  
  clopts.attach_option("zero", &ZERO, ZERO, "support zero edges");  
  clopts.attach_option("scaling", &scaling, scaling, "time scaling factor (optional)");  
  clopts.attach_option("truncating", &truncating, truncating, "time truncation factor (optional)");  
  clopts.attach_option("savegraph", &savegraph, savegraph, "save graphs to file");  
  clopts.attach_option("loadgraph", &loadgraph, loadgraph, "load graphs to file");  
  clopts.attach_option("savefactors", &savefactors, savefactors, "save factors to file");  
  clopts.attach_option("loadfactors", &loadfactors, loadfactors, "load factors from file");  
  clopts.attach_option("stats", &stats, stats, "compute graph statistics");  
  clopts.attach_option("alpha", &alpha, alpha, "BPTF alpha (noise parameter)");  
  clopts.attach_option("regnormal", &regnormal, regnormal, "regular normalization? ");  
  clopts.attach_option("scalerating", &scalerating, scalerating, "scale rating value ");  
  clopts.attach_option("delayalpha", &delayalpha, delayalpha, "start sampling alpha the delayalpha round ");  
  clopts.attach_option("aggregatevalidation", &aggregatevalidation, aggregatevalidation, "aggregate training and validation into one dataset ");  
 
  gl_types::core glcore;
  assert(clopts.parse(argc-2, argv+2));

  if (delayalpha != 0 && (options != BPTF_TENSOR_MULT && options != BPTF_TENSOR))
	logstream(LOG_WARNING) << "Delaying alpha (sampling of noise level) is ignored in non-MCMC methods" << std::endl;

  if (BURN_IN != 10 && (options != BPTF_TENSOR_MULT && options != BPTF_TENSOR))
	logstream(LOG_WARNING) << "Markov chain burn in period is ignored in non-MCMC methods" << std::endl;



  //read the training data
  printf("loading data file %s\n", infile.c_str());
  if (!loadgraph){
    g=&glcore.graph();
    load_pmf_graph(infile.c_str(), g, TRAINING, glcore);

  //read the vlidation data (optional)
    printf("loading data file %s\n", (infile+"e").c_str());
    load_pmf_graph((infile+"e").c_str(),&validation_graph, VALIDATION, glcore);

  //read the test data (optional)
    printf("loading data file %s\n", (infile+"t").c_str());
    load_pmf_graph((infile+"t").c_str(),&test_graph, TEST, glcore);


    if (savegraph){
	printf("Saving .graph files\n");
	char filename[256];
        sprintf(filename, "%s%d.graph", infile.c_str(), D);
        std::ofstream fout(filename, std::fstream::binary);
        graphlab::oarchive oarc(fout);
	oarc << M << N << K << L << Le << Lt << D;
        oarc << *g << validation_graph << test_graph;
        printf("Done!\n");
        fout.close();
	exit(0);
    }

  } else {
    char filename[256];
    sprintf(filename, "%s%d.graph", infile.c_str(), D);
    std::ifstream fin(filename, std::fstream::binary);
    graphlab::iarchive iarc(fin);
    iarc >> M >> N >> K >> L >> Le >> Lt >> D;
    printf("Loading graph from file\n");
    iarc >> glcore.graph() >> validation_graph >> test_graph;
    g=&glcore.graph();
    printf("Matrix size is: USERS %dx MOVIES %dx TIME BINS %d D=%d\n", M, N, K, D);   
    printf("Creating %d edges (observed ratings)...\n", L);
  }
  

  if (loadfactors){
     import_uvt_from_file();
  }


  if (stats){
    calc_stats(TRAINING);
    calc_stats(VALIDATION);
    calc_stats(TEST);
    exit(0);
  }

  printf("setting regularization weight to %g\n", LAMBDA);
  pU=pV=LAMBDA;
  glcore.set_engine_options(clopts); 

  if (tensor)
    dp = GenDiffMat(K)*pT;
  if (debug)
    std::cout<<dp<<std::endl;

  std::vector<vertex_id_t> um;
  for (int i=0; i< M+N; i++)
    um.push_back(i);
 
  // add update function for user and movie nodes (tensor dims 1+2) 
  glcore.add_tasks(um, user_movie_nodes_update_function, 1);
  
  std::vector<vertex_id_t> tv;
  if (tensor){
    for (int i=M+N; i< M+N+K; i++)
      tv.push_back(i);
     // add update function for time nodes (tensor dim 3)
    glcore.add_tasks(tv, time_node_update_function, 1);
  }

  
  printf("%s for %s (%d, %d, %d):%d.  D=%d\n", BPTF?"BPTF":"PTF_ALS", tensor?"tensor":"matrix", M, N, K, L, D);
  
  if (!BPTF)
    printf("pU=%g, pV=%g, pT=%g, muT=%g, D=%d\n", pU, pV, pT, muT,D);  

  //if (BPTF)
  init_self_pot(); //init anyway
  init_pmf();

  double res, res2;
  double rmse =  calc_rmse(g, false, res);
  printf("complete. Obj=%g, TRAIN RMSE=%0.4f VALIDATION RMSE=%0.4f.\n", calc_obj(res), rmse, calc_rmse(&validation_graph, true, res2));


  if (BPTF){
    if (delayalpha < iiter)
    	sample_alpha(L);
    sample_U();
    sample_V();
    if (tensor) 
      sample_T();
  }

  /// Timing
  gt.start();
  g->finalize();  

  /**** START GRAPHLAB *****/
  glcore.start();

  // calculate final RMSE
  rmse =  calc_rmse_q(res);
  printf("Final result. Obj=%g, TRAIN RMSE= %0.4f VALIDATION RMSE= %0.4f.\n", calc_obj(res),  rmse, calc_rmse(&validation_graph, true, res2));
  
  if (infile == "kddcup" || infile == "kddcup2")
  	export_kdd_format(&test_graph, true);

  /**** POST-PROCESSING *****/
  double runtime = gt.current_time();
  printf("Finished in %lf \n", runtime);
      
  //timing counters
  for (int i=0; i<11; i++){
    if (counter[i] > 0)
    	printf("Counters are: %d) %s, %g\n",i, countername[i], counter[i]); 
   }

  export_uvt_to_file();
}


/**
 * read edges from file, with support with multiple edges between the same pair of nodes (in different times)
 */
template<typename edgedata>
int read_mult_edges(FILE * f, int nodes, testtype type, graph_type * _g, bool symmetry = false){
     
  //typedef typename graph::edge_data_type edge_data;
  bool * flags = NULL;
  if (options == BPTF_TENSOR_MULT || options == ALS_TENSOR_MULT){
    flags = new bool[nodes];
    memset(flags, 0, sizeof(bool)*nodes);
  }
 
  unsigned int e;
  int rc = fread(&e,1,4,f);
  assert(rc == 4);
  printf("Creating %d edges (observed ratings)...\n", e);
  assert(e>0);
  int total = 0;
  edgedata* ed = new edgedata[200000];
  int edgecount_in_file = e;
  while(true){
    //memset(ed, 0, 200000*sizeof(edge_float));
    rc = (int)fread(ed, sizeof(edgedata), _min(200000, edgecount_in_file - total), f);
    total += rc;

    for (int i=0; i<rc; i++){
#ifndef GL_NO_MULT_EDGES
      multiple_edges edges;
#endif
      edge_data edge;
      if (!ZERO)
	 assert(ed[i].weight != 0); // && ed[i].weight <= 5);
      assert((int)ed[i].from >= 1 && (int)ed[i].from <= nodes);
      assert((int)ed[i].to >= 1 && (int)ed[i].to <= nodes);
      assert((int)ed[i].to != (int)ed[i].from);
      edge.weight = (double)ed[i].weight;
      if (scalerating)
	edge.weight /= scalerating;
      edge.time = (int)((ed[i].time -1-truncating)/scaling);
 
      std::pair<bool, edge_id_t> ret;
      if (options != BPTF_TENSOR_MULT && options != ALS_TENSOR_MULT){//no support for specific edge returning on different times
        ret.first = false;
      }
      else if (flags[(int)ed[i].from-1] == true && flags[(int)ed[i].to-1] == true){
        ret = _g->find((int)ed[i].from-1, (int)ed[i].to-1);
      }
      else ret.first = false;

      if (ret.first == false){
#ifndef GL_NO_MULT_EDGES
        edges.medges.push_back(edge); 
        _g->add_edge((int)ed[i].from-1, (int)ed[i].to-1, edges); // Matlab export has ids starting from 1, ours start from 0
#else
	_g->add_edge((int)ed[i].from-1, (int)ed[i].to-1, edge);
#endif
        if (type == VALIDATION && aggregatevalidation)//add validation edges into training dataset as well
#ifndef GL_NO_MULT_EDGES          
           g->add_edge((int)ed[i].from-1, (int)ed[i].to-1, edges); // Matlab export has ids starting from 1, ours start from 0
#else
	   g->add_edge((int)ed[i].from-1, (int)ed[i].to-1, edge);
#endif
        if (options == BPTF_TENSOR_MULT || options == ALS_TENSOR_MULT){
          flags[(int)ed[i].from-1] = true;
          flags[(int)ed[i].to-1] = true;
        }
      }
#ifndef GL_NO_MULT_EDGES
      else {
        _g->edge_data(ret.second).medges.push_back(edge);
      }
#endif
    } 
    printf(".");
    fflush(0);
    if (rc == 0 || total >= edgecount_in_file)
      break;
  }
  assert(total == (int)e);
  delete [] ed; ed = NULL;
  if (flags != NULL)
    delete[] flags;
  return e;
}

//go over all ratings and count how ratings for each node (user/movie)
void count_all_edges(graph_type * _g){
    for (int i=0; i<M+N; i++){
        vertex_data &vdata = _g->vertex_data(i);
        if (i < M)
          vdata.num_edges = count_edges(_g->out_edge_ids(i));
        else
          vdata.num_edges = count_edges(_g->in_edge_ids(i));
     }
}
/* function that reads the tensor from file */
void load_pmf_graph(const char* filename, graph_type * _g, testtype data_type,gl_types::core & glcore) {

  printf("Loading %s %s\n", filename, testtypename[data_type]);
  FILE * f = fopen(filename, "r");
  if (data_type!=TRAINING && f == NULL){//skip optional files, if the file is missing
    printf("skipping file\n");
    return;
  }

  if(data_type==TRAINING && f== NULL){
	logstream(LOG_ERROR) << " can not find input file. aborting " << std::endl;
	exit(1);
  }

  int _M,_N,_K;
  fread(&_M,1,4,f);//movies
  fread(&_N,1,4,f);//users/
  fread(&_K,1,4,f);//time
  assert(_K>=1);
  assert(_M>=1 && _N>=1); 
  if (data_type != TRAINING && M != _M)
	logstream(LOG_WARNING) << " wrong number of users: " << _M << " instead of " << M << " in training file" << std::endl;
  if (data_type != TRAINING && N != _N)
	logstream(LOG_WARNING) << " wrong number of movies: " << _N << " instead of " << M << " in training file" << std::endl;
  if (data_type != TRAINING && K != _K)
	logstream(LOG_WARNING) << " wrong number of time bins: " << _K << " instead of " << K << " in training file" << std::endl;

  if (data_type == TRAINING){
  	M=_M; N= _N; K= _K;
	if (infile == "kddcup" || infile == "kddcup2")// DB: ugly - kdd cup data has more time bins for test data than in training data. can fix this buy setting the time bins in training data to 6649.
		K=6649;
        K=ceil((K-truncating)/scaling);
  }

  //if (data_type==TRAINING)
  printf("Matrix size is: USERS %d MOVIES %d TIME BINS %d\n", M, N, K);
 
  vertex_data vdata;

  // add M movie nodes (tensor dim 1)
  for (int i=0; i<M; i++){
    
    vdata.pvec = debug ? (itpp::ones(D)*0.1) : (itpp::randu(D)*0.1);
    _g->add_vertex(vdata);
    if (debug && (i<= 5 || i == M-1))
      debug_print_vec("U: ", vdata.pvec, D);
  }
  
  // add N user node (tensor dim 2) 
  for (int i=0; i<N; i++){
    vdata.pvec = debug ? (itpp::ones(D)*0.1) : (itpp::randu(D)*0.1);
    _g->add_vertex(vdata);
    if (debug && (i<=5 || i==N-1))
      debug_print_vec("V: ", vdata.pvec, D);
  }
  
  if (data_type==TRAINING && tensor){
    //init times
    times = new vertex_data[K];
    vec tones = itpp::ones(D)*(K==1?1:0.1);
    //add T time node (tensor dim 3)
    for (int i=0; i<K; i++){
      times[i].pvec =tones;
      _g->add_vertex(times[i]);
      if (debug && (i <= 5 || i == K-1))
        debug_print_vec("T: ", times[i].pvec, D);
    }
  }
  
  // read tensor non zero edges from file
  int val = 0; 
  if (!FLOAT) 
	val = read_mult_edges<edge_double>(f, M+N, data_type, _g);
  else val = read_mult_edges<edge_float>(f,M+N, data_type, _g);

  switch(data_type){
	case TRAINING: 
		L= val; 
		break;
        case VALIDATION: 
		Le = val; 
		if (aggregatevalidation)
			L+=  Le; //add edges of validation dataset into the training data set as well.
		break; 
        case TEST: 
		Lt = val; 
		break;
  }  

  if (data_type==TRAINING && tensor && K>1) 
    edges = new std::vector<edge_id_t>[K]();

  //verify edges
  for (int i=M; i < M+N; i++){
    foreach(graphlab::edge_id_t eid, _g->in_edge_ids(i)){          
#ifndef GL_NO_MULT_EDGES      
      multiple_edges & tedges= _g->edge_data(eid);
#endif
      int from = _g->source(eid);
      int to = _g->target(eid);
      assert(from < M);
      assert(to >= M && to < M+N);

#ifndef GL_NO_MULT_EDGES
      for (int j=0; j< (int)tedges.medges.size(); j++){
        edge_data & data= tedges.medges[j];
#else
      edge_data data = _g->edge_data(eid);
#endif
	if (!ZERO)
        	assert(data.weight != 0);  
        assert(data.time < K);
  
        if (K > 1 && data_type==TRAINING && tensor)
          edges[(int)data.time].push_back(eid);
#ifndef GL_NO_MULT_EDGES      
        }
#endif
    }
  }
  
 //store number of edges for each node 
 if (data_type == TRAINING || (aggregatevalidation && data_type == VALIDATION)){
 	count_all_edges(g);
  }

 
  //verify correct number of edges encourntered
  if (data_type==TRAINING && tensor && K>1){
    int cnt = 0;
    for (int i=0; i<K; i++){
      cnt+= edges[i].size();
    }
    assert(cnt == L);
  }

  fclose(f);
}


//main function 
int main(int argc,  char *argv[]) {

  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logstream(LOG_INFO)<< "PMF/ALS Code written By Danny Bickson, CMU\nSend bug reports and comments to danny.bickson@gmail.com\n";
#ifdef GL_NO_MULT_EDGES
  logstream(LOG_WARNING)<<"Code compiled with GL_NO_MULT_EDGES flag - this mode does not support multiple edges between user and movie in different times\n";
#endif
#ifdef GL_NO_MCMC
  logstream(LOG_WARNING)<<"Code compiled with GL_NO_MCMC flag - this mode does not support MCMC methods.\n";
#endif


  if (argc < 3){
       logstream(LOG_ERROR) <<  "Not enough input arguments. Usage is ./pmf <input file name> <run mode> \n \tRun mode are: \n\t0 = Matrix factorization using alternating least squares \n\t1 = Matrix factorization using MCMC procedure \n\t2 = Tensor factorization using MCMC procedure, single edge exist between user and movies \n\t3 = Tensor factorization, using MCMC procedure with support for multiple edges between user and movies in different times \n\t4 = Tensor factorization using alternating least squars\n";  
       exit(1);
   }
   
   // select run mode (type of algorithm)
  infile = argv[1];
  options = (runmodes)atoi(argv[2]);
  printf("Setting run mode %s\n", runmodesname[options]);
  switch(options){
  // iterative matrix factorization using alternating least squares
  // or SVD ++
  case ALS_MATRIX:
  case SVD_PLUS_PLUS:
    tensor = false; BPTF = false;
    break;
    // MCMC tensor factorization
  case BPTF_TENSOR:
    tensor = true; BPTF = true;
   break;
    //MCMC matrix factorization
  case BPTF_MATRIX:
    tensor = false; BPTF = true;
    break;
    // tensor factorization , allow for multiple edges between user and movie in different times
  case BPTF_TENSOR_MULT:
    tensor = true; BPTF = true;
    break;
   // tensor factorization
  case ALS_TENSOR_MULT:
    tensor = true; BPTF = false;
    break;
  default:
    assert(0);
  }

  logger(LOG_INFO, "%s starting\n",runmodesname[options]);
#ifdef GL_NO_MCMC
  if (BPTF)
    logstream(LOG_ERROR) << "Can not run MCMC method with GL_NO_MCMC flag. Please comment flag and recompile\n";
#endif

#ifdef GL_NO_MULT_EDGES
  if (options == ALS_TENSOR_MULT || options == BPTF_TENSOR_MULT)
    logstream(LOG_ERROR) << "Can not have support for multiple edges with GL_NO_MULT_EDGES flag. Please comment flag and recompile\n";
#endif  

   start(argc, argv);
}

// calc statistics about matrix/tensor and exit  
void calc_stats(testtype type){
   graph_type * gr = NULL;
   switch(type){
	   case TRAINING: gr = g; break;
	   case VALIDATION: gr = &validation_graph; break;
	    case TEST: gr = &test_graph; break;
   }

   if (gr->num_vertices() == 0){
	printf("%s is missing, skipping data\n", testtypename[type]);
        return;
   } 

   if (tensor && type == TRAINING){
	int firsttimeused=-1;
	int lasttimeused=-1;
	for (int i=0; i<K; i++){
		if (edges[i].size() > 0)
		   firsttimeused = i;
	}
	for (int i=K-1; i>=0; i--){
		if (edges[i].size() > 0)
		   lasttimeused = i;

	}
	printf("Out of total %d time components, first used is %d, last used is %d\n", K, firsttimeused, lasttimeused);
  }	

  double avgval=-1, minval=1e100, maxval=-1e100;
  double avgtime=-1, mintime=1e100, maxtime=-1e100;
  double minV=1e100, maxV=-1e100, minU=1e100, maxU=-1e100;
  int moviewithoutedges = 0;
  int userwithoutedges = 0;
  int numedges = 0;
  for (int i=M; i< M+N; i++){ 
    const vertex_data * data = &gr->vertex_data(i);
         if (itpp::min(data->pvec) < minU)
		minU = itpp::min(data->pvec);
	 if (itpp::max(data->pvec) > maxU)
		maxU = itpp::max(data->pvec);
	if (gr->in_edge_ids(i).size() == 0)
	moviewithoutedges++;
    foreach(edge_id_t iedgeid, gr->in_edge_ids(i)) {
#ifndef GL_NO_MULT_EDGES            
         multiple_edges & edges = gr->edge_data(iedgeid);
#endif
         //vertex_data * pdata = &gr->vertex_data(gr->source(iedgeid)); 
#ifndef GL_NO_MULT_EDGES        
         for (int j=0; j< (int)edges.medges.size(); j++){     
		edge_data & data = edges.medges[j];
#else
		edge_data & data = gr->edge_data(iedgeid);
#endif
		numedges++;
		avgval += data.weight;
		avgtime += data.time;
		if (data.weight<minval)
		   minval=data.weight;
		if (data.time <mintime)
		   mintime = data.time;
		if (data.weight>maxval)
		   maxval=data.weight;
		if (data.time > maxtime)
		   maxtime =data.time;
#ifndef GL_NO_MULT_EDGES	        
	 }  
#endif
	
     }
 }
  for (int i=0; i< M; i++){ 
    const vertex_data * data = &gr->vertex_data(i);
    if (itpp::min(data->pvec) < minV)
	minV = itpp::min(data->pvec);
	if (itpp::max(data->pvec) > maxV)
		maxV = itpp::max(data->pvec);
	 
    if (gr->out_edge_ids(i).size() == 0)
	userwithoutedges++;
 }
 
 avgval /= numedges;
 avgtime /= numedges; 
 printf("%s Avg matrix value %g min val %g max value %g\n", testtypename[type],avgval, minval, maxval);
 printf("%s Avg time value %g min val %g max value %g\n", testtypename[type], avgtime, mintime, maxtime);
 printf("%s User without edges: %d movie without edges: %d\n", testtypename[type], userwithoutedges, moviewithoutedges);
 printf("%s Min V: %g Max V: %g Min U: %g, Max U: %g \n", testtypename[type], minV, maxV, minU, maxU);

 //verify we did not miss any ratings
 switch(type){
 	case TRAINING: assert(numedges==L); break;
	 case VALIDATION: assert(numedges==Le); break;
	 case TEST: assert(numedges==Lt); break;
 }
 }

//
//The input prediction file should contain 6005940 lines, corresponding
//to the 6005940 user-item pairs in the test set.
//Each line contains a predicted score (a real number between 0 and 100).
//The generated output file can be submitted to the KDD-Cup'11 evaluation
//system.


void export_kdd_format(graph_type * _g, bool dosave) {

        bool debugkdd = true;
        assert(_g != NULL);
        if (!dosave)
		assert(BPTF);	

	FILE * outFp = NULL;
        if (dosave){
                printf("Exporting KDD cup test graph: %s\n", (infile+"t.kdd.out").c_str());
		outFp = fopen((infile+"t.kdd.out").c_str(), "w");
		assert(outFp);
	}
	const int ExpectedTestSize = 6005940;

	int lineNum = 0;
	double prediction;
	double sumPreds=0;


     for (int i=0; i< M; i++){ //TODO: optimize to start from N?
       vertex_data & data = g->vertex_data(i);
       foreach(edge_id_t iedgeid, _g->out_edge_ids(i)) {
#ifndef GL_NO_MULT_EDGES            
         multiple_edges & edges = _g->edge_data(iedgeid);
#endif
         vertex_data & pdata = g->vertex_data(_g->target(iedgeid)); 
#ifndef GL_NO_MULT_EDGES
         for (int j=0; j< (int)edges.medges.size(); j++){  
           edge_data & edge = edges.medges[j];
#else
	   edge_data & edge = _g->edge_data(iedgeid);
#endif     
           if (!ZERO)
           	assert(edge.weight != 0);

           prediction = 0;
           rmse(data.pvec, pdata.pvec, tensor? (&times[(int)edge.time].pvec):NULL, D, edge.weight, prediction);
#ifndef GL_NO_MCMC 
          if (BPTF && iiter > BURN_IN){
             edge.avgprd += prediction;
             //add = powf((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
              prediction = (edge.avgprd / (iiter - BURN_IN));
           }
#endif
           if (debugkdd && (i== M || i == M+N-1))
             cout<<lineNum<<") prediction:"<<prediction<<endl; 
  	   if (prediction<0)
		prediction=0;
	   else if (prediction>100)
		prediction=100; 
          unsigned char roundScore = (unsigned char)(2.55*prediction + 0.5); 
          if (dosave)
	  	fwrite(&roundScore,1,1,outFp);
	  sumPreds += prediction;

 	  lineNum++; 
#ifndef GL_NO_MULT_EDGES          
         }
#endif
       }
     }

   assert(lineNum==Lt); 
   if (lineNum!= ExpectedTestSize)
	logstream(LOG_WARNING) << "KDD test data has wrong length." << " current length is: " << Lt << " correct length " << ExpectedTestSize << std::endl;
 
  if (dosave){
    fclose(outFp);
    fprintf(stderr, "**Completed successfully (mean prediction: %lf)**\n",sumPreds/ExpectedTestSize);
  }
}













#include <graphlab/macros_undef.hpp>
