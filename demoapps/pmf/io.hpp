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
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 */


#ifndef _IO_HPP
#define _IO_HPP

#include "stats.hpp"
#include <graphlab/macros_def.hpp>

extern double scalerating;
extern double truncating;
extern double scaling;
extern int L,M,N,L,Lt,Le;
extern mat U,V,T;
extern bool aggregatevalidation;
extern runmodes algorithm;
extern double globalMean[3];
extern string implicitratingtype;

void add_implicit_edges(graph_type * g);


void fill_factors_uvt(){
 if (algorithm != LANCZOS){
   U = zeros(M,D);
   V = zeros(N,D);
   for (int i=0; i< M+N; i++){ 
      vertex_data & data = g->vertex_data(i);
      if (i < M)
          U.set_row(i, data.pvec);
      else
        V.set_row(i-M, data.pvec);
   }

   if (tensor){ 
     T = zeros(K,D);
     for (int i=0; i<K; i++){
        T.set_row(i, times[i].pvec);
     }
    } 
  }
} 


template<typename edgedata>
int read_mult_edges(FILE * f, int nodes, testtype type, graph_type * _g, bool symmetry = false);
 //
//The input prediction file should contain 6005940 lines, corresponding
//to the 6005940 user-item pairs in the test set.
//Each line contains a predicted score (a real number between 0 and 100).
//The generated output file can be submitted to the KDD-Cup'11 evaluation
//system.

//write an output vector to file
void write_vec(FILE * f, int len, double * array){
  assert(f != NULL && array != NULL);
  fwrite(array, len, sizeof(double), f);
}



void export_kdd_format(graph_type * _g, testtype type, bool dosave) {

  bool debugkdd = true;
  assert(_g != NULL);
  if (!dosave)
    assert(BPTF);	

    FILE * outFp = NULL;
    if (dosave){
      printf("Exporting KDD cup %s graph: %s\n", testtypename[type], (infile+"t.kdd.out").c_str());
      outFp = fopen((infile+"t.kdd.out").c_str(), "w");
      assert(outFp);
    }
    const int ExpectedTestSize = 6005940;

    int lineNum = 0;
    float prediction;
    double sumPreds=0;


    for (int i=0; i< M; i++){ //TODO: optimize to start from N?
      vertex_data & data = g->vertex_data(i);


#ifdef GL_SVD_PP
       int n = data.num_edges; //+1.0 ? //regularization
       data.weight = zeros(D);
       foreach(edge_id_t oedgeid, g->out_edge_ids(i)) {
         vertex_data & movie = g->vertex_data(g->target(oedgeid)); 
	 data.weight += movie.weight;
       }
       float usrnorm = float(1.0/sqrt(n));
       data.weight *= usrnorm;

#endif


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
          predict(data, 
                  pdata, 
                  algorithm == WEIGHTED_ALS ? &edge : NULL, 
                  tensor? (&times[(int)edge.time]):NULL, 
                  edge.weight, 
                  prediction);
#ifndef GL_NO_MCMC 
          if (BPTF && iiter > BURN_IN){
             edge.avgprd += prediction;
             //add = powf((edge.avgprd / (iiter - BURN_IN)) - edge.weight, 2);
              prediction = (edge.avgprd / (iiter - BURN_IN));
           }
#endif
          if (prediction<minval)
	     prediction=minval;
	  else if (prediction>maxval)
	     prediction=maxval; 
            
	  if (scalerating != 1)
	    prediction *= scalerating;
	  if (debugkdd && (i== 0 || i == M))
            cout<<lineNum<<") prediction:"<<prediction<<endl; 
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
   switch(type){
     case TEST:
        if (lineNum!= ExpectedTestSize)
  	   logstream(LOG_WARNING) << "KDD test data has wrong length." << " current length is: " << Lt << " correct length " << ExpectedTestSize << std::endl;
           assert(lineNum==Lt); 
        break;
     case VALIDATION:
       assert(lineNum==Le);
       break;
     case TRAINING:
       assert(false);
       break; 
  }

  if (dosave){
    fclose(outFp);
    fprintf(stderr, "**Completed successfully (mean prediction: %lf)**\n",sumPreds/lineNum);
  }
}



//OUTPUT: SAVE FACTORS U,V,T to a binary file

// FORMAT:  M N K D (4 x ints)
// MATRIX U ( M x D doubles)
// MATRIX V ( N x D doubles)
// MATRIX K ( K x D doubles - optional, only for tensor)
// TOTAL FILE SIZE: 4 ints + (M+N+K)*D - for tensor
//                  4 ints + (M+N)*D - for matrix
void export_uvt_to_binary_file(){

  fill_factors_uvt();

  char dfile[256] = {0};
  sprintf(dfile,"%s%d.out",infile.c_str(), D);
  FILE * f = fopen(dfile, "w");
  assert(f!= NULL);

  int rc = fwrite(&M, 1, 4, f);
  assert(rc == 4);
  rc = fwrite(&N, 1, 4, f);
  assert(rc == 4);
  rc = fwrite(&K, 1, 4, f);
  assert(rc == 4);
  rc = fwrite(&D, 1, 4, f);
  assert(rc == 4);

  write_vec(f, M*D, U._data());
  write_vec(f, N*D, V._data());
  if (tensor)
    write_vec(f, K*D, T._data());

  fclose(f); 

}


//OUTPUT: SAVE FACTORS U,V,T TO IT++ FILE
void export_uvt_to_itpp_file(){

  fill_factors_uvt();

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


/* function that reads the tensor from file */
/* Input format is:
 * M - number of users (int)
 * N - number of movies (int)
 * K - number of time bins (int), in case of weighted ALS this value is ignored
 * e - number of edges (int)
 * A list of edges in the format
 * [from] [to] [ time] [weight]  (4 floats)
 * where [from] is an integeter from 1 to M
 * [to] is an interger from 1 to N
 * [time] is an integer from 1 to K (for weighted ALS this is the weight, which is float)
 * [weight] - this is the rating, which is float. Rating is assumed non-zero unless the --zero=true flas is on 
 */
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
  int rc = fread(&_M,1,4,f);//movies
  assert(rc==4); 
  rc=fread(&_N,1,4,f);//users/
  assert(rc==4); 
  rc=fread(&_K,1,4,f);//time
  assert(rc==4); 
  assert(_K>=1);
  assert(_M>=1 && _N>=1); 
  if (data_type == TRAINING){
  	M=_M; N= _N; K= _K;
	if (infile == "kddcup" || infile == "kddcup2")// DB: ugly - kdd cup data has more time bins for test data than in training data. can fix this buy setting the time bins in training data to 6649.
		K=6649;
     K=ceil((K-truncating)/scaling);
  }

 if (data_type != TRAINING && M != _M)
	logstream(LOG_WARNING) << " wrong number of users: " << _M << " instead of " << M << " in " << testtypename[data_type] << std::endl;
  if (data_type != TRAINING && N != _N)
	logstream(LOG_WARNING) << " wrong number of movies: " << _N << " instead of " << N << " in " << testtypename[data_type] << std::endl;
  if (data_type != TRAINING && K != _K)
	logstream(LOG_WARNING) << " wrong number of time bins: " << _K << " instead of " << K << " in " << testtypename[data_type] <<std::endl;

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
  else 
     val = read_mult_edges<edge_float>(f,M+N, data_type, _g);

  switch(data_type){
    case TRAINING: 
      L = val; break;
    
    case VALIDATION: 
      Le = val; 
      if (aggregatevalidation)
	L+=  Le; //add edges of validation dataset into the training data set as well.
	break; 
     
    case TEST: 
      Lt = val; break;
  }  

  if (data_type==TRAINING && tensor && K>1) 
    edges = new std::vector<edge_id_t>[K]();

  //verify edges
  for (int i=M; i < M+N; i++){
    foreach(graphlab::edge_id_t eid, _g->in_edge_ids(i)){          
#ifndef GL_NO_MULT_EDGES      
     const  multiple_edges & tedges= _g->edge_data(eid);
#endif
      int from = _g->source(eid);
      int to = _g->target(eid);
      assert(from < M);
      assert(to >= M && to < M+N);

#ifndef GL_NO_MULT_EDGES
      for (int j=0; j< (int)tedges.medges.size(); j++){
        const edge_data & data= tedges.medges[j];
#else
      const edge_data & data = _g->edge_data(eid);
#endif
	if (!ZERO)
          assert(data.weight != 0);  
        if (algorithm != WEIGHTED_ALS)
          assert(data.time < K);
  
        if (K > 1 && data_type==TRAINING && tensor)
          edges[(int)data.time].push_back(eid);
#ifndef GL_NO_MULT_EDGES      
        }
#endif
    }
  }
  fclose(f);
  
  //add implicit edges if requested
  if (data_type == TRAINING && implicitratingtype != "none")
     add_implicit_edges(_g);

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

}



/**
 * read edges from file, with support with multiple edges between the same pair of nodes (in different times)
 */
template<typename edgedata>
int read_mult_edges(FILE * f, int nodes, testtype type, graph_type * _g, bool symmetry = false){
     
  //typedef typename graph::edge_data_type edge_data;
  bool * flags = NULL;
  if (algorithm == BPTF_TENSOR_MULT || algorithm == ALS_TENSOR_MULT){
    flags = new bool[nodes];
    memset(flags, 0, sizeof(bool)*nodes);
  }

  int matlab_offset_user_movie = 1; //matlab array start from 1
  int matlab_offset_time = 1; //matlab arrays start from 1
  if (algorithm == WEIGHTED_ALS)
    matlab_offset_time = 0; //for weighted ALS there are no time bins which are integers, so no need to convert them

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
    rc = (int)fread(ed, sizeof(edgedata), std::min(200000, edgecount_in_file - total), f);
    total += rc;

    //go over each rating (edges)
    for (int i=0; i<rc; i++){
#ifndef GL_NO_MULT_EDGES
      multiple_edges edges;
#endif
      edge_data edge;
      if (!ZERO) //usually we do not allow zero ratings, unless --zero=true flag is set.
	 assert(ed[i].weight != 0); 
      //verify node ids are in allowed range
      assert((int)ed[i].from >= matlab_offset_user_movie && (int)ed[i].from <= nodes);
      assert((int)ed[i].to >= matlab_offset_user_movie && (int)ed[i].to <= nodes);
      //no self edges
      assert((int)ed[i].to != (int)ed[i].from);
      edge.weight = (double)ed[i].weight;
    
      //if sacling of rating values is requested to it here.
      if (scalerating != 1.0)
	     edge.weight /= scalerating;
      globalMean[type] += edge.weight;
     
      //if scaling of time bins request do it here
      double time  = ((ed[i].time - matlab_offset_time - truncating)/(double)scaling);
      edge.time = time;
      //assert weights in WALS are not zero (no sense to give zero weight)
      if (algorithm == WEIGHTED_ALS && !ZERO)
         assert(edge.time != 0);

      std::pair<bool, edge_id_t> ret;
      if (algorithm != BPTF_TENSOR_MULT && algorithm != ALS_TENSOR_MULT){//no support for multple edges (ratings) of the same user - item pair at different times
        ret.first = false;
      }
      else if (flags[(int)ed[i].from-matlab_offset_user_movie] == true && flags[(int)ed[i].to-matlab_offset_user_movie] == true){
        ret = _g->find((int)ed[i].from-matlab_offset_user_movie, (int)ed[i].to-matlab_offset_user_movie);
      }
      else ret.first = false;

      if (ret.first == false){
#ifndef GL_NO_MULT_EDGES
        edges.medges.push_back(edge); 
        _g->add_edge((int)ed[i].from-matlab_offset_user_movie, (int)ed[i].to-matlab_offset_user_movie, edges); // Matlab export has ids starting from 1, ours start from 0
#else
	_g->add_edge((int)ed[i].from-matlab_offset_user_movie, (int)ed[i].to-matlab_offset_user_movie, edge);
#endif
        if (type == VALIDATION && aggregatevalidation)//add validation edges into training dataset as well
#ifndef GL_NO_MULT_EDGES          
           g->add_edge((int)ed[i].from-matlab_offset_user_movie, (int)ed[i].to-matlab_offset_user_movie, edges); // Matlab export has ids starting from 1, ours start from 0
#else
	   g->add_edge((int)ed[i].from-matlab_offset_user_movie, (int)ed[i].to-matlab_offset_user_movie, edge);
#endif
        if (algorithm == BPTF_TENSOR_MULT || algorithm == ALS_TENSOR_MULT){
          flags[(int)ed[i].from-matlab_offset_user_movie] = true;
          flags[(int)ed[i].to-matlab_offset_user_movie] = true;
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
  globalMean[type] /= e;
  delete [] ed; ed = NULL;
  if (flags != NULL)
    delete[] flags;
  return e;
}



#include <graphlab/macros_undef.hpp>
#endif
