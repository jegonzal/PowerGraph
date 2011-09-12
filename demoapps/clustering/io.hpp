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
#include "../../libs/matrixmarket/mmio.h" //matrix market format support
#include <graphlab/macros_def.hpp>

extern advanced_config config;
extern problem_setup ps;
extern const char * inittypenames[];
void init();

template<typename edgedata>
int read_edges(FILE * f, int nodes, graph_type * _g);


void fill_output(){
  
   if (ac.algorithm == LDA)
	return;
  
   ps.output_clusters = zeros(ps.K, ps.N);
   for (int i=0; i<ps.K; i++)
     ps.output_clusters.set_row(i, ps.clusts.cluster_vec[i].location);
     
   int cols = 1;
   if (ac.algorithm == K_MEANS_FUZZY)
	cols = ac.K;
   ps.output_assignements = zeros(ps.M, cols);
     for (int i=0; i< ps.M; i++){ 
        const vertex_data & data = ps.g->vertex_data(i);
        if (ac.algorithm == K_MEANS){
          ps.output_assignements.set(i,0, data.current_cluster);
        } 
	else if (ac.algorithm == K_MEANS_FUZZY) 
          ps.output_assignements.set_row(i, data.distances);
     }
} 



//write an output vector to file
void write_vec(FILE * f, int len, double * array){
  assert(f != NULL && array != NULL);
  fwrite(array, len, sizeof(double), f);
}





//OUTPUT: SAVE FACTORS U,V,T to a binary file

// FORMAT:  M N K D (4 x ints)
// MATRIX U ( M x D doubles)
// MATRIX V ( N x D doubles)
// MATRIX K ( K x D doubles - optional, only for ps.tensor)
// TOTAL FILE SIZE: 4 ints + (M+N+K)*D - for ps.tensor
//                  4 ints + (M+N)*D - for matrix
void export_to_binary_file(){

  fill_output();

  char dfile[256] = {0};
  sprintf(dfile,"%s%d.out",ac.datafile.c_str(),ps.K);
  FILE * f = fopen(dfile, "w");
  assert(f!= NULL);

  int rc = fwrite(&ps.M, 1, 4, f);
  assert(rc == 4);
  rc = fwrite(&ps.N, 1, 4, f);
  assert(rc == 4);
  rc = fwrite(&ps.K, 1, 4, f);
  assert(rc == 4);
  write_vec(f, ps.K*ps.N, ps.output_clusters._data());
  write_vec(f, ps.M, ps.output_assignements._data());
  fclose(f); 

}


//OUTPUT: SAVE FACTORS U,V,T TO IT++ FILE
void export_to_itpp_file(){

  fill_output();

  char dfile[256] = {0};
  sprintf(dfile,"%s%d.out",ac.datafile.c_str(), ac.D);
  it_file output(dfile);
  output << Name("Clusters") << ps.output_clusters;
  output << Name("Assignments") << ps.output_assignements;
  output.close();
}


void export_to_matrixmarket(){
  fill_output();
  char dfile[256] = {0};
  sprintf(dfile,"%s%d",ac.datafile.c_str(), ps.K);
  save_matrix_market_format(dfile);  
}

//LOAD FACTORS FROM FILE
void import_from_file(){

 mat clusters;
 vec assignments;
 char dfile[256] = {0};
 sprintf(dfile,"%s%d.out",ac.datafile.c_str(), ac.D);
 printf("Loading clusters from file\n");
 it_file input(dfile);
 input >> Name("Clusters") >> clusters;
 input >> Name("Assignments") >> assignments;
 input.close();
 //saving output to file 
 for (int i=0; i< ps.M; i++){ 
    vertex_data & data = ps.g->vertex_data(i);
    data.current_cluster = assignments[i]; 
 }
 for (int i=0; i<ps.K; i++){
    ps.clusts.cluster_vec[i].location = clusters.get_row(i);
 }
}


void add_vertices(graph_type * _g){
  assert(ps.K > 0);
  vertex_data vdata;
  // add M movie nodes (ps.tensor dim 1)
  for (int i=0; i<ps.M; i++){
    switch (ps.init_type){
       case INIT_RANDOM:
         vdata.current_cluster = randi(0, ps.K-1);
         break;
       
       case INIT_ROUND_ROBIN:
	 vdata.current_cluster = i % ps.K;
         break;

       case INIT_KMEANS_PLUS_PLUS: //is done later
       case INIT_RANDOM_CLUSTER:
	 vdata.current_cluster = -1;
	 break;
    }

    //vdata.datapoint.set_size(ps.N);
    if (ps.algorithm == K_MEANS_FUZZY)
	vdata.distances = zeros(ps.K);

    _g->add_vertex(vdata);
    if (ac.debug && (i<= 5 || i == ps.M-1))
       std::cout<<"node " << i <<" initial assignment is: " << vdata.current_cluster << std::endl; 
 }
  
}


/* function that reads the problem from file */
/* Input format is:
 * M - number of matrix rows
 * N - number of matrix cols
 * zero - unused
 * e - number of edges (int)
 * A list of edges in the format
 * [from] [to] [weight]  (3 floats)
 * where [from] is an integeter from 1 to M
 * [to] is an interger from 1 to N
 * [weight] float
 */
void load_graph(const char* filename, graph_type * _g, gl_types::core & glcore) {


  if (ac.matrixmarket){
      printf("Loading Matrix Market file %s\n", filename);
      load_matrix_market(filename, _g);
      return;
  }

  printf("Loading %s\n", filename);
  FILE * f = fopen(filename, "r");
  if(f== NULL){
	logstream(LOG_ERROR) << " can not find input file. aborting " << std::endl;
	exit(1);
  }

  int _M,_N,_K;
  int rc = fread(&_M,1,4,f);//matrix rows
  assert(rc==4); 
  rc=fread(&_N,1,4,f);//matrix cols
  assert(rc==4); 
  rc=fread(&_K,1,4,f);//unused
  assert(rc==4);
  if (!ac.supportgraphlabcf) 
    assert(_K == 0); 
  else assert(_K >= 1);
  assert(_M>=1 && _N>=1); 
  


  ps.M=_M; ps.N= _N;
  init();
  add_vertices(_g);
 
  // read tensor non zero edges from file
  
 int val;
  if (!ac.supportgraphlabcf)
    val = read_edges<edge_float>(f,ps.N, _g);
  else {
      if (ac.FLOAT)
    val = read_edges<edge_float_cf>(f,ps.N,_g);
    else val = read_edges<edge_double_cf>(f,ps.N,_g);
  }assert(val == ps.L);

  fclose(f);
}


/**
 * read edges from file, with support with multiple edges between the same pair of nodes (in different times)
 */
template<typename edgedata>
int read_edges(FILE * f, int column_dim, graph_type * _g){
     
  int matlab_offset = 1; //matlab array start from 1

  unsigned int e;
  int rc = fread(&e,1,4,f);
  assert(rc == 4);
  printf("Creating %d edges (observed ratings)...\n", e);
  assert(e>0);

  int total = 0;
  edgedata* ed = new edgedata[200000];
  int edgecount_in_file = e;
  while(true){
    rc = (int)fread(ed, sizeof(edgedata), std::min(200000, edgecount_in_file - total), f);
    total += rc;

    //go over each rating (edges)
    for (int i=0; i<rc; i++){
      if (!ac.zero) //usually we do not allow zero entries, unless --zero=true flag is set.
	 assert(ed[i].weight != 0); 
      //verify node ids are in allowed range
      assert((int)ed[i].from >= matlab_offset && (int)ed[i].from <= ps.M);
      if (ac.supportgraphlabcf)
        ed[i].to -= ps.M;
      assert((int)ed[i].to >= matlab_offset && (int)ed[i].to <= ps.N);
      //no self edges
      //assert((int)ed[i].to != (int)ed[i].from);
    
      //if sacling of rating values is requested to it here.
      if (ac.scalerating != 1.0)
	     ed[i].weight /= ac.scalerating;
  
      vertex_data & vdata = _g->vertex_data(ed[i].from - matlab_offset);
      vdata.datapoint.add_elem(ed[i].to - matlab_offset, ed[i].weight);  
      if (ps.algorithm == K_MEANS){ //compute mean for each cluster by summing assigned points
         ps.clusts.cluster_vec[vdata.current_cluster].cur_sum_of_points[ed[i].to - matlab_offset] += ed[i].weight;  
      }
      if (! vdata.reported){
         vdata.reported = true;
         if (ps.algorithm == K_MEANS) 
           ps.clusts.cluster_vec[vdata.current_cluster].num_assigned_points++;
         ps.total_assigned++; //count the total number of non-zero rows we encountered
       }
    }
      printf(".");
      fflush(0);
      if (rc == 0 || total >= edgecount_in_file)
       break;

  }
  if (total != (int)e){
      logstream(LOG_ERROR) << "Missing edges in file. Should be " << e << " but in file we counted only " << total << " edges. Please check your conversion script and verify the file is not truncated and edges are not missing. " << std::endl;
  }
  assert(total == (int)e);
  delete [] ed; ed = NULL;
  ps.L = e;
  return e;
}



#include <graphlab/macros_undef.hpp>
#endif
