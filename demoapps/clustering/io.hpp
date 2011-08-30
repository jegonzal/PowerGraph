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

template<typename edgedata>
int read_edges(FILE * f, int nodes, graph_type * _g);


void fill_output(){
   ps.output_clusters = CLUSTER_LOCATIONS.get_val();
   ps.output_assignements = zeros(ps.M);
   for (int i=0; i< ps.M; i++){ 
      vertex_data & data = ps.g->vertex_data(i);
      ps.output_assignements[i] = data.current_cluster;
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

 CLUSTER_LOCATIONS.set(clusters);  
}


void add_vertices(graph_type * _g){
  assert(ps.K > 0);
  vertex_data vdata;
  // add M movie nodes (ps.tensor dim 1)
  for (int i=0; i<ps.M; i++){
    vdata.current_cluster = randi(0, ps.K-1);
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
  assert(_K == 0); 
  assert(_M>=1 && _N>=1); 


  ps.M=_M; ps.N= _N;
  add_vertices(_g);
 
  // read tensor non zero edges from file
  int val = read_edges<edge_float>(f,ps.M, _g);
  assert(val == ps.L);

  fclose(f);
}


/**
 * read edges from file, with support with multiple edges between the same pair of nodes (in different times)
 */
template<typename edgedata>
int read_edges(FILE * f, int nodes, graph_type * _g){
     
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
      assert((int)ed[i].from >= matlab_offset && (int)ed[i].from <= nodes);
      assert((int)ed[i].to >= matlab_offset && (int)ed[i].to <= nodes);
      //no self edges
      assert((int)ed[i].to != (int)ed[i].from);
    
      //if sacling of rating values is requested to it here.
      if (ac.scalerating != 1.0)
	     ed[i].weight /= ac.scalerating;
  
      vertex_data & vdata = _g->vertex_data(ed[i].from);
      vdata.datapoint.set_new(ed[i].to, ed[i].weight);   
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
