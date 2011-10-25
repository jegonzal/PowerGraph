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


#ifndef __SVD_HPP
#define __SVD_HPP

#include <stdio.h>
#include <stdlib.h>
#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>


/**
 *
 *  Implementation of the Lanczos algorithm, as given in:
 *  http://en.wikipedia.org/wiki/Lanczos_algorithm
 * 
 *  Code written by Danny Bickson, CMU, June 2011
 * */

extern advanced_config ac;
extern problem_setup ps;


using namespace graphlab;

void last_iter();
double predict(const vertex_data& user, const vertex_data &movie, flt_dbl rating, flt_dbl & prediction);
void verify_result(double a, double b, double c);

void debug_print_vec2(const char * name, flt_dbl_vec & v, int size, int i){
   if (!ac.reduce_mem_consumption)
      return debug_print_vec(name, v, size);
   else printf("%s: %g %g\n", name, *find_pvec(ps.iiter-1,i,NULL),*find_pvec(ps.iiter,i,NULL)); 
}

//LANCZOS VARIABLES
extern flt_dbl_vec lancbeta;
flt_dbl_vec lancbeta2;
extern flt_dbl_vec lancalpha;
flt_dbl_vec lancalpha2;
extern int offset, offset2, offset3;

struct global_pvec{
   float** pvec;
   int size;
   bool swap;
   global_pvec(int _size){
     size = _size;
     swap = false;
     pvec = new float*[2];
     for (int i=0; i<2; i++){
      pvec[i] = new float[size];
      memset(pvec[i], 0, size*sizeof(float));
     }
   }
};
global_pvec * pglobal_pvec;

void set_pvec(flt_dbl* pos, double val){
  *pos = val;
}
flt_dbl* find_pvec(int pos, int i, vertex_data* data){
  if (!ac.reduce_mem_consumption || (pos == 0 && data!= NULL))
    return &data->pvec[pos];
  else {
    int offset = 0;
    if (pglobal_pvec->swap)
      offset = 1;
      assert(pos - ps.iiter - offset == -1 || pos - ps.iiter -offset == 0);
    assert(i >= 0 && i < pglobal_pvec->size);
    return &pglobal_pvec->pvec[1+(pos-ps.iiter-offset)][i];
  }
}
void swap_global_pvec(){
   FILE * pfile = open_file(((ac.datafile + "swap") + boost::lexical_cast<std::string>(ps.iiter)).c_str(), "wb");
   write_vec(pfile, pglobal_pvec->size, pglobal_pvec->pvec[0]);
   fclose(pfile); 
   memcpy(pglobal_pvec->pvec[0], pglobal_pvec->pvec[1], pglobal_pvec->size*sizeof(flt_dbl));
   memset(pglobal_pvec->pvec[1], 0, pglobal_pvec->size*sizeof(flt_dbl));
   pglobal_pvec->swap = true;
}

void save_vec(const char * name, flt_dbl_vec & v){
  FILE * pfile = open_file(name,"wb");
  write_vec(pfile, v.size(), data(v));
  fclose(pfile);
}
void load_vec(const char * name, flt_dbl_vec & v, int size){
  FILE * pfile = open_file(name, "r");
  v = zeros(size);
  read_vec(pfile, size, (flt_dbl*)data(v));
  fclose(pfile);
}

void read_lanc_alpha_beta(int size){
  assert(ac.svd_compile_eigenvectors && ac.reduce_mem_consumption);
  logstream(LOG_INFO) << "Loading alpha beta from file" << std::endl;
  load_vec((ac.datafile + "lancalpha").c_str(), lancalpha, size);
  load_vec((ac.datafile + "lancalpha2").c_str(), lancalpha2, size);
  load_vec((ac.datafile + "lancbeta").c_str(), lancbeta, size);
  load_vec((ac.datafile + "lancbeta2").c_str(), lancbeta2, size);
}


void write_lanc_alpha_beta(){
   assert(ac.reduce_mem_consumption);
   logstream(LOG_INFO) << "Saving alpha beta to file" << std::endl;
   save_vec((ac.datafile + "lancalpha").c_str(), lancalpha);
   save_vec((ac.datafile + "lancalpha2").c_str(),lancalpha2);
   save_vec((ac.datafile + "lancbeta").c_str(), lancbeta);
   save_vec((ac.datafile + "lancbeta2").c_str(),lancbeta2);
}

void end_iter(){
   assert(ac.reduce_mem_consumption);
   pglobal_pvec->swap = false;
   write_lanc_alpha_beta();
}
/**
 *
 * [n,k] = size(A);
   V = zeros(k,m+1);
   V(:,2) = rand(k,1);
   V(:,2)=V(:,2)/norm(V(:,2),2);
   beta(2)=0;
 *
 * */
void init_svd(){
   int m = ac.iter;
   ps.iiter = 1;
   lancbeta = zeros(m+3);
   lancbeta2 = zeros(m+3);
   lancalpha = zeros(m+3);
   lancalpha2 = zeros(m+3);
   flt_dbl sum = 0;

  if (ac.reduce_mem_consumption)
    pglobal_pvec = new global_pvec(ps.M+ps.N);

  const graph_type *g = ps.g<graph_type>(TRAINING);

  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    if (!ac.reduce_mem_consumption)
      data->pvec = zeros(m+3);
    else data->pvec = zeros(1); 
    set_pvec( find_pvec(1, i, data), ac.debug? 0.5: rand());
    //data->pvec[1] = ac.debug ? 0.5 : rand();
    
    //sum += data->pvec[1]*data->pvec[1];
    sum += pow(*find_pvec(1, i, data), 2);
  }

  sum = sqrt(sum);
  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    set_pvec( find_pvec(1, i, data), *find_pvec(1,i,data)/ sum);
    //data->pvec[1] /= sum;
    if (ac.debug && i- ps.M < 20)
      cout<<"Initial V(:,2) is " << *find_pvec(1,i,data) << endl;
  }

  sum = 0;
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    if (!ac.reduce_mem_consumption)
      data->pvec = zeros(m+3);
    else data->pvec = zeros(1);
    //data->pvec[1] = ac.debug ? 0.5 : rand();
    set_pvec( find_pvec(1, i, data), ac.debug? 0.5: rand());
    //sum += data->pvec[1]*data->pvec[1];
    sum += pow(*find_pvec(1, i, data),2);
  }

  sum = sqrt(sum);
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    //data->pvec[1] /= sum;
    set_pvec( find_pvec(1, i, data), *find_pvec(1,i,data)/ sum);
    if (ac.debug && i < 20)
      cout<<"Initial V2(:,2) is " << *find_pvec(1,i,data) << endl;
  }
}

/***
 * UPDATE FUNCTION (ROWS)
 */
void svd_Axb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  
  /* print statistics */
  if (ac.debug&& ((int)scope.vertex() == 0 || ((int)scope.vertex() == ps.M-1))){
    printf("svd_Axb: entering  node  %u \n",  (int)scope.vertex());   
  }

  //user.rmse = 0;
  user.pvec[0] = 0;

  if (!user.num_edges){
    return; //if this user/movie have no ratings do nothing
  }


  //gl_types::edge_list outs = scope.out_edge_ids();
  timer t;
  t.start(); 

  // foreach(gl_types::edge_id oedgeid, outs) {
  FOR_ITERATOR(i,user.datapoint){
      //edge_data & edge = scope.edge_data(oedgeid);
      double weight = get_nz_data(user.datapoint, i);
      int index = get_nz_index(user.datapoint, i);
      assert(index>= 0 && index < ps.N);
      //vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      vertex_data & movie = ps.g<graph_type>(TRAINING)->vertex_data(index+ps.M);
      //user.pvec[0] += weight * movie.pvec[offset];
      user.pvec[0] += weight * *find_pvec( offset, index+ps.M, &movie);
  }
 
  ps.counter[SVD_MULT_A] += t.current_time();
  if (ac.debug&& ((int)scope.vertex() == 0 || ((int)scope.vertex() == ps.M-1))){
    printf("svd_Axb: computed value  %u %g \n",  (int)scope.vertex(), user.pvec[0]);   
  }


}
/***
 * UPDATE FUNCTION2 (ROWS)
 */
void svd_Axb2(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  
  /* print statistics */
  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("svd_Axb2: entering  node  %u \n",  (int)scope.vertex());   
  }

  //user.rmse = 0;
  user.pvec[0]= 0;

  if (!user.num_edges){
    return; //if this user/movie have no ratings do nothing
  }


  //gl_types::edge_list ins = scope.in_edge_ids();
  timer t;
  t.start(); 

  FOR_ITERATOR(i, user.datapoint){
      double weight = get_nz_data(user.datapoint, i);
      int index = get_nz_index(user.datapoint, i);
      assert(index>=0 && index < ps.M);
      vertex_data & movie = ps.g<graph_type>(TRAINING)->vertex_data(index);
      //user.pvec[0] += weight * movie.pvec[offset];
      user.pvec[0] += weight * *find_pvec( offset, index, &movie);
  }
 
  ps.counter[SVD_MULT_A] += t.current_time();
  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("svd_Axb2: computed value  %u %g \n",  (int)scope.vertex(), user.pvec[0]);   
  }


}




/***
 * UPDATE FUNCTION (COLS)
 */
void svd_ATxb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  int m = ac.iter; 
  
  /* print statistics */
  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("svd_ATxb: entering  node  %u \n",  (int)scope.vertex());   
    debug_print_vec2("V" , user.pvec, m, scope.vertex());
  }

  //user.rmse = 0;
  user.pvec[0] = 0;

  //if (user.num_edges == 0){
  //  return; //if this user/movie have no ratings do nothing
 // }

  gl_types::edge_list ins = scope.in_edge_ids();
  timer t;
  t.start(); 

   //foreach(gl_types::edge_id iedgeid, ins) {
   //   edge_data & edge = scope.edge_data(iedgeid);
   //   vertex_data  & movie = scope.neighbor_vertex_data(scope.source(iedgeid));
  FOR_ITERATOR(i, user.datapoint){
      double weight = get_nz_data(user.datapoint, i);
      int index = get_nz_index(user.datapoint, i);
      assert(index >= 0 && index < ps.M);
      vertex_data & movie = ps.g<graph_type>(TRAINING)->vertex_data(index);
      user.pvec[0] += weight * movie.rmse;
      }

   assert(offset2 < m+2 && offset3 < m+2);
   //user.pvec[0] -=  lancbeta[offset2] * user.pvec[offset3];
   user.pvec[0] -= lancbeta[offset2] * *find_pvec( offset3, scope.vertex(), &user);
 
   ps.counter[SVD_MULT_A_TRANSPOSE] += t.current_time();

  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("svd_ATxb: computed value  %u %g beta: %g v %g \n",  (int)scope.vertex(), user.pvec[0],lancbeta[offset2],  *find_pvec(offset3, scope.vertex(), &user));   
  }



}


/***
 * UPDATE FUNCTION2 (COLS)
 */
void svd_ATxb2(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  int m = ac.iter; 
  
  /* print statistics */
  if (ac.debug&& ((int)scope.vertex() == 0 || ((int)scope.vertex() == ps.M-1))){
    printf("svd_ATxb: entering  node  %u \n",  (int)scope.vertex());   
    debug_print_vec2("V" , user.pvec, m, scope.vertex());
  }

  //user.rmse = 0;
  user.pvec[0] = 0;
  //if (user.num_edges == 0){
  //  return; //if this user/movie have no ratings do nothing
 // }

  gl_types::edge_list out = scope.out_edge_ids();
  timer t;
  t.start(); 

  FOR_ITERATOR(i, user.datapoint){
      double weight = get_nz_data(user.datapoint, i);
      int index = get_nz_index(user.datapoint, i);
      assert(index >=0 && index< ps.N);
      vertex_data & movie = ps.g<graph_type>(TRAINING)->vertex_data(index+ps.M);
  //  foreach(gl_types::edge_id oedgeid, out) {
  //    edge_data & edge = scope.edge_data(oedgeid);
  //    vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      user.pvec[0] += weight * movie.rmse;
      }

     assert(offset2 < m+2 && offset3 < m+2);
     //user.pvec[0] -= lancbeta2[offset2] * user.pvec[offset3];
     user.pvec[0] -= lancbeta2[offset2] * *find_pvec(offset3, scope.vertex(), &user);
 
   ps.counter[SVD_MULT_A_TRANSPOSE] += t.current_time();

  if (ac.debug&& ((int)scope.vertex() ==  0 || ((int)scope.vertex() == ps.M-1))){
    printf("svd_ATxb2: computed value  %u %g beta: %g v %g \n",  (int)scope.vertex(), user.pvec[0],lancbeta2[offset2],  *find_pvec(offset3, scope.vertex(), &user));   
  }



}


void set_rmse(){
  const graph_type *g = ps.g<graph_type>(TRAINING);
#pragma omp parallel for
  for (int i=0; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->rmse = data->pvec[0]; 
 }
}
 
double wTV(int j);

double wTV2(int j){
  graph_type *g = ps.g<graph_type>(TRAINING);

  double lancalpha = 0;
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = &g->vertex_data(i);
    //lancalpha+= data->rmse*data->pvec[j];
    lancalpha+= data->rmse * *find_pvec(j, i, data);
  }
  if (ac.debug)
	cout<<"alpha2: " << lancalpha<<endl;

  return lancalpha;

}
double w_lancalphaV(int j);

double w_lancalphaV2(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING);
  
  double norm = 0;
  if (ac.debug)
	cout << "w: " ;
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->rmse -= lancalpha2[j]* *find_pvec(j, i, data);
    if (ac.debug && i <20)
	cout<<data->rmse<<" ";
    norm += data->rmse*data->rmse;
  }
  if (ac.debug){
	cout<<endl;
        cout<<"Beta2: " << sqrt(norm) << endl;
  }
  return sqrt(norm);
}


void update_V(int j);

void update_V2(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING); 

  if (ac.debug)
	cout<<"V2: ";

#pragma omp parallel for
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    *find_pvec(j, i, data) = data->rmse/ lancbeta2[j];
    if (ac.debug && i <20)
        cout << *find_pvec(j, i, data) << " ";
  }
  if (ac.debug)
	cout<<endl;
}

flt_dbl_mat calc_V(bool other_side, const flt_dbl_mat & eigenvectors){ 
 
  int start = ps.M;
  int end = ps.M+ps.N;
  if (other_side){
    start = 0;
    end = ps.M;
  }

  
 if (!ac.reduce_mem_consumption){ 
     flt_dbl_mat V = zeros(end-start,ac.iter+1);
     const graph_type *g = ps.g<graph_type>(TRAINING); 
     for (int i=start; i< end; i++){ 
       const vertex_data * data = (vertex_data*)&g->vertex_data(i);
      set_row(V, i-start, mid(data->pvec, 1, ac.iter+1));
     }
     return V*eigenvectors;
  }
  else {
    int block_size = std::min(end-start, ac.svd_compile_eigenvectors_block_size);
    int howmany = (end-start)/block_size;
    int reminder = (end-start)%block_size;
    if (reminder > 0)
       howmany++;

    for (int cnt=0; cnt < howmany; cnt++){
      logstream(LOG_INFO) << "Processing block number " << cnt << " at time " << ps.gt.current_time() << std::endl;
      flt_dbl_mat V = zeros(((cnt==howmany-1 && reminder>0) ? reminder : block_size), ac.iter+1);
      for (int i=1; i<= ac.iter+1; i++){
        FILE * pfile = open_file(((ac.datafile + "swap") + boost::lexical_cast<std::string>(i+1)).c_str(), "r");
        read_vec(pfile, ps.M+ps.N, pglobal_pvec->pvec[0]);
        fclose(pfile);
        assert(start+cnt*block_size < end);
        flt_dbl_vec col(pglobal_pvec->pvec[0] +start+cnt*block_size, ((cnt==howmany-1 && reminder>0) ? reminder : block_size));
        set_col(V, i-1, col);
      }
      it_file output((ac.datafile + (other_side ? ".U" : ".V") + boost::lexical_cast<std::string>(cnt)).c_str());
      output << Name(other_side ? "U" : "V");
      output << fmat2mat(V*eigenvectors);
      output.close();
 
      if (cnt == howmany-1)
        return V*eigenvectors;
   }
 }
 return zeros(1,1);
}
/* 
mat calc_V(bool other_side){

  int start = ps.M;
  int end = ps.M+ps.N;
  if (other_side){
    start = 0;
    end = ps.M;
  }

  if (ac.debug)
    logstream(LOG_INFO) << "Allocating a matrix of size: " << ((end-start)*ac.iter+1) <<  " time: " << ps.gt.current_time() << std::endl;
  mat V = dbl_fzeros(end-start,ac.iter+1);
  if (ac.debug)
    logstream(LOG_INFO) << "Done! in time" << ps.gt.current_time() << std::endl; 

    const graph_type *g = ps.g<graph_type>(TRAINING); 
    for (int i=start; i< end; i++){ 
      const vertex_data * data = (vertex_data*)&g->vertex_data(i);
      set_row(V, i-start, fvec2vec(mid(data->pvec, 1, ac.iter+1)));
    }
  return V;
}*/


/*
void calc_V_block(mat &V, int start, int end, bool other_side){
  assert(start <= end);
  int offset = ps.M;
  int total = ps.N;
  if (other_side){
    offset = 0;
    total = ps.M;
  }


  if (ac.debug)
    logstream(LOG_INFO) << "Allocating a matrix of size: " << (total*(end-start)) <<  " time: " << ps.gt.current_time() << std::endl;
  V = dbl_fzeros(total,end-start);

  if (ac.debug)
    logstream(LOG_INFO) << "Done! in time" << ps.gt.current_time() << std::endl; 

  int size = ac.iter/10;
  for (int i=start; i< end; i++){
      FILE * pfile = open_file(((ac.datafile + "swap") + boost::lexical_cast<std::string>(i+1)).c_str(), "r");
      read_vec(pfile, ps.M+ps.N, pglobal_pvec->pvec[0]);
      fclose(pfile);
      flt_dbl_vec col(pglobal_pvec->pvec[0] + offset, total);
      assert(i - start >= 0 && i - start < V.cols());
      set_col(V, i-start, fvec2vec(col));
      if (ac.debug && (i%size) == 0)
        logstream(LOG_INFO) << "Read total of " << i << " columns at time " << ps.gt.current_time() << std::endl;
  }
}*/

/*
mat calc_V2(){
    
  mat V = fmat2mat(zeros(ps.M,ac.iter+1));

 const graph_type *g = ps.g<graph_type>(TRAINING); 
 for (int i=0; i< ps.M; i++){ 
   const vertex_data * data = (vertex_data*)&g->vertex_data(i);
   set_row(V, i, fvec2vec(mid(data->pvec, 1, ac.iter+1)));
  }
  return V;
}*/
vec calc_eigenvalues(mat & T, bool other_side){
 vec eigenvalues; 
 mat eigenvectors;
 assert(::eig_sym(T, eigenvalues, eigenvectors));
 cout << "Here are the computed eigenvalues" << endl;
 for (int i=0; i< std::min((int)eigenvalues.size(),20); i++)
	cout<<"eigenvalue " << i << " val: " << eigenvalues[i] << endl;
 
 flt_dbl_mat V = calc_V(other_side,mat2fmat(eigenvectors));    

 if (!ac.reduce_mem_consumption)
     (other_side ? ps.output_clusters : ps.output_assignements) = V;

 if (ac.debug && V.size() < 1000)
     cout<<"V is: " << V << endl;
 return eigenvalues;
}


void print_w(bool rows);

void extract_eigenvectors(){
 int m = ac.iter;
 mat T=fmat2mat(zeros(m+1,m+1));

 if (ac.svd_compile_eigenvectors && ac.reduce_mem_consumption){
    read_lanc_alpha_beta(m+3);
    pglobal_pvec = new global_pvec(ps.M+ps.N);
 }

 if (ac.debug){
   debug_print_vec("lancalpha", lancalpha, 20);
   debug_print_vec("lancbeta", lancbeta, 20);
 } 

for (int i=1; i<=m; i++){
   set_val(T,i-1,i-1,lancalpha[i]);
   set_val(T,i-1,i,lancbeta[i+1]);
   set_val(T,i,i-1,lancbeta[i+1]);
 }
 set_val(T,m,m,lancalpha[m+1]);

 vec eigenvalues = calc_eigenvalues(T, false);

 ps.T=zeros(T.rows(),2);
 set_col(ps.T,0,vec2fvec(eigenvalues)); 

 mat T2=fmat2mat(zeros(m+1,m+1));
 for (int i=1; i<=m; i++){
   set_val(T2,i-1,i-1,lancalpha2[i]);
   set_val(T2,i-1,i,lancbeta2[i+1]);
   set_val(T2,i,i-1,lancbeta2[i+1]);
 }
 set_val(T2,m,m,lancalpha2[m+1]);
  if (ac.debug && m < 100){
    cout<<"Matrix T is: " << T << " size of T: " << T.rows() << ":" << T.cols() << endl;
    cout<<"Matrix T2 is: " << T2 << endl;
 }

 vec eigenvalues2 = calc_eigenvalues(T2, true);

 set_col(ps.T,1,vec2fvec(eigenvalues2)); 
 if (ac.svd_compile_eigenvectors && ac.reduce_mem_consumption){
    save_matrix((ac.datafile + ".D").c_str(), "D", fmat2mat(ps.T));
 }
}

template<typename core>
void svd(core & glcore){
  assert(false);
}

template<>
void svd<>(gl_types::core & glcore){
  
   init_svd();

   std::vector<vertex_id_t> rows,cols;
   for (int i=0; i< ps.M; i++)
      rows.push_back(i);
   for (int i=ps.M; i< ps.M+ps.N; i++)
      cols.push_back(i);
 
   //for j=2:m+2
   for (ps.iiter=1; ps.iiter<= ac.iter+1; ps.iiter++){
        //w = A*V(:,j) 
        offset = ps.iiter;
        offset2 = offset3 = -1;
	glcore.add_tasks(rows, svd_Axb, 1);
        glcore.add_tasks(cols, svd_Axb2, 1);
        glcore.start();
        set_rmse();
	if (ac.debug){
           print_w(true);
           print_w(false);
        }
        //w = w - lancbeta(j)*V(:,j-1);
        offset2 = ps.iiter;
        offset3 = ps.iiter-1;
        glcore.add_tasks(rows, svd_ATxb2, 1);
        glcore.add_tasks(cols, svd_ATxb, 1);
	glcore.start();

        set_rmse();

        if (ac.reduce_mem_consumption)
          swap_global_pvec();

        if (ac.debug){
          print_w(false);
          print_w(true);
          logstream(LOG_INFO) <<"Middle iteration " << ps.iiter << " in time: " << ps.gt.current_time() << std::endl;
        }
        //lancalpha(j) = w'*V(:,j);
	lancalpha[ps.iiter] = wTV(ps.iiter);
	lancalpha2[ps.iiter] = wTV2(ps.iiter);

        //w =  w - lancalpha(j)*V(:,j);
        //lancbeta(j+1)=norm(w,2);
        lancbeta[ps.iiter+1] = w_lancalphaV(ps.iiter);
        lancbeta2[ps.iiter+1] = w_lancalphaV2(ps.iiter);

        //V(:,j+1) = w/lancbeta(j+1);
        update_V(ps.iiter+1); 
        update_V2(ps.iiter+1); 
        logstream(LOG_INFO) << "Finished iteration " << ps.iiter << " in time: " << ps.gt.current_time() << std::endl;

        if (ac.reduce_mem_consumption)
          end_iter();
   } 
  /* 
 * T=sparse(m+1,m+1);
 * for i=2:m+1
 *     T(i-1,i-1)=lancalpha(i);
 *     T(i-1,i)=lancbeta(i+1);
 *     T(i,i-1)=lancbeta(i+1);
 * end 
 * T(m+1,m+1)=lancalpha(m+2);
 * V = V(:,2:end-1);
 */
 if (ac.reduce_mem_consumption)
   swap_global_pvec();

 if (ac.svd_finalize){ 
   extract_eigenvectors();
 }

 if (ac.unittest > 0){
   verify_result(0, 0, 0);
 }

 
}

#include "graphlab/macros_undef.hpp"
#endif //__SVD_HPP
