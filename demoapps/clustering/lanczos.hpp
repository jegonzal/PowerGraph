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


#ifndef __LANCZOS_HPP
#define __LANCZOS_HPP

#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>

#define pvec distances
#define rmse min_distance
#define num_edges reported
#define U output_clusters
#define V output_assignements
using namespace std;
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
double predict(const vertex_data& user, const vertex_data &movie, float rating, float & prediction);
void verify_result(double a, double b, double c);
flt_dbl* find_pvec(int pos, int i, vertex_data* data);


//LANCZOS VARIABLES
flt_dbl_vec lancbeta;
flt_dbl_vec lancalpha;
int offset, offset2, offset3;



/**
 *
 * [n,k] = size(A);
   V = zeros(k,m+1);
   V(:,2) = rand(k,1);
   V(:,2)=V(:,2)/norm(V(:,2),2);
   beta(2)=0;
 *
 * */
void init_lanczos(){
   int m = ac.iter;
   lancbeta = zeros(m+3);
   lancalpha = zeros(m+3);
   double sum = 0;

  const graph_type *g = ps.g<graph_type>(TRAINING);

  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec = zeros(m+3);
    data->pvec[1] = ac.debug ? 0.5 : rand();
    sum += data->pvec[1]*data->pvec[1];
  }

  sum = sqrt(sum);
  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec[1] /= sum;
    if (ac.debug && i- ps.M < 20)
      cout<<"Initial V(:,2) is " << data->pvec[1] << endl;
  }
}

/***
 * UPDATE FUNCTION (ROWS)
 */
void Axb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
 
  
  /* print statistics */
  if (ac.debug&& (scope.vertex() == 0 || ((int)scope.vertex() == ps.M-1))){
    printf("Axb: entering  node  %u \n",  (int)scope.vertex());   
    //debug_print_vec("V" , user.pvec, D);
  }

  user.rmse = 0;

  if (!user.num_edges ){
    return; //if this user/movie have no ratings do nothing
  }


  //gl_types::edge_list outs = scope.out_edge_ids();
  timer t;
  t.start(); 

  //   foreach(gl_types::edge_id oedgeid, outs) {
  FOR_ITERATOR(i,user.datapoint){
      //edge_data & edge = scope.edge_data(oedgeid);
      double weight = get_nz_data(user.datapoint, i);
      int index = get_nz_index(user.datapoint, i);
      //vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      assert(index >= ps.M && index < ps.M+ps.N);
      vertex_data& movie = ps.g<graph_type>(TRAINING)->vertex_data(index);
      //user.rmse += edge.weight * movie.pvec[offset];
      user.rmse += weight * movie.pvec[offset];
  }

  ps.counter[SVD_MULT_A] += t.current_time();
  if (ac.debug&& (scope.vertex() == 0 || ((int)scope.vertex() == ps.M-1))){
    printf("Axb: computed value  %u %g \n",  (int)scope.vertex(), user.rmse);   
  }


}
/***
 * UPDATE FUNCTION (COLS)
 */
void ATxb(gl_types::iscope &scope, 
	 gl_types::icallback &scheduler) {
    


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
  int m = ac.iter; 
  
  /* print statistics */
  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("Axb: entering  node  %u \n",  (int)scope.vertex());   
    debug_print_vec("V" , user.pvec, m);
  }

  user.rmse = 0;

  if (!user.num_edges){
    return; //if this user/movie have no ratings do nothing
  }

  //gl_types::edge_list ins = scope.in_edge_ids();
  timer t;
  t.start(); 

   //foreach(gl_types::edge_id iedgeid, ins) {
   FOR_ITERATOR(i, user.datapoint){
      //edge_data & edge = scope.edge_data(iedgeid);
      double weight = get_nz_data(user.datapoint, i);
      int index = get_nz_index(user.datapoint, i);
     
      //vertex_data  & movie = scope.neighbor_vertex_data(scope.source(iedgeid));
      vertex_data & movie = ps.g<graph_type>(TRAINING)->vertex_data(index);
      //user.rmse += edge.weight * movie.rmse;
      user.rmse += weight * movie.rmse;
   }
   

   assert(offset2 < m+2 && offset3 < m+2);
   user.rmse -= lancbeta[offset2] * user.pvec[offset3];
   ps.counter[SVD_MULT_A_TRANSPOSE] += t.current_time();

  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("Axb: computed value  %u %g beta: %g v %g \n",  (int)scope.vertex(), user.rmse,lancbeta[offset2],  user.pvec[offset3]);   
  }



}

 
double wTV(int j){
  graph_type *g = ps.g<graph_type>(TRAINING);

  timer t; t.start();
  double lancalpha = 0;
  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = &g->vertex_data(i);
    //lancalpha+= data->rmse*data->pvec[j];
    lancalpha+= data->rmse* *find_pvec(j, i, data);
  }
  ps.counter[CALC_RMSE_Q] += t.current_time();
  if (ac.debug)
	cout<<"alpha: " << lancalpha<<endl;

  return lancalpha;

}
double w_lancalphaV(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING);
  
  timer t; t.start();
  double norm = 0;
  if (ac.debug)
	cout << "w: " ;
  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    //data->rmse -= lancalpha[j]*data->pvec[j];
    data->rmse -= lancalpha[j]* *find_pvec(j, i, data);
    if (ac.debug && i-ps.M<20)
	cout<<data->rmse<<" ";
    norm += data->rmse*data->rmse;
  }
  ps.counter[CALC_RMSE_Q] += t.current_time();
  if (ac.debug){
	cout<<endl;
        cout<<"Beta: " << sqrt(norm) << endl;
  }
  return sqrt(norm);
}

void update_V(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING); 

  timer t; t.start();
  if (ac.debug)
	cout<<"V: ";
  for (int i=ps.M; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    //data->pvec[j] = data->rmse / lancbeta[j];
    *find_pvec(j,i,data) = data->rmse/ lancbeta[j];
    if (ac.debug && i-ps.M<20)
	cout<<*find_pvec(j,i,data)<<" ";
  }
  if (ac.debug)
	cout<<endl;
  ps.counter[CALC_RMSE_Q] += t.current_time();
}


void print_w(bool rows){

  const graph_type *g = ps.g<graph_type>(TRAINING); 
  
  int size =rows? ps.M : ps.N;
  int start=rows? 0 : ps.M;
  int end = rows?ps.M : ps.M+ ps.N;
  flt_dbl_vec v = zeros(size);
  for (int i=start; i< end; i++){ 
    const vertex_data * data = (vertex_data*)&g->vertex_data(i);
    v[i - start] = data->rmse;
  }
  cout<<"w is: " << mid(v,0,20) << endl;
  if (end - start > 40)
    cout<<"w end is: " << mid(v, v.size()-20, 20) << endl;
}
mat calc_V();

template<typename core>
void lanczos(core & glcore){
  assert(false);
}

template<>
void lanczos<>(gl_types::core & glcore){
   
  std::vector<vertex_id_t> rows,cols;
   for (int i=0; i< ps.M; i++)
      rows.push_back(i);
   for (int i=ps.M; i< ps.M+ps.N; i++)
      cols.push_back(i);
 
   ac.iter--;

   //for j=2:m+2
   for (int j=1; j<= ac.iter+1; j++){
        //w = A*V(:,j) 
        offset = j;
	glcore.add_tasks(rows, Axb, 1);
        glcore.start();
	if (ac.debug)
           print_w(true);
        //w = w - lancbeta(j)*V(:,j-1);
        offset2 = j;
        offset3 = j-1;
        glcore.add_tasks(cols, ATxb, 1);
	glcore.start();
        if (ac.debug)
          print_w(false);
       
        //lancalpha(j) = w'*V(:,j);
	lancalpha[j] = wTV(j);

        //w =  w - lancalpha(j)*V(:,j);
        //lancbeta(j+1)=norm(w,2);
        lancbeta[j+1] = w_lancalphaV(j);

        //V(:,j+1) = w/lancbeta(j+1);
        update_V(j+1); 
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
 int m = ac.iter;
 mat T=fmat2mat(zeros(m+1,m+1));
 for (int i=1; i<=m; i++){
   set_val(T,i-1,i-1,lancalpha[i]);
   set_val(T,i-1,i,lancbeta[i+1]);
   set_val(T,i,i-1,lancbeta[i+1]);
 }
 set_val(T,m,m,lancalpha[m+1]);
 if (ac.debug && m < 100){
    cout<<"Matrix T is: " << T << endl;
 }

 mat Vectors=calc_V(); 
   
 vec eigenvalues; 
 mat eigenvectors;
 assert(::eig_sym(T, eigenvalues, eigenvectors));
 cout << "Here are the computed eigenvalues, from larger to smaller" << endl;
 for (int i=0; i< std::min((int)eigenvalues.size(),20); i++)
	cout<<"eigenvalue " << i << " val: " << eigenvalues[i] << endl;


 ps.U=mat2fmat(Vectors*eigenvectors);
 if (ac.debug)
   cout<<"Eigen vectors are:" << ps.U << endl << "V is: " << Vectors << endl << " Eigenvectors (u) are: " << eigenvectors;
 ps.V=zeros(eigenvalues.size(),1);
 set_col(ps.V,0,vec2fvec(eigenvalues)); 

 if (ac.unittest > 0){
   verify_result(0, 0, 0);
 }
}

#include "graphlab/macros_undef.hpp"
#endif //__LANCZOS_HPP
