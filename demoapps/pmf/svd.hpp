/**  
            violation = -Gp;   
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
double predict(const vertex_data& user, const vertex_data &movie, float rating, float & prediction);
void verify_result(double a, double b, double c);


//LANCZOS VARIABLES
extern vec lancbeta;
vec lancbeta2;
extern vec lancalpha;
vec lancalpha2;
extern int offset, offset2, offset3;



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
   lancbeta = zeros(m+3);
   lancbeta2 = zeros(m+3);
   lancalpha = zeros(m+3);
   lancalpha2 = zeros(m+3);
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

  sum = 0;
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec = zeros(m+3);
    data->pvec[1] = ac.debug ? 0.5 : rand();
    sum += data->pvec[1]*data->pvec[1];
  }

  sum = sqrt(sum);
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec[1] /= sum;
    if (ac.debug && i < 20)
      cout<<"Initial V2(:,2) is " << data->pvec[1] << endl;
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

  if (user.num_edges == 0){
    return; //if this user/movie have no ratings do nothing
  }


  gl_types::edge_list outs = scope.out_edge_ids();
  timer t;
  t.start(); 

   foreach(gl_types::edge_id oedgeid, outs) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      user.pvec[0] += edge.weight * movie.pvec[offset];
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

  if (user.num_edges == 0){
    return; //if this user/movie have no ratings do nothing
  }


  gl_types::edge_list ins = scope.in_edge_ids();
  timer t;
  t.start(); 

   foreach(gl_types::edge_id iedgeid, ins) {
      edge_data & edge = scope.edge_data(iedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.source(iedgeid));
      user.pvec[0] += edge.weight * movie.pvec[offset];
  }
 
  ps.counter[SVD_MULT_A] += t.current_time();
  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("svd_Axb: computed value  %u %g \n",  (int)scope.vertex(), user.pvec[0]);   
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
    debug_print_vec("V" , user.pvec, m);
  }

  //user.rmse = 0;
  user.pvec[0] = 0;

  //if (user.num_edges == 0){
  //  return; //if this user/movie have no ratings do nothing
 // }

  gl_types::edge_list ins = scope.in_edge_ids();
  timer t;
  t.start(); 

   foreach(gl_types::edge_id iedgeid, ins) {
      edge_data & edge = scope.edge_data(iedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.source(iedgeid));
      user.pvec[0] += edge.weight * movie.rmse;
      }

   assert(offset2 < m+2 && offset3 < m+2);
   user.pvec[0] -=  lancbeta[offset2] * user.pvec[offset3];
 
   ps.counter[SVD_MULT_A_TRANSPOSE] += t.current_time();

  if (ac.debug&& ((int)scope.vertex() == ps.M || ((int)scope.vertex() == ps.M+ps.N-1))){
    printf("svd_ATxb: computed value  %u %g beta: %g v %g \n",  (int)scope.vertex(), user.pvec[0],lancbeta[offset2],  user.pvec[offset3]);   
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
    debug_print_vec("V" , user.pvec, m);
  }

  //user.rmse = 0;
  user.pvec[0] = 0;
  //if (user.num_edges == 0){
  //  return; //if this user/movie have no ratings do nothing
 // }

  gl_types::edge_list out = scope.out_edge_ids();
  timer t;
  t.start(); 

   foreach(gl_types::edge_id oedgeid, out) {
      edge_data & edge = scope.edge_data(oedgeid);
      vertex_data  & movie = scope.neighbor_vertex_data(scope.target(oedgeid));
      user.pvec[0] += edge.weight * movie.rmse;
      }

     assert(offset2 < m+2 && offset3 < m+2);
     user.pvec[0] -= lancbeta2[offset2] * user.pvec[offset3];
 
   ps.counter[SVD_MULT_A_TRANSPOSE] += t.current_time();

  if (ac.debug&& ((int)scope.vertex() ==  0 || ((int)scope.vertex() == ps.M-1))){
    printf("svd_ATxb2: computed value  %u %g beta: %g v %g \n",  (int)scope.vertex(), user.pvec[0],lancbeta2[offset2],  user.pvec[offset3]);   
  }



}


void set_rmse(){
  const graph_type *g = ps.g<graph_type>(TRAINING);
  for (int i=0; i< ps.M+ps.N; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->rmse = data->pvec[0]; 
 }
}
 
double wTV(int j);

double wTV2(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING);

  timer t; t.start();
  double lancalpha = 0;
  for (int i=0; i< ps.M; i++){ 
    const vertex_data * data = &g->vertex_data(i);
    lancalpha+= data->rmse*data->pvec[j];
  }
  ps.counter[CALC_RMSE_Q] += t.current_time();
  if (ac.debug)
	cout<<"alpha2: " << lancalpha<<endl;

  return lancalpha;

}
double w_lancalphaV(int j);

double w_lancalphaV2(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING);
  
  timer t; t.start();
  double norm = 0;
  if (ac.debug)
	cout << "w: " ;
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->rmse -= lancalpha2[j]*data->pvec[j];
    if (ac.debug && i <20)
	cout<<data->rmse<<" ";
    norm += data->rmse*data->rmse;
  }
  ps.counter[CALC_RMSE_Q] += t.current_time();
  if (ac.debug){
	cout<<endl;
        cout<<"Beta2: " << sqrt(norm) << endl;
  }
  return sqrt(norm);
}


void update_V(int j);

void update_V2(int j){
  const graph_type *g = ps.g<graph_type>(TRAINING); 

  timer t; t.start();
  if (ac.debug)
	cout<<"V2: ";
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = (vertex_data*)&g->vertex_data(i);
    data->pvec[j] = data->rmse / lancbeta2[j];
    if (ac.debug && i <20)
	cout<<data->pvec[j]<<" ";
  }
  if (ac.debug)
	cout<<endl;
  ps.counter[CALC_RMSE_Q] += t.current_time();
}


mat calc_V();

mat calc_V2(){

  const graph_type *g = ps.g<graph_type>(TRAINING); 
  
  mat V = zeros(ps.M,ac.iter+1);
  for (int i=0; i< ps.M; i++){ 
    const vertex_data * data = (vertex_data*)&g->vertex_data(i);
    set_row(V, i, mid(data->pvec, 1, ac.iter+1));
  }
  return V;
}


void print_w(bool rows);


template<typename core>
void svd(core & glcore){
  assert(false);
}

template<>
void svd<>(gl_types::core & glcore){
   
  std::vector<vertex_id_t> rows,cols;
   for (int i=0; i< ps.M; i++)
      rows.push_back(i);
   for (int i=ps.M; i< ps.M+ps.N; i++)
      cols.push_back(i);
 
   //for j=2:m+2
   for (int j=1; j<= ac.iter+1; j++){
        //w = A*V(:,j) 
        offset = j;
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
        offset2 = j;
        offset3 = j-1;
        glcore.add_tasks(rows, svd_ATxb2, 1);
        glcore.add_tasks(cols, svd_ATxb, 1);
	glcore.start();

        set_rmse();

        if (ac.debug){
          print_w(false);
          print_w(true);
        }
        //lancalpha(j) = w'*V(:,j);
	lancalpha[j] = wTV(j);
	lancalpha2[j] = wTV2(j);

        //w =  w - lancalpha(j)*V(:,j);
        //lancbeta(j+1)=norm(w,2);
        lancbeta[j+1] = w_lancalphaV(j);
        lancbeta2[j+1] = w_lancalphaV2(j);

        //V(:,j+1) = w/lancbeta(j+1);
        update_V(j+1); 
        update_V2(j+1); 
        logstream(LOG_INFO) << "Finished iteration " << j << " in time: " << ps.gt.current_time() << std::endl;
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
 mat T=zeros(m+1,m+1);
 for (int i=1; i<=m; i++){
   set_val(T,i-1,i-1,lancalpha[i]);
   set_val(T,i-1,i,lancbeta[i+1]);
   set_val(T,i,i-1,lancbeta[i+1]);
 }
 set_val(T,m,m,lancalpha[m+1]);
 mat T2=zeros(m+1,m+1);
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

 mat Vectors=calc_V(); 
 vec eigenvalues; 
 mat eigenvectors;
 assert(::eig_sym(T, eigenvalues, eigenvectors));
 cout << "Here are the computed eigenvalues" << endl;
 for (int i=0; i< std::min((int)eigenvalues.size(),20); i++)
	cout<<"eigenvalue " << i << " val: " << eigenvalues[i] << endl;

 mat Vectors2=calc_V2(); 
 vec eigenvalues2; 
 mat eigenvectors2;
 assert(::eig_sym(T2, eigenvalues2, eigenvectors2));
 cout << "Here are the computed eigenvalues: other side" << endl;
 for (int i=0; i< std::min((int)eigenvalues2.size(),20); i++)
	cout<<"eigenvalue2 " << i << " val: " << eigenvalues2[i] << endl;


 ps.U=Vectors*eigenvectors;
 if (ac.debug && ps.U.size() < 1000)
   cout<<"Eigen vectors are:" << ps.U << endl << "V is: " << Vectors << endl << " Eigenvectors (u) are: " << eigenvectors;

 ps.V=Vectors2*eigenvectors2;
 if (ac.debug && ps.V.size() < 1000)
   cout<<"Eigen vectors2 are:" << ps.V << endl << "V is: " << Vectors2 << endl << " Eigenvectors (u) are: " << eigenvectors2;

 ps.T=zeros(eigenvalues.size(),2);
 set_col(ps.T,0,::sqrt(fabs(eigenvalues))); //should be all positive eigenvalues, but because of some accuracy error sometimes
 set_col(ps.T,1,::sqrt(fabs(eigenvalues2))); //we get small negative numbers

 if (ac.unittest > 0){
   verify_result(0, 0, 0);
 }
}

#include "graphlab/macros_undef.hpp"
#endif //__SVD_HPP
