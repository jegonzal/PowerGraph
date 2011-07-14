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


#ifndef _GL_NMF
#define _GL_NMF

#include "graphlab.hpp"
#include <graphlab/macros_def.hpp>

extern int D,M,N;
gl_types::core * glcore;
double * x1;
double * x2;
void last_iter();

/*
function [w,h]=nmf(v,r,verbose)
%
% Jean-Philippe Brunet
% Cancer Genomics 
% The Broad Institute
% brunet@broad.mit.edu
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
% or functionality. 
%
% NMF divergence update equations :
% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix 
% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
%
% v (n,m) : N (genes) x M (samples) original matrix 
%           Numerical data only. 
%           Must be non negative. 
%           Not all entries in a row can be 0. If so, add a small constant to the 
%           matrix, eg.v+0.01*min(min(v)),and restart.
%           
% r       : number of desired factors (rank of the factorization)
%
% verbose : prints iteration count and changes in connectivity matrix elements
%           unless verbose is 0 
%           
% Note : NMF iterations stop when connectivity matrix has not changed 
%        for 10*stopconv interations. This is experimental and can be
%        adjusted.
%
% w    : N x r NMF factor
% h    : r x M NMF factor


% test for negative values in v
if min(min(v)) < 0
error('matrix entries can not be negative');
return
end
if min(sum(v,2)) == 0
error('not all entries in a row can be zero');
return
end


[n,m]=size(v);
stopconv=40;      % stopping criterion (can be adjusted)
niter = 2000;     % maximum number of iterations (can be adjusted)

cons=zeros(m,m);
consold=cons;
inc=0;
j=0;

%
% initialize random w and h
%
w=rand(n,r);
h=rand(r,m); 


for i=1:niter

% divergence-reducing NMF iterations

x1=repmat(sum(w,1)',1,m);
h=h.*(w'*(v./(w*h)))./x1;
x2=repmat(sum(h,2)',n,1);
w=w.*((v./(w*h))*h')./x2;

% test convergence every 10 iterations

if(mod(i,10)==0)  
j=j+1;

% adjust small values to avoid undeflow
h=max(h,eps);w=max(w,eps);

% construct connectivity matrix
[y,index]=max(h,[],1);   %find largest factor
mat1=repmat(index,m,1);  % spread index down
mat2=repmat(index',1,m); % spread index right
cons=mat1==mat2;

if(sum(sum(cons~=consold))==0) % connectivity matrix has not changed
inc=inc+1;                     %accumulate count 
else
inc=0;                         % else restart count
end
if verbose                     % prints number of changing elements 
fprintf('\t%d\t%d\t%d\n',i,inc,sum(sum(cons~=consold))), 
end

if(inc>stopconv)
break,                % assume convergence is connectivity stops changing 
end 

consold=cons;

end
end
*/

using namespace std;

void nmf_init(){
  x1 = new double[D];
  x2 = new double[D];
  for (int i=0; i<D; i++){
     x1[i] = x2[i] = 0;
  }
}
const gl_types::edge_list get_edges(bool isuser, gl_types::iscope & scope){
     return isuser ? scope.out_edge_ids(): scope.in_edge_ids();
}
const vertex_data& get_neighbor(bool isuser, gl_types::iscope & scope, edge_id_t oedgeid){
     return isuser ? scope.const_neighbor_vertex_data(scope.target(oedgeid)) :  scope.const_neighbor_vertex_data(scope.source(oedgeid));
}
     
void nmf_update_function(gl_types::iscope & scope, 
      gl_types::icallback & scheduler){


  /* GET current vertex data */
  vertex_data& user = scope.vertex_data();
 
  bool isuser = ((int)scope.vertex() < M);
  int id = scope.vertex();

  
  /* print statistics */
  if (debug&& ((id == 0) || (id == M-1) || (id  == M) || (id == M+N-1))){
    printf("NMF: entering %s node  %u \n", (isuser ? "user":"movie"), id);   
    debug_print_vec((isuser ? "V " : "U") , user.pvec, D);
  }

  user.rmse = 0;

  if (user.num_edges == 0){
    return; //if this user/movie have no ratings do nothing
  }

  double * buf = new double[isuser? N:M];
  double * ret = new double[D];
  for (int i=0; i<D; i++)
     ret[i] = 0;
  for (int i=0; i< (isuser ? N:M); i++)
     buf[i] = 0;

  gl_types::edge_list outs = get_edges(isuser, scope);
  timer t; t.start(); 
    
   foreach(gl_types::edge_id oedgeid, outs) {

#ifdef GL_NO_MULT_EDGES
      const edge_data & edge = scope.const_edge_data(oedgeid);
#else
      const multiple_edges edges = scope.const_edge_data(oedgeid);
      for (int j=0; j< (int)edges.medges.size(); j++){  
         const edge_data & edge = edges.medges[j];
#endif
      double v = edge.weight;
      const vertex_data  & movie = get_neighbor(isuser, scope, oedgeid);
      float prediction;     
      float sq_err = predict(user, movie, NULL, v , prediction);
      int pos = isuser ? (scope.target(oedgeid) - M) : scope.source(oedgeid); 
      buf[pos] = v/prediction;
      user.rmse += sq_err;
#ifndef GL_NO_MULT_EDGES
        }
#endif
   }
  counter[EDGE_TRAVERSAL] += t.current_time();

  t.start();
    foreach(graphlab::edge_id_t oedgeid, outs){
      
       const vertex_data  & movie = get_neighbor(isuser,scope, oedgeid);
       for (int j=0; j<D; j++){
         int pos = isuser ? (scope.target(oedgeid) - M) : scope.source(oedgeid); 
         ret[j] += movie.pvec[j] * buf[pos];
       }
    }
  counter[SVD_MULT_A] += t.current_time();   


  double * px;
  if (isuser)
     px = x1;
  else 
     px = x2;

  for (int i=0; i<D; i++){
     assert(px[i] != 0);
     user.pvec[i] *= ret[i] / px[i];
  }
 /* print statistics */
  if (debug&& (id == 0 || (id == M-1) || (id == M) || (id == M+N-1))){
    printf("NMF: exiting %s node  %u \n", (isuser ? "user":"movie"), id);   
    debug_print_vec((isuser ? "V " : "U") , user.pvec, D);
  }

  delete[] buf;
}

void pre_user_iter(){

   for (int i=0; i<D; i++)
      x1[i] = 0;

   for (int i=M; i<M+N; i++){
    const vertex_data & data = glcore->graph().vertex_data(i);
    for (int i=0; i<D; i++){
      x1[i] += data.pvec[i];
    }
  }
}
void pre_movie_iter(){
   for (int i=0; i<D; i++)
      x2[i] = 0;

   for (int i=0; i<M; i++){
    const vertex_data & data = glcore->graph().vertex_data(i);
    for (int i=0; i<D; i++){
      x2[i] += data.pvec[i];
    }
  }
}

void nmf(gl_types::core * _glcore){

   glcore = _glcore;
  
   nmf_init();
   std::vector<vertex_id_t> rows,cols;
   for (int i=0; i< M; i++)
      rows.push_back(i);
   for (int i=M; i< M+N; i++)
      cols.push_back(i);
 

   for (int j=1; j<= svd_iter+1; j++){

     pre_user_iter();
     if (debug)
        cout<<"x1: " << x1[0] <<endl;
     glcore->add_tasks(rows, nmf_update_function, 1);
     glcore->start();
     
     pre_movie_iter();
     if (debug)
        cout<<"x2:" << x2[0] << endl;
     glcore->add_tasks(cols, nmf_update_function, 1);
     glcore->start();
     

     //last_iter();  //TODO
   }
}

#include "graphlab/macros_undef.hpp"
#endif
