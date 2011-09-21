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


#ifndef __BPTF_H
#define __BPTF_H

#include "pmf.h"
#include "../gabp/advanced_config.h"
#include "io.hpp"

//external global variables defined at pmf.cpp
extern advanced_config ac;
extern problem_setup ps;

/* variables for BPTF */
double nuAlpha = 1;
double Walpha = 1;
double nu0 = ac.D;
double alpha = 0;
double beta = 1;
vec beta0 = init_vec("1", 1);
vec mu0T = init_vec("1", 1);
mat W0;
mat W0T;
double iWalpha;
mat iW0;
mat iW0T;
mat A_U, A_V, A_T;
vec mu_U, mu_V, mu_T;

using namespace graphlab;

mat GenDiffMat(int K){
    mat ret = zeros(K,K); 
    for (int i=0; i<K; i++){
        ret(i,i) = 2;
    }
    for (int i=1; i<K; i++){
	ret(i-1,i) = -1;
    }
    for (int i=1; i<K; i++){
	ret(i,i-1) = -1;
    }
   return ret;
}


/**
* calc the A'A part of the least squares solution inv(A'A)*A'
* for time nodes
*/
template<typename graph_type>
mat calc_MMT(int start_pos, int end_pos, vec &Umean){

  int batchSize = 1000;
  mat U = zeros(batchSize,ac.D);
  mat MMT = zeros(ac.D,ac.D);
  int cnt = 0;
  timer t;

  for (int i=start_pos; i< end_pos; i++){
    if ((i-start_pos) % batchSize == 0){
      U=zeros(batchSize, ac.D);
      cnt = 1;
    }

    const vertex_data * data= &ps.g<graph_type>(TRAINING)->vertex_data(i);
     
    vec mean = data->pvec;
    Umean += mean;
    //Q.set_col(k, ret);
    t.start(); 
    for (int s=0; s<ac.D; s++)
      U(i%batchSize,s)=mean(s);
    if (ac.debug && (i==start_pos || i == end_pos-1))
      std::cout<<" clmn "<<i<< " vec: " << mean <<std::endl;

    if ((cnt  == batchSize) || (cnt < batchSize && i == end_pos-1)){
      MMT = MMT+transpose(U)*U;
    }
    ps.counter[8] += t.current_time();
    cnt++;
  }
  Umean /= (end_pos-start_pos);
  if (ac.debug)
    cout<<"mean: "<<Umean<<endl;

  assert(MMT.rows() == ac.D && MMT.cols() == ac.D);
  assert(Umean.size() == ac.D);
  return MMT;
}
template<>
mat calc_MMT<graph_type_svdpp>(int start_pos, int end_pos, vec &Umean){
  assert(false);
}
template<>
mat calc_MMT<graph_type>(int start_pos, int end_pos, vec &Umean){
  assert(false);
}



void init_self_pot(){
  //assert(BPTF);

  W0 = eye(ac.D);
  W0T = eye(ac.D);
  iWalpha = 1.0/Walpha;
  iW0 = inv(W0);
  iW0T = inv(W0T);
  nu0 = ac.D;

  A_U = eye(ac.D); //cov prior for users
  A_V = eye(ac.D); //cov prior for movies
  A_T = eye(ac.D); //cov prior for time nodes

  mu_U = zeros(ac.D); mu_V = zeros(ac.D); mu_T = zeros(ac.D);
  printf("nuAlpha=%g, Walpha=%g, mu0=%d, muT=%g, nu=%g, "
         "beta=%g, W=%g, WT=%g bptf_burn_in=%d\n", nuAlpha, Walpha, 0, 
         mu0T[0], nu0, beta0[0], W0(1,1), W0T(1,1), ac.bptf_burn_in);


  //test_randn(); 
  //test_wishrnd();
  //test_wishrnd2(); 
  //test_chi2rnd();
  //test_wishrnd3();
  //test_mvnrndex();
}

/**
 * sample the noise level (PMF/BPTF only)
 * Euqation A.2 in Xiong paper
 */
void sample_alpha(double res2){
  
  if (ac.debug)
  printf("res is %g\n", res2); 
  
  double res = res2;
  assert(ps.BPTF);
  if (nuAlpha > 0){
    double nuAlpha_ =nuAlpha+ ps.L;
    mat iWalpha_(1,1);
    set_val(iWalpha_, 0,0,iWalpha + res);
    mat iiWalpha_ = zeros(1,1);
    iiWalpha_ = inv(iWalpha_);
    alpha = get_val(wishrnd(iiWalpha_, nuAlpha_),0,0);
    assert(alpha != 0);

    if (ac.debug)
      cout<<"Sampling from alpha" <<nuAlpha_<<" "<<iWalpha<<" "<< iiWalpha_<<" "<<alpha<<endl;
    printf("sampled alpha is %g\n", alpha); 
  }
}



// sample movie nodes hyperprior
// according to equation A.3 in Xiong paper.
template<typename graph_type>
void sample_U(){
  assert(ps.BPTF);

  vec Umean;
  mat UUT = calc_MMT<graph_type>(0,ps.M,Umean);
  
  double beta0_ = beta0[0] + ps.M;
  vec mu0_ = (ps.M*Umean)/beta0_;
  double nu0_ = nu0 +ps.M;
  vec dMu = - Umean;
  if (ac.debug)
    cout<<"dMu:"<<dMu<<"beta0: "<<beta0[0]<<" beta0_ "<<beta0_<<" nu0_ " <<nu0_<<" mu0_ " << mu0_<<endl;
  mat UmeanT = ps.M*(outer_product(Umean, Umean));
  assert(UmeanT.rows() == ac.D && UmeanT.cols() == ac.D);
  mat dMuT = (beta0[0]*ps.M/beta0_)*(outer_product(dMu, dMu));
  mat iW0_ = iW0 + UUT - UmeanT + dMuT;
  mat W0_; 
  bool ret =inv(iW0_, W0_);
  assert(ret);
  mat tmp = (W0_+transpose(W0_))*0.5;
  if (ac.debug)
    cout<<iW0<<UUT<<UmeanT<<dMuT<<W0_<<tmp<<nu0_<<endl;
  A_U = wishrnd(tmp, nu0_);
  mat tmp2;  
  ret =  inv(beta0_ * A_U, tmp2);
  assert(ret);
  mu_U = mvnrndex(mu0_, tmp2, ac.D);
  if (ac.debug)
    cout<<"Sampling from U" <<A_U<<" "<<mu_U<<" "<<Umean<<" "<<W0_<<tmp<<endl;
}

// sample user nodes hyperprior
// according to equation A.4 in Xiong paper
template<typename graph_type>
void sample_V(){

  assert(ps.BPTF);
  vec Vmean;
  mat VVT = calc_MMT<graph_type>(ps.M, ps.M+ps.N, Vmean);   

  double beta0_ = beta0[0] + ps.N;
  vec mu0_ = (ps.N*Vmean)/beta0_;
  double nu0_ = nu0 +ps.N;
  vec dMu = - Vmean;
  if (ac.debug)
    cout<<"dMu:"<<dMu<<"beta0: "<<beta0[0]<<" beta0_ "<<beta0_<<" nu0_ " <<nu0_<<endl;
  mat VmeanT = ps.N*(outer_product(Vmean, Vmean));
  assert(VmeanT.rows() == ac.D && VmeanT.cols() == ac.D);
  mat dMuT =  (beta0[0]*ps.N/beta0_)*outer_product(dMu, dMu);
  mat iW0_ = iW0 + VVT - VmeanT + dMuT;
  mat W0_;
  bool ret = inv(iW0_, W0_);
  assert(ret);
  mat tmp = (W0_+transpose(W0_))*0.5;
  if (ac.debug)
    cout<<"iW0: "<<iW0<<" VVT: "<<VVT<<" VmeanT: "<<VmeanT<<" dMuT: " <<dMuT<<"W0_"<< W0_<<" tmp: " << tmp<<" nu0_: "<<nu0_<<endl;
  A_V = wishrnd(tmp, nu0_);
  mat tmp2; 
  ret = inv(beta0_*A_V, tmp2);
  assert(ret);
  mu_V = mvnrndex(mu0_, tmp2, ac.D);
  if (ac.debug)
    cout<<"Sampling from V: A_V" <<A_V<<" mu_V: "<<mu_V<<" Vmean: "<<Vmean<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
}


//calc laplacian matrix for a chain connection between time nodes (for tensor)
//according to equation A.5 in Xiong paper
mat calc_DT(){

  assert(ps.tensor);

  mat T = zeros(ac.D, ps.K);
  for (int i=0; i<ps.K; i++){
    set_col(T,i,ps.times[i].pvec);
  }
  
  mat diff = zeros(ac.D,ps.K-1);
  for (int i=0; i<ps.K-1; i++){
    set_col(diff, i , get_col(T,i) - get_col(T,i+1));
  }
  if (ac.debug)
    cout<<"T:"<<T<<" diff: " << diff<<endl;
  
  return diff;

}

// sample from time nodes
void sample_T(){
  assert(ps.BPTF);
  assert(ps.tensor);

  double beta0_ = beta0[0] + 1;
  vec pvec = ps.times[0].pvec; 
  vec mu0_ = (pvec + beta0*mu0T[0])/beta0_;
  double nu0_ = nu0 +ps.K;
  //vec dMu = mu0 - Umean;
  if (ac.debug){
    cout<<"beta0_ " << beta0_ << " beta0: " << beta0[0] << " nu0_ " << nu0_ << endl;
  } 

  mat dT = calc_DT();
  vec dTe = pvec - mu0T[0];
  mat iW0_ = iW0T + dT*transpose(dT) + (beta0[0]/beta0_)*(outer_product(dTe,dTe));
  
  mat W0_;
  bool ret =inv(iW0_, W0_);
  assert(ret);
  mat tmp = W0_+transpose(W0_)*0.5;
  A_T = wishrnd(tmp, nu0_);

  mat tmp2 ;
  ret = inv(beta0_*A_T, tmp2);
  assert(ret);
  mu_T = mvnrndex(mu0_, tmp2, ac.D);
  if (ac.debug)
    cout<<"Sampling from T: A_T" <<A_T<<" mu_V: "<<mu_T<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
   
}

/**
 * for BPTF: sample hyperprior and noise level at the end of each round
 */
template<typename graph_type, typename vertex_data, typename edge_data>
void last_iter_bptf(double res){
    if (ps.iiter == ac.bptf_burn_in){
      printf("Finished burn-in period. starting to aggregate samples\n");
    timer t;
    t.start();
    if (ps.iiter > ac.bptf_delay_alpha)
    	sample_alpha(res);
    sample_U<graph_type>();
    sample_V<graph_type>();
    if (ps.tensor) 
      sample_T();
    ps.counter[BPTF_SAMPLE_STEP] += t.current_time();
    if (ac.datafile == "kddcup" || ac.datafile == "kddcup2")
	export_kdd_format<graph_type, vertex_data, edge_data>(*ps.g<graph_type>(TEST), TEST, false);
    }
}



#endif //__BPTF_H
