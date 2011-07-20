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

//external global variables defined at pmf.cpp
extern bool debug,tensor;
extern bool BPTF;
extern int M,N,K,L,BURN_IN;
extern vertex_data * times;
extern double counter[20];
extern graph_type * g;
extern graph_type test_graph;
extern int iiter;
extern int delayalpha;
extern string infile;

/* variables for BPTF */
double nuAlpha = 1;
double Walpha = 1;
double mu0 = 0;
double mu0T = 1;
double nu0 = D;
double alpha = 0;
double beta = 1;
double beta0 = 1; //TODO
mat W0;
mat W0T;
double iWalpha;
mat iW0;
mat iW0T;
mat A_U, A_V, A_T;
vec mu_U, mu_V, mu_T;

using namespace itpp;
using namespace graphlab;

void export_kdd_format(graph_type * _g, testtype type, bool dosave);
/**
* calc the A'A part of the least squares solution inv(A'A)*A'
* for time nodes
*/
mat calc_MMT(int start_pos, int end_pos, vec &Umean){

  int batchSize = 1000;
  mat U = zeros(batchSize,D);
  mat MMT = zeros(D,D);
  int cnt = 0;
  timer t;

  for (int i=start_pos; i< end_pos; i++){
    if ((i-start_pos) % batchSize == 0){
      U=zeros(batchSize, D);
      cnt = 1;
    }

    const vertex_data * data= &g->vertex_data(i);
     
    vec mean = data->pvec;
    Umean += mean;
    //Q.set_col(k, ret);
    t.start(); 
    for (int s=0; s<D; s++)
      U(i%batchSize,s)=mean(s);
    if (debug && (i==start_pos || i == end_pos-1))
      std::cout<<" clmn "<<i<< " vec: " << mean <<std::endl;

    if ((cnt  == batchSize) || (cnt < batchSize && i == end_pos-1)){
      MMT = MMT+transpose(U)*U;
    }
    counter[8] += t.current_time();
    cnt++;
  }
  Umean /= (end_pos-start_pos);
  if (debug)
    cout<<"mean: "<<Umean<<endl;

  assert(MMT.rows() == D && MMT.cols() == D);
  assert(Umean.size() == D);
  return MMT;
}


void init_self_pot(){
  //assert(BPTF);

  W0 = eye(D);
  W0T = eye(D);
  iWalpha = 1.0/Walpha;
  iW0 = inv(W0);
  iW0T = inv(W0T);
  nu0 = D;

  A_U = eye(D); //cov prior for users
  A_V = eye(D); //cov prior for movies
  A_T = eye(D); //cov prior for time nodes

  mu_U = zeros(D); mu_V = zeros(D); mu_T = zeros(D);
  printf("nuAlpha=%g, Walpha=%g, mu=%g, muT=%g, nu=%g, "
         "beta=%g, W=%g, WT=%g BURN_IN=%d\n", nuAlpha, Walpha, mu0, 
         mu0T, nu0, beta0, W0(1,1), W0T(1,1), BURN_IN);


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
  
  if (debug)
  printf("res is %g\n", res2); 
  
  double res = res2;
  assert(BPTF);
  if (nuAlpha > 0){
    double nuAlpha_ =nuAlpha+ L;
    mat iWalpha_(1,1);
    iWalpha_.set(0,0, iWalpha + res);
    mat iiWalpha_ = zeros(1,1);
    iiWalpha_ = inv(iWalpha_);
    alpha = wishrnd(iiWalpha_, nuAlpha_).get(0,0);
    assert(alpha != 0);

    if (debug)
      cout<<"Sampling from alpha" <<nuAlpha_<<" "<<iWalpha<<" "<< iiWalpha_<<" "<<alpha<<endl;
    printf("sampled alpha is %g\n", alpha); 
  }
}



// sample movie nodes hyperprior
// according to equation A.3 in Xiong paper.
void sample_U(){
  assert(BPTF);

  vec Umean;
  mat UUT = calc_MMT(0,M,Umean);
  
  double beta0_ = beta0 + M;
  vec mu0_ = (beta0*mu0 + M*Umean)/beta0_;
  double nu0_ = nu0 +M;
  vec dMu = mu0 - Umean;
  if (debug)
    cout<<"dMu:"<<dMu<<"beta0: "<<beta0<<" beta0_ "<<beta0_<<" nu0_ " <<nu0_<<" mu0_ " << mu0_<<endl;
  mat UmeanT = M*(itpp::outer_product(Umean, Umean));
  assert(UmeanT.rows() == D && UmeanT.cols() == D);
  mat dMuT = (beta0*M/beta0_)*(itpp::outer_product(dMu, dMu));
  mat iW0_ = iW0 + UUT - UmeanT + dMuT;
  mat W0_; 
  bool ret =inv(iW0_, W0_);
  assert(ret);
  mat tmp = (W0_+transpose(W0_))*0.5;
  if (debug)
    cout<<iW0<<UUT<<UmeanT<<dMuT<<W0_<<tmp<<nu0_<<endl;
  A_U = wishrnd(tmp, nu0_);
  mat tmp2;  
  ret =  inv(beta0_ * A_U, tmp2);
  assert(ret);
  mu_U = mvnrndex(mu0_, tmp2, D);
  if (debug)
    cout<<"Sampling from U" <<A_U<<" "<<mu_U<<" "<<Umean<<" "<<W0_<<tmp<<endl;
}

// sample user nodes hyperprior
// according to equation A.4 in Xiong paper
void sample_V(){

  assert(BPTF);
  vec Vmean;
  mat VVT = calc_MMT(M, M+N, Vmean);   

  double beta0_ = beta0 + N;
  vec mu0_ = (beta0*mu0 + N*Vmean)/beta0_;
  double nu0_ = nu0 +N;
  vec dMu = mu0 - Vmean;
  if (debug)
    cout<<"dMu:"<<dMu<<"beta0: "<<beta0<<" beta0_ "<<beta0_<<" nu0_ " <<nu0_<<endl;
  mat VmeanT = N*(itpp::outer_product(Vmean, Vmean));
  assert(VmeanT.rows() == D && VmeanT.cols() == D);
  mat dMuT =  (beta0*N/beta0_)*itpp::outer_product(dMu, dMu);
  mat iW0_ = iW0 + VVT - VmeanT + dMuT;
  mat W0_;
  bool ret = inv(iW0_, W0_);
  assert(ret);
  mat tmp = (W0_+transpose(W0_))*0.5;
  if (debug)
    cout<<"iW0: "<<iW0<<" VVT: "<<VVT<<" VmeanT: "<<VmeanT<<" dMuT: " <<dMuT<<"W0_"<< W0_<<" tmp: " << tmp<<" nu0_: "<<nu0_<<endl;
  A_V = wishrnd(tmp, nu0_);
  mat tmp2; 
  ret = inv(beta0_*A_V, tmp2);
  assert(ret);
  mu_V = mvnrndex(mu0_, tmp2, D);
  if (debug)
    cout<<"Sampling from V: A_V" <<A_V<<" mu_V: "<<mu_V<<" Vmean: "<<Vmean<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
}


//calc laplacian matrix for a chain connection between time nodes (for tensor)
//according to equation A.5 in Xiong paper
mat calc_DT(){

  assert(tensor);

  mat T = zeros(D, K);
  for (int i=0; i<K; i++){
    T.set_col(i,times[i].pvec);
  }
  
  mat diff = zeros(D,K-1);
  for (int i=0; i<K-1; i++){
    diff.set_col(i , T.get_col(i) - T.get_col(i+1));
  }
  if (debug)
    cout<<"T:"<<T<<" diff: " << diff<<endl;
  
  return diff;

}

// sample from time nodes
void sample_T(){
  assert(BPTF);
  assert(tensor);

  double beta0_ = beta0 + 1;
  vec pvec = times[0].pvec; 
  vec mu0_ = (pvec + beta0*mu0T)/beta0_;
  double nu0_ = nu0 +K;
  //vec dMu = mu0 - Umean;
  if (debug){
    cout<<"beta0_ " << beta0_ << " beta0: " << beta0 << " nu0_ " << nu0_ << endl;
  } 

  mat dT = calc_DT();
  vec dTe = pvec - mu0T;
  mat iW0_ = iW0T + dT*transpose(dT) + (beta0/beta0_)*(itpp::outer_product(dTe,dTe));
  
  mat W0_;
  bool ret =inv(iW0_, W0_);
  assert(ret);
  mat tmp = W0_+transpose(W0_)*0.5;
  A_T = wishrnd(tmp, nu0_);

  mat tmp2 ;
  ret = inv(beta0_*A_T, tmp2);
  assert(ret);
  mu_T = mvnrndex(mu0_, tmp2, D);
  if (debug)
    cout<<"Sampling from T: A_T" <<A_T<<" mu_V: "<<mu_T<<" W0_: "<<W0_<<" tmp: "<<tmp<<endl;
   
}

/**
 * for BPTF: sample hyperprior and noise level at the end of each round
 */
void last_iter_bptf(double res){
    if (iiter == BURN_IN){
      printf("Finished burn-in period. starting to aggregate samples\n");
    }
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
	export_kdd_format(&test_graph, TEST, false);
}



#endif //__BPTF_H
