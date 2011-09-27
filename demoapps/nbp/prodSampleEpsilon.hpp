/* Copyright (c) 2003 Alexander Ihler
 * Original code from: http://www.ics.uci.edu/~ihler/code/index.html
 *
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
 ***********************************************************************
 ** multi-tree approximate sampling MEX code
 **
 **
 ***********************************************************************/
//
// Written by Alex Ihler and Mike Mandel
// Copyright (C) 2003 Alexander Ihler
// Converted to C++ by Danny Bickson, CMU, 2010

#ifndef PROD_SAMPL_EPS
#define PROD_SAMPL_EPS


#include "fakemex.h"
#include <itpp/itstat.h>
#include <itpp/itbase.h>
#include <itpp/base/sort.h>
#include "kde.h"

#include "cpp/BallTreeDensity.h"
#include "prob.hpp"

class prodSampleEpsilon{

public:


  // a little addressing formula: 
  //   to access a^th dimension of density pair (b,c)'s constant
#define SIGVALSMAX(a,b,c) (SigValsMax + a+Ndim*b+Ndim*Ndens*c)
#define SIGVALSMIN(a,b,c) (SigValsMin + a+Ndim*b+Ndim*Ndens*c)
  double *SigValsMax, *SigValsMin;

  //BallTreeDensity *trees;    // structure of all trees
  std::vector<BallTreeDensity> trees;
  BallTree::index *ind;      // indices of this level of the trees

  double *C,*sC,*M;

  double *randunif1, *randunif2, *randnorm;  // required random numbers
  double *samples;
  BallTree::index* indices;    // return data

  double maxErr;                 // epsilon tolerance (%) of algorithm
  double total, soFar, soFarMin; // partition f'n and accumulation

  unsigned int Ndim,Ndens;   // useful constants
  unsigned long Nsamp;
  bool bwUniform ;

  prodSampleEpsilon(){
    SigValsMin = SigValsMax = 0;
    ind = 0; 
    C = sC = M = 0;
    randunif2 = randunif1 = randnorm = 0;
    samples = 0; indices = 0;
    maxErr = 0; total = 0; soFarMin =0; soFar = 0;
    Ndim = 0; Ndens = 0; Nsamp = 0; bwUniform = true;
  }

  ~prodSampleEpsilon(){

    mxFree(C); mxFree(sC); mxFree(M); mxFree(SigValsMin); mxFree(SigValsMax);
  }

#ifdef MEX
  //////////////////////////////////////////////////////////////////////
  // MEX WRAPPER
  //////////////////////////////////////////////////////////////////////
  void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
  {
    mxArray *rNorm, *rUnif1, *rUnif2, *rsize;
    unsigned int i,j;
  
    /*********************************************************************
     ** Verify arguments and initialize variables
     *********************************************************************/

    if (nrhs != 3)
      mexErrMsgTxt("Takes 3 input arguments");
    if (nlhs >  2)
      mexErrMsgTxt("Outputs 2 results");

    Ndens = mxGetN(prhs[0]);                               // get # of densities
    //  trees = new BallTreeDensity[Ndens];
    trees = (BallTreeDensity*) mxMalloc(Ndens*sizeof(BallTreeDensity));
    bwUniform = true;
    bool allGaussians = true;
    for (i=0;i<Ndens;i++) {                               // load densities
      trees[i] = BallTreeDensity( mxGetCell(prhs[0],i) );  
      if (trees[i].getType() != BallTreeDensity::Gaussian) allGaussians = false;
      bwUniform = bwUniform && trees[i].bwUniform();
    }
    if (!allGaussians)
      mexErrMsgTxt("Sorry -- only Gaussian kernels supported");

    Ndim  = trees[0].Ndim();                      // more accessible dimension variable
    Nsamp = (unsigned long) mxGetScalar(prhs[1]); // # of requested samples
    maxErr= 2*mxGetScalar(prhs[2]);               // epsilon (we always use 2*epsilon)

    // Obtain enough random numbers for the sampling algorithm
    //
    rsize = mxCreateDoubleMatrix(1,2,mxREAL);
    double* rsizeP= mxGetPr(rsize); rsizeP[0] = 1; rsizeP[1] = Nsamp+1;
    rUnif1 = mxCreateDoubleMatrix(1,Nsamp+1,mxREAL);
    mexCallMATLAB(1, &rNorm, 1, &rsize, "rand");   randunif1 = mxGetPr(rNorm);
    randunif1[Nsamp] = 100;
    mexCallMATLAB(1, &rUnif1, 1, &rNorm, "sort");  randunif1 = mxGetPr(rUnif1);
    mxDestroyArray(rNorm);
    rsizeP[0] = Ndens; rsizeP[1] = Nsamp;
    mexCallMATLAB(1, &rUnif2, 1, &rsize, "rand");  randunif2 = mxGetPr(rUnif2);
    rsizeP[0] = Ndim; rsizeP[1] = Nsamp;
    mexCallMATLAB(1, &rNorm, 1, &rsize, "randn");  randnorm  = mxGetPr(rNorm);

    plhs[0] = mxCreateDoubleMatrix(Ndim,Nsamp,mxREAL);
    samples = (double*) mxGetData(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(Ndens,Nsamp,mxUINT32_CLASS,mxREAL);
    indices = (BallTree::index*) mxGetData(plhs[1]);

    SigValsMax = (double*) mxMalloc(Ndim*Ndens*Ndens*sizeof(double));  // precalc'd constants
    SigValsMin = (double*) mxMalloc(Ndim*Ndens*Ndens*sizeof(double));  // precalc'd constants
    C       = (double*) mxMalloc(Ndim*sizeof(double));
    sC      = (double*) mxMalloc(Ndim*sizeof(double));
    M       = (double*) mxMalloc(Ndim*sizeof(double));
  
    total =    -1; soFar = soFarMin = 0;   multiEval();  // calculate total weight
    total = soFar; soFar = soFarMin = 0;   multiEval();  //   then sample

    //  delete[] trees;
    mxFree(trees);

    mxFree(C); mxFree(sC); mxFree(M); mxFree(SigValsMin); mxFree(SigValsMax);

    mxDestroyArray(rUnif1); mxDestroyArray(rUnif2); 
    mxDestroyArray(rNorm); mxDestroyArray(rsize);
  }
#else
  //////////////////////////////////////////////////////////////////////
  // MEX WRAPPER
  //////////////////////////////////////////////////////////////////////
  kde prodSampleEpsilonRun(unsigned int _Ndens, //number of densities to product
                           unsigned int _Nsamp,  //number of samples
                           double _maxErr,  //epsilon
                           std::vector<kde>& kdes)
  {
    unsigned int i;//,j;
  
 

    /*********************************************************************
     ** Verify arguments and initialize variables
     *********************************************************************/

    bool debug = false;

    Ndens = _Ndens;
    Nsamp = _Nsamp;
    maxErr = _maxErr;
    assert(Ndens >= 1);
    assert(Nsamp>= 1);
    // get # of densities
    //trees = new BallTreeDensity[Ndens];
    //trees = (BallTreeDensity*) mxMalloc(Ndens*sizeof(BallTreeDensity));
    bwUniform = true;
    bool allGaussians = true;
    for (i=0;i<Ndens;i++) {                               // load densities
      //trees[i] = BallTreeDensity( mxGetCell(prhs[0],i) ); 
      trees.push_back(BallTreeDensity( kdes[i] ));  
      if (trees[i].getType() != BallTreeDensity::Gaussian) allGaussians = false;
      bwUniform = bwUniform && trees[i].bwUniform();
      //assert(kdes[i].getPoints() < 13);
    }
    if (!allGaussians)
      mexErrMsgTxt("Sorry -- only Gaussian kernels supported");

    Ndim  = trees[0].Ndim();                      // more accessible dimension variable
 
    itpp::vec rUnif1 = zeros(Nsamp+1);
    randv(Nsamp+1, rUnif1);
    rUnif1.set(Nsamp,100);
    itpp::Sort<double> mysort;
    mysort.sort(0,Nsamp,rUnif1);
    randunif1 = vec2vec(&rUnif1); 
    
    itpp::vec rUnif2 = zeros(Nsamp*Nsamp);
    randv(Ndens*Nsamp, rUnif2);
    randunif2 = vec2vec(&rUnif2);
  
    itpp::mat rNorm = zeros(Ndim*Nsamp);
    rNorm = randn1(Ndim, Nsamp);
    randnorm = vec2vec(&rNorm);

    kde out;
    out.centers = itpp::zeros(Ndim, Nsamp); 
    samples = out.centers._data();
    out.indices = itpp::imat(Ndens, Nsamp);
    indices = (BallTree::index*) out.indices._data();

    SigValsMax = (double*) mxMalloc(Ndim*Ndens*Ndens*sizeof(double));  // precalc'd constants
    SigValsMin = (double*) mxMalloc(Ndim*Ndens*Ndens*sizeof(double));  // precalc'd constants
    C       = (double*) mxMalloc(Ndim*sizeof(double));
    sC      = (double*) mxMalloc(Ndim*sizeof(double));
    M       = (double*) mxMalloc(Ndim*sizeof(double));
  
    total =    -1; soFar = soFarMin = 0;   multiEval();  // calculate total weight
    total = soFar; soFar = soFarMin = 0;   multiEval();  //   then sample

    out.ROT();
    out.weights = itpp::ones(1, out.getPoints())/(double)out.getPoints();
    if (debug)
      out.matlab_print();
    //else {printf("."); fflush(NULL);}
    out.indices = out.indices - 1; //c++ count starts from zero
    out.verify();
 
 
    for (size_t i=0; i < trees.size(); i++)
      trees[i].clean();
    trees.clear(); 
 
    return out;
  }

#endif


  double normConstant(void) {
    unsigned int i,j;
    double tmp, normConst;
    //const double pi = 3.141592653589;
  
    normConst = 1;                               // precalculate influence of normalization
    tmp = pow(2*pi,((double)Ndim)/2);
    for (i=0;i<Ndens;i++) {                      // divide by norm fact of each indiv. gauss.
      normConst /= tmp;
      if (bwUniform) for (j=0;j<Ndim;j++) {
          normConst /= sqrt(trees[i].bwMin(0)[j]);
        }
    }
    normConst *= tmp;                            // times norm factor of resulting gaussian
    for (j=0;j<Ndim;j++) {
      tmp = 0;
      if (bwUniform) {
        for (i=0;i<Ndens;i++) tmp += 1/trees[i].bwMin(0)[j];  // compute result bandwidth
        normConst /= sqrt(tmp);                               // and its norm factor
      }
    }
    return normConst;
  }


  //////////////////////////////////////////////////////////////////////////
  // calculate bounds on the min/max distance possible between two ball-trees
  //   return un-exponentiated values
  //
  double minDistProd(const BallTreeDensity& bt1, BallTree::index i,
                     const BallTreeDensity& bt2, BallTree::index j,
                     const double* SigValIJ,const double* SigNIJ) //  precomp'd weighting factors
  {
    double result=0;
    const double *center1, *center2;

    center1 = bt1.center(i); center2 = bt2.center(j);
    for (unsigned int k=0;k<Ndim;k++) {
      double tmp = fabs( center1[k] - center2[k] );
      tmp-= bt1.range(i)[k] + bt2.range(j)[k];
      if (tmp < 0) tmp = 0;
      result -= (tmp*tmp) * SigValIJ[k];
      if (!bwUniform) result += log(SigNIJ[k]);
    }
    result /= 2;
    return result;
  }

  double maxDistProd(const BallTreeDensity& bt1, BallTree::index i,
                     const BallTreeDensity& bt2, BallTree::index j,
                     const double* SigValIJ,const double* SigNIJ) //  precomp'd weighting factors
  {
    double result=0;
    const double *center1, *center2;

    center1 = bt1.center(i); center2 = bt2.center(j);
    for (unsigned int k=0;k<Ndim;k++) {
      double tmp = fabs( center1[k] - center2[k] );
      tmp+= bt1.range(i)[k] + bt2.range(j)[k];
      result -= (tmp*tmp) * SigValIJ[k];
      if (!bwUniform) result += log(SigNIJ[k]);
    }
    result /= 2;
    return result;
  }

  // Compute (1 over) the \Lambda_(i,j) values needed for distance-weight computations
  // 
  void computeSigVals(void) {
    unsigned int i,j,k;
    assert(Ndim > 0);
    double *SigNormMin = (double*) mxMalloc(Ndim*sizeof(double));
    double *SigNormMax = (double*) mxMalloc(Ndim*sizeof(double));
    for (i=0;i<Ndim;i++) {
      SigNormMin[i] = SigNormMax[i] = 0;
      for (j=0;j<Ndens;j++) SigNormMin[i]+=1/trees[j].bwMin(ind[j])[i]; // compute \Lambda_L 
      for (j=0;j<Ndens;j++) SigNormMax[i]+=1/trees[j].bwMax(ind[j])[i]; //
      SigNormMax[i] = 1/SigNormMax[i]; SigNormMin[i] = 1/SigNormMin[i];
    }
    for (i=0;i<Ndim;i++) {
      for (j=0;j<Ndens;j++)                                    //  then compute pairwise leave-
        for (k=j;k<Ndens;k++) {                                //  two-out normalized values
          *SIGVALSMIN(i,k,j) = SigNormMax[i] / (trees[j].bwMin(ind[j])[i]*trees[k].bwMin(ind[k])[i]);
          *SIGVALSMAX(i,k,j) = SigNormMin[i] / (trees[j].bwMax(ind[j])[i]*trees[k].bwMax(ind[k])[i]);
          *SIGVALSMIN(i,j,k) = *SIGVALSMIN(i,k,j);             //  make symmetric
          *SIGVALSMAX(i,j,k) = *SIGVALSMAX(i,k,j);
        }
    }
    //  delete[] SigNorm;  //(don't need this anymore)
    mxFree(SigNormMin);
    mxFree(SigNormMax);

  }

  void multiEvalRecursive(void) {
    unsigned int i,j;
    double minVal=0, maxVal=0;                    // for computing bounds and 
    unsigned int maxInd0, maxInd1;  //  determining which tree to split

    //
    // find min/max values of product
    //
    if (!bwUniform) computeSigVals();

    double maxDiscrep = -1;
    bool allLeaves = true;
    for (i=0; i<Ndens; i++) {                       // For each pair of densities, bound
      for (j=i+1;j<Ndens;j++) {                     //   the total weight of their product:
        double maxValT = minDistProd(trees[i],ind[i],trees[j],ind[j],SIGVALSMAX(0,i,j),SIGVALSMIN(0,i,j));  // compute min & max
        double minValT = maxDistProd(trees[i],ind[i],trees[j],ind[j],SIGVALSMIN(0,i,j),SIGVALSMAX(0,i,j));  // dist = max/min values
        maxVal += maxValT; minVal += minValT;

        if ((maxValT - minValT) > maxDiscrep) {           // also find which pair
          maxDiscrep = maxValT - minValT;                 //   has the largest
          maxInd0=i; maxInd1=j;                           //   discrepancy (A/B)
        }
      }
      allLeaves = allLeaves && trees[i].isLeaf(ind[i]);
    }
    maxVal = exp(maxVal); minVal = exp(minVal);

    // If the approximation is good enough,
    if (allLeaves || fabs(maxVal - minVal) <= maxErr * (soFarMin+minVal) ) {  // APPROXIMATE
      double add = (maxVal + minVal)/2;                   // compute contribution
      for (i=0;i<Ndens;i++) add *= trees[i].weight(ind[i]);
      soFar += add;
      add = minVal; for (i=0;i<Ndens;i++) add *= trees[i].weight(ind[i]);
      soFarMin += add;

      while (*randunif1 <= soFar/total) {                 // for all the samples coming from this block
        randunif1++;
        for (j=0;j<Ndim;j++) M[j] = 0;                    // clear out M
        if (!bwUniform) for (j=0;j<Ndim;j++) C[j] = 0;    // clear out C if necc.

        for (i=0;i<Ndens;i++) {                           // find an index within this block
          double SumTmp = 0;
          BallTree::index index = trees[i].leafFirst(ind[i]);  // start with 1st leaf and
          for (;index <= trees[i].leafLast(ind[i]);index++) {
            SumTmp += trees[i].weight(index) / trees[i].weight(ind[i]);
            if (SumTmp > *randunif2) break;
          }
          randunif2++;
          for (j=0;j<Ndim;j++)                                 // compute product mean:
            M[j] += trees[i].center(index)[j] / trees[i].bw(index)[j];
          *(indices++) = trees[i].getIndexOf(index)+1;         // and save selected indices
          //assert(trees[i].getIndexOf(index)+1 <= Ndens*Nsamp);
          if (!bwUniform) for (j=0;j<Ndim;j++)                 // compute covariance
                            C[j] += 1/trees[i].bw(index)[j];                 //  contribution of each dens.
        }
        if (!bwUniform) for (j=0;j<Ndim;j++) {                 // finish computing covar and
            C[j] = 1/C[j];                                     //  std dev. of product kernel
            sC[j] = sqrt(C[j]);
          }

        for (j=0;j<Ndim;j++) M[j] *= C[j];
        for (j=0;j<Ndim;j++)                              // sample from the product dist.
          *(samples++) = M[j] + sC[j] * (*(randnorm++));
      }

      // Otherwise, we need to subdivide at least one tree:
    } else {                                              // RECURSION  
      unsigned int split;
      double size0 = trees[maxInd0].range(ind[maxInd0])[0];  // from the pair with the largest
      double size1 = trees[maxInd1].range(ind[maxInd1])[0];  // pairwise max-min discrepancy term,

      for(BallTree::index k=0; k<trees[maxInd0].Ndim(); k++)
        if(trees[maxInd0].range(ind[maxInd0])[k] > size0)
          size0 = trees[maxInd0].range(ind[maxInd0])[k];
      for(BallTree::index k=0; k<trees[maxInd1].Ndim(); k++)
        if(trees[maxInd1].range(ind[maxInd1])[k] > size1)
          size1 = trees[maxInd1].range(ind[maxInd1])[k];    

      split = (size0 > size1) ? maxInd0 : maxInd1;        // take the largest.
    
      BallTree::index current = ind[split];
      if (!trees[split].isLeaf(current)) {
        ind[split] = trees[split].left(current);  
        multiEvalRecursive();                             // recurse left 
        ind[split] = trees[split].right(current);         //   and right tree
        multiEvalRecursive();                             // restore indices 
        ind[split] = current;                             //   for calling f'n
      }                                                   
    }
  }


  void multiEval(void) {
    unsigned int i,j;//,k;
    //  ind = new BallTree::index[Ndens];               // construct index array  

    assert(Ndens>0);
    ind = (BallTree::index*) mxMalloc(Ndens*sizeof(BallTree::index));    // construct index array  
    memset(ind, 0, Ndens * sizeof(BallTree::index));
    for (i=0;i<Ndens;i++) ind[i] = trees[i].root(); //  & init to root node

    if (bwUniform) {                                     // if all one kernel size, do this in
      computeSigVals();                                  //   one operation.
      for (i=0;i<Ndim;i++) {                             // compute covariance and
        double tmp = 0;                                  //   std. deviation of a
        for (j=0;j<Ndens;j++)                            // resulting product kernel 
          tmp += 1/trees[j].bw(trees[j].leafFirst(trees[j].root()))[i]; 
        C[i] = 1/tmp;
        sC[i] = sqrt(C[i]);
      }
    }

    multiEvalRecursive();

    //  delete[] ind;
    mxFree(ind);
  }
}; //class
inline void test_product(){
  printf("testing product..\n");
  kde k = kde("3 1", "1 1", "1 2");
  kde j = kde("2", ".5", "1");
  std::vector<kde> vecs;
  vecs.push_back(k);
  vecs.push_back(j);
  prodSampleEpsilon prod;
  kde out = prod.prodSampleEpsilonRun(2,48,1e-5,vecs);
  out.matlab_print();
  out.verify();
}




#endif
