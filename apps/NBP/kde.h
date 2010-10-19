#ifndef __KDE_H
#define __KDE_H


#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include "assert.h"
#include <vector>

typedef itpp::Mat<unsigned int> uimat;

inline double square(const double v) {
  return v * v;
}

inline double log_likelihood(double u, double sigma, double x){
  return log(1.0/sigma) - square(x - u) / (2 * sigma * sigma);
}

class kde{
public:
	itpp::mat centers;
	itpp::mat bw;
        itpp::vec weights;
        itpp::mat indices;

        kde(){};

	kde(itpp::mat &_centers, itpp::mat &_bw, itpp::vec &_weights){
	    centers = _centers;
            bw = _bw;
            weights = _weights;
            normalize_weights();
        }
	kde(const char * _centers, const char * _bw, const char * _weights){
	    centers = itpp::mat(_centers);
            bw = itpp::mat(_bw);
            weights = itpp::vec(_weights);
            normalize_weights();
        }
		kde(itpp::mat &_centers, itpp::mat &_bw, itpp::mat &_weights){
	    centers = _centers;
            bw = _bw;
            assert(_weights.rows() == 1);
            weights = _weights.get_row(0);
            normalize_weights();
        }

	kde(itpp::mat &_centers, itpp::mat &_bw){
	    centers = _centers;
            bw = _bw;
        }

        kde marginal(int dim){
	    assert(dim < centers.rows());
            assert(dim >= 0);
            itpp::mat slice = centers(dim, dim, 0, centers.cols()-1);
            itpp::mat bwslice = bw;
            if (bw.rows() > 1)
                 bwslice = bw(dim, dim, 0, bw.cols() - 1); 
            return kde(slice, bwslice, weights); 
        }

       void normalize_weights(){
           double sum = itpp::sum(weights);
           assert(sum > 0);
           weights = weights/ sum;
       }


       //L = evaluate(kde,getPoints(kde));
       //%[mx,mxind] = max(L);
       //X = getPoints(kde,mxind(1));
        double max() const{
           assert(getDim() == 1);
           double max = -1e100;
           int pos = -1;
           for (int i=0; i<getPoints(); i++){
                double like = likelihood(centers(i));
                if (max < like){
		   max = like;
                   pos = i;
                }
                  
           }
           assert(pos >= 0);
           return centers(pos);
        }

       

        double likelihood(const double &d) const{
           double ret = 0;
           for (size_t j = 0;j < (size_t)getPoints(); ++j) {
             double curll = log_likelihood(centers(j), bw(j), d);
             ret += weights(j) * exp(curll);
           }
           return ret;
       }

 
        void verify() const{
	    assert(sum(weights) > 0);
            assert(sumsum(bw) > 0);
            assert(itpp::min(itpp::min(bw))>0);
            assert(centers.rows() > 0);
            assert(centers.cols() > 0);
            assert(itpp::max(itpp::max(centers)) < 1e10);
            assert(itpp::max(itpp::max(bw)) < 1e10);
            assert(weights.size() == getPoints());
            assert(bw.cols() == centers.cols()); 
       }

        int getDim() const{
            return centers.rows();
        }

        int getPoints() const{
            return centers.cols();
        }

        // points = pts(:,ind) + getBW(npd,ind).*randKernel(getDim(npd),length(ind),getType(npd));
         kde sample(itpp::mat & ind,itpp::vec & weights){
             assert(sum(weights)>0);
             assert(itpp::max(itpp::max(ind)) < centers.cols());
             assert(min(min(ind)) >= 0);
             itpp::mat randN; 
             itpp::randn(getDim(), ind.size(), randN);
             itpp::mat pts = itpp::zeros(centers.rows(), ind.size());
             itpp::mat pbw = itpp::zeros(centers.rows(), ind.size());
             for (int i=0; i< centers.rows(); i++){
                for (int j=0; j< ind.size(); j++){
                   pts.set(i,j,centers(i,ind(j)));
                   pbw.set(i,j,bw(ind(j)));
                }
             }
             itpp::mat points = pts + elem_mult(pbw, randN);
             return kde(points, pbw, weights);
                    
         }

      static void matlab_print(const itpp::vec & data){
          for (int i=0; i< data.size(); i++)
           	std::cout<<" "<<data(i);
      }

      static void matlab_print(const itpp::mat & data){
          std::cout<<"[";
          for (int i=0; i< data.rows(); i++){
             matlab_print(data.get_row(i));
             if (i < data.rows() -1 )
                 std::cout<<";";
          }
          std::cout<<"]";
      }



      void matlab_print(){
          std::cout<<"kde("; 
          matlab_print(centers);
          std::cout<< ",";
          matlab_print(bw); 
          std::cout<< ",[";
          matlab_print(weights);
          std::cout << "]);" << std::endl;
      }
      
/*function h = ksizeROT(npd,noIQR)
% "Rule of Thumb" estimate (Silverman)
%    Estimate is based on assumptions of Gaussian data and kernel
%    Actually the multivariate version in Scott ('92) 
%  Use ksizeROT(X,1) to force use of stddev. instead of min(std,C*iqr)
%       (iqr = interquartile range, C*iqr = robust stddev estimate)
%

% Copyright (C) 2003 Alexander Ihler; distributable under GPL -- see README.txt

  X = getPoints(npd);
  N = size(X,2);  dim = size(X,1);
  if (nargin<2) noIQR=0; end;

  Rg = .282095; Mg=1;                     % See ksizeCalcUseful for derivation
  Re = .6;      Me = .199994;             %   this is the canonical kernel adjustment
  Rl = .25;     Ml = 1.994473;            %   for product kernels of these types
  switch(npd.type),
      case 0, prop = 1.0;                 % Approximate; 1D prop = 1.059224; % Gaussian
      case 1, prop = ((Re/Rg)^dim / (Me/Mg)^2 )^(1/(dim+4)); % 1D prop = 2.344944; % Epanetchnikov
      case 2, prop = ((Rl/Rg)^dim / (Ml/Mg)^2 )^(1/(dim+4)); % 1D prop = 0.784452; % Laplacian
  end;
  
  sig = std(X,0,2);            % estimate sigma (standard)
  if (noIQR)
    h = prop*sig*N^(-1/(4+dim));
  else  
    iqrSig = .7413*iqr(X')';     % find interquartile range sigma est.
    if (max(iqrSig)==0) iqrSig=sig; end;
    h = prop * min(sig,iqrSig) * N^(-1/(4+dim));
  end;

%
*/

        void ROT(){
            assert(getDim() == 1);
            assert(getPoints() > 1); //no meaning to compute variance over one point
            double prop = 1.0;
            double sig = sqrt(itpp::variance(centers.get_row(0)));
            assert(!isnan(sig));
            assert(sig > 0);
            double h = prop*sig*powf(getPoints(),(-1.0/(4.0+getDim())));
            bw = itpp::ones(1,getPoints()) * h; 
        }

  };

kde prodSampleEpsilon(unsigned int Ndens, //number of densities to product
		       unsigned int Nsamp,  //number of samples
                       double maxErr,  //epsilon
                       std::vector<kde>& kdes);//kdes to product


      inline void test_marginal(){ 
          printf("testing marginal..\n");
          itpp::mat mcenters = "1 2 3; 3 2 1";
          itpp::mat mbw = "0.5 0.5 0.2; 0.5 0.5 0.2";
          itpp::vec weights = "0.2 0.3 0.4";
          kde k = kde(mcenters, mbw, weights);
          assert(k.getDim() == 2);
          assert(k.getPoints() == 3); 
          k.verify();
          k.matlab_print();
          kde k1 = k.marginal(0);
          k1.matlab_print();
          assert(k1.centers.get_row(0) == itpp::vec(" 1 2 3"));
          assert(k1.bw.get_row(0) == itpp::vec(".5 .5 .2"));
          assert(square(k1.weights(0) - 0.22222) < 1e-8);
 
          kde k2 = k.marginal(1);
          assert(k2.centers.get_row(0) == itpp::vec(" 3 2 1"));
          assert(k2.bw.get_row(0) == itpp::vec(".5 .5 .2"));
          assert(square(k2.weights(0) - 0.22222) < 1e-8);
          k2.matlab_print();
      }

       inline void test_max(){ 
          printf("testing max..\n");
          itpp::mat mcenters = " 1  2   3    -1  -2  3    2   1";
          itpp::mat mbw =      "0.5 0.5 0.2  0.5 0.5 0.2  3   2";
          itpp::vec weights =  "0.2 0.3 0.4  0.1 0.05 0.05 0.05 0.05";
          kde k = kde(mcenters, mbw, weights);
          k.matlab_print();
          std::cout<<k.max()<<std::endl;
          assert(k.max() == 3);
	  for (int i=0; i< mcenters.cols(); i++)
          	std::cout<<"i:"<<i<<" "<< k.likelihood(mcenters(i))<<std::endl;
          kde k1 = kde("3 3 2", "3 3 2", "1 1 1");
          assert(k1.max() == 2);
          kde k2 = kde("3 3 2", " 1 1 2", "1 1 1");
          assert(k2.max() == 3);
        }

        inline void test_sample(){
          printf("testing sample..\n");
           itpp::mat cent = "0"; itpp::mat bw= "1"; itpp::mat weight = "1";
           kde k(cent, bw, weight);
           k.matlab_print();
           double sum = 0;
           itpp::mat ind = "0";
           itpp::vec vweight = "1";
           for (int i=0; i< 10000; i++){
              kde out = k.sample(ind, vweight);
              sum += out.centers(0);
           }
           std::cout<<" mean is: " << sum/10000 << " should be: 0 "<< std::endl;

        }
       inline void test_sample2(){
          printf("testing sample2..\n");
           itpp::mat cent = "1 2 3 1 -1 2"; itpp::mat bw= "1 0.5 0.1 0.01 3 2"; itpp::mat weight = "0.5 0.2 0.1 0.1 0.1 0.1";
           kde k(cent, bw, weight);
           k.matlab_print();
           double sum = 0;
           itpp::mat ind =     "0 1 3 2 3 2 3 3 2 1 4 5";
           itpp::vec vweight = ".5 .5 .2 .1 .2 .1 .2 .1 .1 .05 .05";
           for (int i=0; i< 10000; i++){
              kde out = k.sample(ind, vweight);
              sum += itpp::sum(itpp::sum(out.centers));
           }
           std::cout<<" mean is: " << sum/(11*01000) << std::endl;

        }

        inline void test_ROT(){
          printf("testing ROT..\n");
          kde k1 = kde("3 3 2", "3 3 2", "1 1 1");
          k1.matlab_print();
          k1.ROT();
          k1.matlab_print();
          assert(square(k1.bw(0) - 0.4635)<1e-8);
        }

 
        inline void test_product(){
          printf("testing product..\n");
	   kde k = kde("3 1", "1 1", "1 2");
	   kde j = kde("2", ".5", "1");
           std::vector<kde> vecs;
           vecs.push_back(k);
           vecs.push_back(j);
           kde out = prodSampleEpsilon(2,48,1e-5,vecs);
           out.matlab_print();
           out.verify();
         }



#endif
