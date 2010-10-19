#ifndef __KDE_H
#define __KDE_H


#include <itpp/itbase.h>
#include "assert.h"
#include <vector>

typedef itpp::Mat<unsigned int> uimat;

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
        }
	kde(itpp::mat &_centers, itpp::mat &_bw){
	    centers = _centers;
            bw = _bw;
        }

        kde marginal(int dim){
	    assert(dim < centers.rows());
            assert(dim >= 0);
            itpp::mat slice = centers(dim, dim, 0, centers.cols());
            return kde(slice, bw, weights); 
        }

        void verify(){
	    assert(sum(weights) > 0);
            assert(sumsum(bw) > 0);
        }

        int getDim(){
            return centers.rows();
        }

        // points = pts(:,ind) + getBW(npd,ind).*randKernel(getDim(npd),length(ind),getType(npd));
         kde sample(itpp::mat & ind,itpp::vec & weights){
             assert(sum(weights)>0);
             assert(max(max(ind)) < centers.cols());
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
      
      void test_marginal(){ 
          itpp::mat mcenters = "[ 1 2 3; 3 2 1]";
          itpp::mat mbw = "[0.5 0.5 0.2; 0.5 0.5 0.2]";
          itpp::vec weights = "[0.2 0.3 0.4]";
          kde k = kde(centers, mbw, weights);
          kde k1 = k.marginal(1);
          k1.matlab_print();
          kde k2 = k.marginal(2);
          k2.matlab_print();
      }
 
};

kde prodSampleEpsilon(unsigned int Ndens, //number of densities to product
		       unsigned int Nsamp,  //number of samples
                       double maxErr,  //epsilon
                       std::vector<kde>& kdes);//kdes to product
#endif
