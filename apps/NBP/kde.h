#ifndef __KDE_H
#define __KDE_H


#include <itpp/itbase.h>
#include "assert.h"

typedef itpp::Mat<unsigned int> uimat;

class kde{
public:
	itpp::mat centers;
	itpp::mat bw;
        itpp::vec weights;
        uimat indices;

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
         kde sample(itpp::vec ind){
             assert(max(ind) < centers.cols());
             assert(min(ind) >= 0);
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
             return kde(points, pbw);
                    
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
       
};


#endif
