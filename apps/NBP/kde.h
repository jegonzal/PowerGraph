#ifndef __KDE_H
#define __KDE_H


#include <itpp/itbase.h>
#include "assert.h"

typedef itpp::Mat<unsigned int> uimat;

class kde{
public:
	itpp::mat centers;
	itpp::vec bw;
        itpp::vec weights;
        uimat indices;

        kde(){};

	kde(itpp::mat &_centers, itpp::vec &_bw, itpp::vec &_weights){
	    centers = _centers;
            bw = _bw;
            weights = _weights;
        }

        kde marginal(int dim){
	    assert(dim < centers.rows());
            assert(dim >= 0);
            itpp::mat slice = centers(dim, dim, 0, centers.cols());
            return kde(slice, bw, weights); 
        }

        void verify(){
	    assert(sum(weights) > 0);
            assert(sum(bw) > 0);
        }

        int getDim(){
            return centers.rows();
        }

        // points = pts(:,ind) + getBW(npd,ind).*randKernel(getDim(npd),length(ind),getType(npd));
         kde sample(itpp:vec ind){
             assert(max(ind) < centers.cols());
             assert(min(ind) >= 0);
             mat randN; 
             itpp::randn(getDim(), ind.size(), randN);
             mat pts = zeros(centers.rows(), ind.size());
             mat pbw = zeros(centers.rows(), ind.size());
             for (int i=0; i< centers.rows(); i++){
                for (int j=0; j< ind.size(); j++){
                   pts.set(i,j,centers(i,ind(j)));
                   pbw.set(i,j,bw(centers.bw(ind(j))));
                }
             }
             for (int j=0; j< ind.size(); j++){
             }
             mat points = pts + pbw .* randN;
             return kde(points, bpw);
                    
         }
       
};


#endif
