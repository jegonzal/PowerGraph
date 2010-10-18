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

};


#endif
