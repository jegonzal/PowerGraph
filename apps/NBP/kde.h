#ifndef __KDE_H
#define __KDE_H


#include <itpp/itbase.h>

class kde{
public:
	itpp::mat centers;
	itp::vec bw;
        itpp::vec weights;

	kde(itpp::mat &_centers, itpp::vec &_bw, itpp::vec &_weights){
	    centers = _centers;
            bw = _bw;
            weights = _weights;
        }

};


#endif
