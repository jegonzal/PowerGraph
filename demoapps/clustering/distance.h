#ifndef _DISTANCE_H
#define _DISTANCE_H

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "graphlab.hpp"


using namespace itpp;



enum distance_measure{
   EUCLIDEAN = 0,
   CHEBYCHEV = 1,
   MANAHOLIS = 2,
   MANHATTAN = 3,
   MINKOWSKI = 4,
   TONIMOTO = 5,
   WEIGHTED = 6,
   WEIGHTED_MANAHOLIS = 7,
   COSINE = 8
};


double calc_distance(itpp::sparse_vec & datapoint, itpp::vec &cluster, double sqr_sum = 0);


#endif //_DISTANCE_H
