#ifndef _DISTANCE_H
#define _DISTANCE_H

#include "graphlab.hpp"
#include "clustering.h"
#include "../pmf/mathlayer.hpp"
#include "mathlayerf.hpp"

enum distance_measure{
   EUCLIDEAN = 0,
   CHEBYCHEV = 1,
   MANAHOLIS = 2,
   MANHATTAN = 3,
   MINKOWSKI = 4,
   TANIMOTO = 5,
   WEIGHTED = 6,
   WEIGHTED_MANAHOLIS = 7,
   COSINE = 8,
   LOGLIKELIHOOD = 9
};


flt_dbl calc_distance(sparse_flt_dbl_vec & datapoint, flt_dbl_vec &cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0);
flt_dbl calc_distance(sparse_flt_dbl_vec & datapoint, sparse_flt_dbl_vec &cluster, flt_dbl sqr_sum , flt_dbl sqr_sum0);


#endif //_DISTANCE_H
