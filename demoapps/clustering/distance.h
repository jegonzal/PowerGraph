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

static const char * distance_measure_name[] = {"EUCLIDEAN", "CHEBYCHEV", "MANAHOLIS", "MANHATTAN", "MINKOWSKI", "TONIMOTO", "WEIGTED", "WEIGHTED_MANAHOLIS", "COSINE"};

double calc_distance(itpp::sparse_vec & datapoint, itpp::vec &cluster);


#endif //_DISTANCE_H
