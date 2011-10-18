#include "distance.h"
#include "clustering.h"
#include "../gabp/advanced_config.h"

extern advanced_config ac;
const char * distance_measure_name[] = {"EUCLIDEAN", "CHEBYCHEV", "MANAHOLIS", "MANHATTAN", "MINKOWSKI", "TONIMOTO", "WEIGTED", "WEIGHTED_MANAHOLIS", "COSINE"};

sparse_vec minus(sparse_vec & dvec1, sparse_vec & dvec2);
vec minus(sparse_vec & dvec1, vec & dvec2);
double sum_sqr(sparse_vec & dvec1);
sparse_vec fabs(sparse_vec& dvec1);
double sum(sparse_vec & dvec);

double calc_euclidian_distance( sparse_vec & datapoint,  sparse_vec &cluster, double sqr_sum, double sqr_sum0){
  //sparse_vec diff = minus(datapoint , cluster);
  //return sqrt(sum_sqr(diff));
  sparse_vec mult = elem_mult(datapoint, cluster);
  double diff = (sqr_sum + sqr_sum0 - 2*sum(mult));
  return sqrt(fabs(diff)); //because of numerical errors, diff may be negative
}


double calc_euclidian_distance( sparse_vec & datapoint,  vec &cluster, double sqr_sum, double sqr_sum0){
  double dist = sqr_sum + sqr_sum0;
  //for (int i=0; i< datapoint.nnz(); i++){
  FOR_ITERATOR(i, datapoint){
      double val = get_nz_data(datapoint, i);
      int pos = get_nz_index(datapoint, i);
      dist -= 2*val*cluster[pos];
   }
  if (fabs(dist) <1e-10)
     dist = 0;
  return sqrt(dist);
}

double calc_chebychev_distance( sparse_vec & datapoint,  sparse_vec &cluster){
   sparse_vec diff = minus(datapoint , cluster);
   double ret = 0;
   FOR_ITERATOR(i, diff){
      ret = std::max(ret, fabs(get_nz_data(diff, i)));
   }
   return ret;

}
double calc_chebychev_distance( sparse_vec & datapoint,  vec &cluster){
   vec diff = minus(datapoint , cluster);
   double ret = 0;
   for (int i=0; i< diff.size(); i++)
      ret = std::max(ret, fabs(diff[i]));

   return ret;

}

double calc_manhatten_distance( sparse_vec & datapoint,  sparse_vec &cluster){
   sparse_vec diff = minus(datapoint , cluster);
   sparse_vec absvec = fabs(diff);
   double ret = sum(absvec);
   return ret;

}
double calc_manhatten_distance( sparse_vec & datapoint,  vec &cluster){
   vec diff = minus(datapoint , cluster);
   double ret = sum(abs(diff));
   return ret;

}

double calc_cosine_distance( sparse_vec & datapoint,  sparse_vec & cluster, double sum_sqr, double sum_sqr0){
   double dotprod = dot_prod(datapoint,cluster);
   double denominator = sqrt(sum_sqr0)*sqrt(sum_sqr);
   return 1.0 - dotprod / denominator; 
}

double calc_cosine_distance( sparse_vec & datapoint,  vec & cluster, double sum_sqr, double sum_sqr0){
   double dotprod = dot_prod(datapoint,cluster);
   double denominator = sqrt(sum_sqr0)*sqrt(sum_sqr);
   return 1.0 - dotprod / denominator; 
}


double calc_distance(sparse_vec &datapoint,  sparse_vec & cluster, double sqr_sum, double sqr_sum0){
   switch(ac.distance_measure){
      case EUCLIDEAN:          
          return calc_euclidian_distance(datapoint, cluster, sqr_sum, sqr_sum0);
      case CHEBYCHEV:
          return calc_chebychev_distance(datapoint, cluster);
      case COSINE:
	  return calc_cosine_distance(datapoint, cluster, sqr_sum, sqr_sum0);  
      case MANHATTAN:
          return calc_manhatten_distance(datapoint, cluster);
      case MANAHOLIS:
      case WEIGHTED_MANAHOLIS:
      case WEIGHTED:
      default:
          logstream(LOG_ERROR)<< "distance measure " << ac.distance_measure<< "  not implemented yet" << std::endl;
    }
    return -1;

}


double calc_distance(sparse_vec &datapoint,  vec & cluster, double sqr_sum, double sqr_sum0){
   switch(ac.distance_measure){
      case EUCLIDEAN:          
          return calc_euclidian_distance(datapoint, cluster, sqr_sum, sqr_sum0);
      case CHEBYCHEV:
          return calc_chebychev_distance(datapoint, cluster);
      case COSINE:
	  return calc_cosine_distance(datapoint, cluster, sqr_sum, sqr_sum0);  
      case MANHATTAN:
          return calc_manhatten_distance(datapoint, cluster);
      case MANAHOLIS:
      case WEIGHTED_MANAHOLIS:
      case WEIGHTED:
      default:
          logstream(LOG_ERROR)<< "distance measure " << ac.distance_measure<< "  not implemented yet" << std::endl;
    }
    return -1;

}

/**
 *

v1 =
         0    1.0000   -3.5000         0   0
v2 =

         1    2        3               4   5
v3 =
         0    0.5000   0               4   0


 *
 *
 *
 * */
void test_distance(){
  sparse_vec v1;
  set_size(v1, 5);
  set_new(v1,1,1.0);
  set_new(v1,2,-3.5);

  sparse_vec v3;
  set_size(v3, 5);
  set_new(v3,1,0.5);
  set_new(v3,3,4);
  vec v2 = init_vec("1 2 3 4 5", 5);
  ac.distance_measure = EUCLIDEAN;
  double ret = calc_distance(v1, v2, sum_sqr(v2), sum_sqr(v1));
  assert(powf(ret - 9.233092656309694,2) < 1e-10);

  ret = calc_distance(v1, v1, sum_sqr(v1), sum_sqr(v1));
  assert(powf(ret - 0, 2) < 1e-10);

  ret = calc_distance(v1, v3, sum_sqr(v1), sum_sqr(v3));
  assert(powf(ret - 5.3385,2) <1e-8);

  ac.distance_measure = COSINE;
  ret = calc_distance(v1, v2, sum_sqr(v1), sum_sqr(v2));
  assert(powf(ret - 1.3149, 2)<1e-8);

  ret = calc_distance(v1, v3, sum_sqr(v1), sum_sqr(v3));
  assert(powf(ret - (1 - .0341), 2)<1e-8);


}

