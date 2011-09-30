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

double calc_euclidian_distance( sparse_vec & datapoint,  sparse_vec &cluster, double sqr_sum){
  sparse_vec diff = minus(datapoint , cluster);
  return sqrt(sum_sqr(diff));
}


double calc_euclidian_distance( sparse_vec & datapoint,  vec &cluster, double sqr_sum){
  double dist = sqr_sum;
  //for (int i=0; i< datapoint.nnz(); i++){
  FOR_ITERATOR(i, datapoint){
      double val = get_nz_data(datapoint, i);
      int pos = get_nz_index(datapoint, i);
      dist += (((val - cluster[pos])*(val - cluster[pos])) - cluster[pos]*cluster[pos]);
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

double calc_cosine_distance( sparse_vec & datapoint,  sparse_vec & cluster){
   double len_sqr1 = sum_sqr(datapoint);
   double len_sqr2 = sum_sqr(cluster);
   double dotprod = datapoint*cluster;
   double denominator = sqrt(len_sqr1)*sqrt(len_sqr2);
   return 1.0 - dotprod / denominator; 
}

double calc_cosine_distance( sparse_vec & datapoint,  vec & cluster){
   double len_sqr1 = sum_sqr(datapoint);
   double len_sqr2 = sum_sqr(cluster);
   double dotprod = datapoint*cluster;
   double denominator = sqrt(len_sqr1)*sqrt(len_sqr2);
   return 1.0 - dotprod / denominator; 
}


double calc_distance(sparse_vec &datapoint,  vec & cluster, double sqr_sum){
   switch(ac.distance_measure){
      case EUCLIDEAN:          
          return calc_euclidian_distance(datapoint, cluster, sqr_sum);
      case CHEBYCHEV:
          return calc_chebychev_distance(datapoint, cluster);
      case COSINE:
	  return calc_cosine_distance(datapoint, cluster);  
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


void test_distance(){
  sparse_vec v1;
  v1.add_elem(1,1.0);
  v1.add_elem(2,-3.5);
  vec v2("1 2 3 4 5");
  ac.distance_measure = EUCLIDEAN;
  double ret = calc_distance(v1, v2);
  assert(powf(ret - 9.233092656309694,2) < 1e-10);

}

