#include "distance.h"
#include "clustering.h"
#include "../gabp/advanced_config.h"

extern advanced_config ac;
const char * distance_measure_name[] = {"EUCLIDEAN", "CHEBYCHEV", "MANAHOLIS", "MANHATTAN", "MINKOWSKI", "TANIMOTO", "WEIGTED", "WEIGHTED_MANAHOLIS", "COSINE"};

sparse_flt_dbl_vec minus(sparse_flt_dbl_vec & dvec1, sparse_flt_dbl_vec & dvec2);
vec minus(sparse_flt_dbl_vec & dvec1, vec & dvec2);
flt_dbl sum_sqr(sparse_flt_dbl_vec & dvec1);
sparse_flt_dbl_vec fabs(sparse_flt_dbl_vec& dvec1);
flt_dbl sum(sparse_flt_dbl_vec & dvec);



flt_dbl calc_tanimoto_distance( sparse_flt_dbl_vec & datapoint, sparse_flt_dbl_vec & cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0){ 
  flt_dbl a_mult_b = datapoint * cluster;
  return a_mult_b / (sqr_sum + sqr_sum0 - a_mult_b);
}

flt_dbl calc_tanimoto_distance( sparse_flt_dbl_vec & datapoint,  flt_dbl_vec &cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0){
  flt_dbl a_mult_b = datapoint * cluster;
  return a_mult_b / (sqr_sum + sqr_sum0 - a_mult_b);
}

flt_dbl calc_euclidian_distance( sparse_flt_dbl_vec & datapoint,  sparse_flt_dbl_vec &cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0){
  //sparse_flt_dbl_vec diff = minus(datapoint , cluster);
  //return sqrt(sum_sqr(diff));
  sparse_flt_dbl_vec mult = elem_mult(datapoint, cluster);
  flt_dbl diff = (sqr_sum + sqr_sum0 - 2*sum(mult));
  return sqrt(fabs(diff)); //because of numerical errors, diff may be negative
}


flt_dbl calc_euclidian_distance( sparse_flt_dbl_vec & datapoint,  flt_dbl_vec &cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0){
  flt_dbl dist = sqr_sum + sqr_sum0;
  //for (int i=0; i< datapoint.nnz(); i++){
  FOR_ITERATOR(i, datapoint){
      flt_dbl val = get_nz_data(datapoint, i);
      int pos = get_nz_index(datapoint, i);
      dist -= 2*val*cluster[pos];
   }
  if (fabs(dist) <1e-10)
     dist = 0;
  return sqrt(dist);
}

flt_dbl calc_chebychev_distance( sparse_flt_dbl_vec & datapoint,  sparse_flt_dbl_vec &cluster){
   sparse_flt_dbl_vec diff = minus(datapoint , cluster);
   flt_dbl ret = 0;
   FOR_ITERATOR(i, diff){
      ret = std::max(ret, (flt_dbl)fabs(get_nz_data(diff, i)));
   }
   return ret;

}
flt_dbl calc_chebychev_distance( sparse_flt_dbl_vec & datapoint,  flt_dbl_vec &cluster){
   flt_dbl_vec diff = minus(datapoint , cluster);
   flt_dbl ret = 0;
   for (int i=0; i< diff.size(); i++)
      ret = std::max(ret, (flt_dbl)fabs(diff[i]));

   return ret;

}

flt_dbl calc_manhatten_distance( sparse_flt_dbl_vec & datapoint,  sparse_flt_dbl_vec &cluster){
   sparse_flt_dbl_vec diff = minus(datapoint , cluster);
   sparse_flt_dbl_vec absvec = fabs(diff);
   flt_dbl ret = sum(absvec);
   return ret;

}
flt_dbl calc_manhatten_distance( sparse_flt_dbl_vec & datapoint,  flt_dbl_vec &cluster){
   flt_dbl_vec diff = minus(datapoint , cluster);
   flt_dbl ret = sum(fabs(diff));
   return ret;

}

flt_dbl calc_cosine_distance( sparse_flt_dbl_vec & datapoint,  sparse_flt_dbl_vec & cluster, flt_dbl sum_sqr, flt_dbl sum_sqr0){
   flt_dbl dotprod = dot_prod(datapoint,cluster);
   flt_dbl denominator = sqrt(sum_sqr0)*sqrt(sum_sqr);
   return 1.0 - dotprod / denominator; 
}

flt_dbl calc_cosine_distance( sparse_flt_dbl_vec & datapoint,  flt_dbl_vec & cluster, flt_dbl sum_sqr, flt_dbl sum_sqr0){
   flt_dbl dotprod = dot_prod(datapoint,cluster);
   flt_dbl denominator = sqrt(sum_sqr0)*sqrt(sum_sqr);
   return 1.0 - dotprod / denominator; 
}


flt_dbl calc_distance(sparse_flt_dbl_vec &datapoint,  sparse_flt_dbl_vec & cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0){
   switch(ac.distance_measure){
      case EUCLIDEAN:          
          return calc_euclidian_distance(datapoint, cluster, sqr_sum, sqr_sum0);
      case CHEBYCHEV:
          return calc_chebychev_distance(datapoint, cluster);
      case COSINE:
	  return calc_cosine_distance(datapoint, cluster, sqr_sum, sqr_sum0);  
      case MANHATTAN:
          return calc_manhatten_distance(datapoint, cluster);
      case TANIMOTO:
          return calc_tanimoto_distance(datapoint, cluster, sqr_sum , sqr_sum0);
       case MANAHOLIS:
      case WEIGHTED_MANAHOLIS:
      case WEIGHTED:
      default:
          logstream(LOG_ERROR)<< "distance measure " << ac.distance_measure<< "  not implemented yet" << std::endl;
    }
    return -1;

}


flt_dbl calc_distance(sparse_flt_dbl_vec &datapoint,  flt_dbl_vec & cluster, flt_dbl sqr_sum, flt_dbl sqr_sum0){
   switch(ac.distance_measure){
      case EUCLIDEAN:          
          return calc_euclidian_distance(datapoint, cluster, sqr_sum, sqr_sum0);
      case CHEBYCHEV:
          return calc_chebychev_distance(datapoint, cluster);
      case COSINE:
	  return calc_cosine_distance(datapoint, cluster, sqr_sum, sqr_sum0);  
      case MANHATTAN:
          return calc_manhatten_distance(datapoint, cluster);
      case TANIMOTO:
          return calc_tanimoto_distance(datapoint, cluster, sqr_sum , sqr_sum0);
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
  sparse_flt_dbl_vec v1;
  set_size(v1, 5);
  set_new(v1,1,1.0);
  set_new(v1,2,-3.5);

  sparse_flt_dbl_vec v3;
  set_size(v3, 5);
  set_new(v3,1,0.5);
  set_new(v3,3,4);
  flt_dbl_vec v2 = init_vec("1 2 3 4 5", 5);
  ac.distance_measure = EUCLIDEAN;
  flt_dbl ret = calc_distance(v1, v2, sum_sqr(v2), sum_sqr(v1));
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

  ac.distance_measure = TANIMOTO;
  sparse_flt_dbl_vec v4;
  set_size(v4,5);
  set_new(v4, 1, 1);
  set_new(v4, 3, 1);
  flt_dbl_vec v5 = init_vec("0 1 1 1 1", 5);
  ret = calc_distance(v4, v5, sum_sqr(v4), sum_sqr(v5));
  assert(powf(ret - 0.5,2) <1e-8);
  

}

