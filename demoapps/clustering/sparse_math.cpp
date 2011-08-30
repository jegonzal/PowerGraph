#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "clustering.h"

using namespace itpp;
extern problem_setup ps;


double get(sparse_vec & v1, int pos){
  for (int i=0; i< v1.nnz(); i++){
     if (v1.get_nz_index(i) < pos)
	continue;
     else if (v1.get_nz_index(i) > pos)
	break;
     else if (v1.get_nz_index(i) == pos)
	return v1.get_nz_data(i);
  }
  return 0;
}


double min( sparse_vec & dvec){
 
  double dmin = 1e100;
  for (int i=0; i< ((sparse_vec)dvec).nnz(); i++){
     dmin = std::min(dmin, ((sparse_vec)dvec).get_nz_data(i));
  }
  return dmin;
}

double max( sparse_vec & dvec){
 
  double dmax = -1e100;
  for (int i=0; i< ((sparse_vec)dvec).nnz(); i++){
     dmax = std::max(dmax, ((sparse_vec)dvec).get_nz_data(i));
  }
  return dmax;
}

double sum( sparse_vec & dvec){
 
  double sum = 0;
  for (int i=0; i< ((sparse_vec)dvec).nnz(); i++){
     sum += ((sparse_vec)dvec).get_nz_data(i);
  }
  return sum;
}
double sum_sqr( sparse_vec & dvec){
 
  double sum = 0;
  for (int i=0; i< ((sparse_vec)dvec).nnz(); i++){
     sum += (((sparse_vec)dvec).get_nz_data(i)*((sparse_vec)dvec).get_nz_data(i));
  }
  return sum;
}

sparse_vec fabs( sparse_vec & dvec1){
   sparse_vec ret = dvec1;
   for (int i=0; i< ((sparse_vec)ret).nnz(); i++){
       ret.set(((sparse_vec)ret).get_nz_index(i), fabs(((sparse_vec)ret).get_nz_data(i)));
   }
   return ret;
	
};

sparse_vec minus(sparse_vec &v1,sparse_vec &v2){
  sparse_vec ret(ps.N, v1.nnz() + v2.nnz());
  for (int i=0; i< v1.nnz(); i++){
      ret.set_new(v1.get_nz_index(i), v1.get_nz_data(i) - get(v2, v1.get_nz_index(i)));
  }
  for (int i=0; i< v2.nnz(); i++){
      ret.set_new(v2.get_nz_index(i), get(v1, v2.get_nz_index(i)) - v2.get_nz_data(i));
  }
  return ret;
}
vec minus( sparse_vec &v1,  vec &v2){
  vec ret = zeros(v2.size());
  for (int i=0; i< v2.size(); i++){
      ret.set(i, get(v1, i) - v2[i]);
  }
  return ret;
}


