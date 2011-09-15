#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "clustering.h"

using namespace itpp;
extern problem_setup ps;

//assign sparse vector value v2 into v1
void assign(vec & v1, sparse_vec & v2){
  v1 = zeros(ps.N);
  for (int i=0; i< (v2).nnz(); i++){
     v1[v2.get_nz_index(i)] = v2.get_nz_data(i);
  }

}

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
  for (int i=0; i< (dvec).nnz(); i++){
     dmin = std::min(dmin, (dvec).get_nz_data(i));
  }
  return dmin;
}

double max( sparse_vec & dvec){
 
  double dmax = -1e100;
  for (int i=0; i< (dvec).nnz(); i++){
     dmax = std::max(dmax, (dvec).get_nz_data(i));
  }
  return dmax;
}

double sum( sparse_vec & dvec){
 
  double sum = 0;
  for (int i=0; i< (dvec).nnz(); i++){
     sum += (dvec).get_nz_data(i);
  }
  return sum;
}
double sum_sqr( sparse_vec & dvec){
 
  double sum = 0;
  for (int i=0; i< (dvec).nnz(); i++){
     sum += ((dvec).get_nz_data(i)*(dvec).get_nz_data(i));
  }
  return sum;
}

sparse_vec fabs( sparse_vec & dvec1){
   sparse_vec ret = dvec1;
   for (int i=0; i< ret.nnz(); i++){
       ret.set(ret.get_nz_index(i), fabs(ret.get_nz_data(i)));
   }
   return ret;
	
};

sparse_vec minus(sparse_vec &v1,sparse_vec &v2){
  sparse_vec ret(ps.N, v1.nnz() + v2.nnz());
  for (int i=0; i< v1.nnz(); i++){
      ret.add_elem(v1.get_nz_index(i), v1.get_nz_data(i) - get(v2, v1.get_nz_index(i)));
  }
  for (int i=0; i< v2.nnz(); i++){
      ret.add_elem(v2.get_nz_index(i), get(v1, v2.get_nz_index(i)) - v2.get_nz_data(i));
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
void plus( vec &v1,  sparse_vec &v2){
  for (int i=0; i< v2.nnz(); i++){
      v1[v2.get_nz_index(i)] += v2.get_nz_data(i);
  }
}
void plus_mul( vec &v1,  sparse_vec &v2, double factor){
  for (int i=0; i< v2.nnz(); i++){
      v1[v2.get_nz_index(i)] += factor*v2.get_nz_data(i);
  }
}
void minus( vec &v1, sparse_vec &v2){
  for (int i=0; i< v2.nnz(); i++){
      v1[v2.get_nz_index(i)] -= v2.get_nz_data(i);
  }
}





void test_math(){
   sparse_vec v1;
   sparse_vec v2;
   v1.add_elem(1,1.0);
   v2.add_elem(2,2.0);
   sparse_vec v3 = minus(v1, v2);
   assert(v3.get_nz_data(0) == 1.0);
   assert(v3.get_nz_data(1) == - 2.0);
   assert(v3.nnz() == 2);

   sparse_vec v4 = fabs(v3); 
   assert(v4.get_nz_data(0) == 1.0);
   assert(v4.get_nz_data(1) == 2.0);
   assert(v4.nnz() == 2);

   double sqr = sum_sqr(v4);
   assert(sqr == 5);
   assert(v4.nnz() == 2);

   double mmin = min(v4);
   assert(mmin = 1);
   double mmax = max(v4);
   assert(mmax = 2);

   vec v5 = vec("1 2 3 4");
   plus(v5, v1);
   assert(v5[0] == 1);
   assert(v5[1] == 3);
   assert(v5[2] == 5);
   assert(v5[3] == 4);
   assert(v5.size() == 4);

}
