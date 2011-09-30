#include "clustering.h"

extern problem_setup ps;

//assign sparse vector value v2 into v1
void assign(vec & v1, sparse_vec & v2){
  v1 = zeros(ps.N);
  for (int i=0; i< (v2).nnz(); i++){
     v1[get_nz_index(vw, i)] = get_nz_data(v2, i);
  }

}

double get(sparse_vec & v1, int pos){
  FOR_ITERATOR(i, v1){
     if (get_nz_index(v1, i) < pos)
	continue;
     else if (get_nz_index(v1, i) > pos)
	break;
     else if (get_nz_index(v1, i) == pos)
	return get_nz_data(v1, i);
  }
  return 0;
}


double min( sparse_vec & dvec){
 
  double dmin = 1e100;
  FOR_ITERATOR(i, dvec){
     dmin = std::min(dmin, get_nz_data(dvec, i));
  }
  return dmin;
}

double max( sparse_vec & dvec){
 
  double dmax = -1e100;
  FOR_ITERATOR(i, dvec){
     dmax = std::max(dmax, get_nz_data(dvec, i));
  }
  return dmax;
}

double sum( sparse_vec & dvec){
 
  double sum = 0;
  FOR_ITERATOR(i, dvec){
     sum += get_nz_data(dvec, i);
  }
  return sum;
}
double sum_sqr( sparse_vec & dvec){
 
  double sum = 0;

  FOR_ITERATOR(i, dvec){
     sum += (get_nz_data(dvec, i)*get_nz_data(dvec, i));
  }
  return sum;
}

sparse_vec fabs( sparse_vec & dvec1){
   sparse_vec ret = dvec1;
   FOR_ITERATOR(i, ret){
       ret.set(get_nz_index(ret, i), fabs(get_nz_data(ret, i)));
   }
   return ret;
	
};

sparse_vec minus(sparse_vec &v1,sparse_vec &v2){
  sparse_vec ret(ps.N, v1.nnz() + v2.nnz());
  FOR_ITERATOR(i, v1){
      ret.add_elem(get_nz_index(v1, i), get_nz_data(v1, i) - get_val(v1, get_nz_index(v2, i)));
  }
  FOR_ITERATOR(i, v2){
      ret.add_elem(v2.get_nz_index(i), get_val(v1, get_nz_index(v2, i)) - get_nz_data(v2, i));
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
  FOR_ITERATOR(i, v2){ 
     v1[get_nz_index(v2, i)] += get_nz_data(v2, i);
  }
}
void plus_mul( vec &v1,  sparse_vec &v2, double factor){
  FOR_ITERATOR(i, v2){  
    v1[get_nz_index(v2, i)] += factor*get_nz_data(v2, i);
  }
}
void minus( vec &v1, sparse_vec &v2){
  FOR_ITERATOR(i, v2){ 
     v1[get_nz_index(v2, i)] -= get_nz_data(v2, i);
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
