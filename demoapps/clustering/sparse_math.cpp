#include "clustering.h"
#include <fstream>

extern problem_setup ps;

//assign sparse vector value v2 into v1
void assign(vec & v1, sparse_vec & v2){
  v1 = zeros(ps.N);
  FOR_ITERATOR(i, v2){
     v1[get_nz_index(v2, i)] = get_nz_data(v2, i);
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
void plus_mul( vec &v1,  sparse_vec &v2, double factor){
  FOR_ITERATOR(i, v2){  
    v1[get_nz_index(v2, i)] += factor*get_nz_data(v2, i);
  }
}

double sum( sparse_vec & dvec){
  double sum = 0;
  FOR_ITERATOR(i, dvec){
     sum += get_nz_data(dvec, i);
  }
  return sum;
}




void test_math(){
   sparse_vec v1(4);
   sparse_vec v2(4);
   set_new(v1,1,1.0);
   set_new(v2,2,2.0);
   sparse_vec v3 = minus(v1, v2);
   assert(get_nz_data(v3,0) == 1.0);
   assert(get_nz_data(v3,1) == - 2.0);
   assert(nnz(v3)== 2);
   assert(get_val(v1, 1) == 1.0);
   assert(get_val(v2, 2) == 2.0);

   assert(sum(v1) == 1.0);
   assert(sum(v2) == 2.0);
   assert(sum(v3) == -1.0);

   sparse_vec v4 = fabs(v3); 
   assert(nnz(v4) == 2);
   assert(get_nz_data(v4,0) == 1.0);
   assert(get_nz_data(v4,1) == 2.0);

   double sqr = sum_sqr(v4);
   assert(sqr == 5);
   assert(nnz(v4) == 2);

   double mmin = min(v4);
   assert(mmin = 1);
   double mmax = max(v4);
   assert(mmax = 2);
   vec v5 = init_vec("1 2 3 4",4);
   plus(v5, v1);
   assert(v5[0] == 1);
   assert(v5[1] == 3);
   assert(v5[2] == 3);
   assert(v5[3] == 4);
   assert(v5.size() == 4);
   minus(v5, v1);
   assert(v5[0] == 1);
   assert(v5[1] == 2);
   assert(v5[2] == 3);
   assert(v5[3] == 4);
   assert(v5.size() == 4);
 
   vec vcum = cumsum(v5);
   assert(vcum[0] == 1);
   assert(vcum[1] == 3);
   assert(vcum[2] == 6);
   assert(vcum[3] == 10);
    
   double dot = dot_prod(v5, v1);
   assert(dot == 2);

   dot = dot_prod(v5,v5);
   assert(30 == dot);

   dot = dot_prod(v1, v1);
   assert(dot == 1);

   vec powv = pow(v5, 2);
   assert(powv[0] == 1);
   assert(powv[1] == 4);
   assert(powv[2] == 9);
   assert(powv[3] == 16);
   assert(get_nz_data(v1, 0) == 1.0);
   set_new(v1, 18, 3.0);
   assert(get_nz_data(v1,1) == 3.0);
   //set_size(v1, 19);
 
   assert(get_val(v5,0) == 1.0);
   assert(get_val(v1,0) == 0); 
   assert(get_val(v1,1) == 1.0); 
   assert(get_val(v1,18) == 3.0); 

   sparse_vec v6;
   set_new(v6, 3, 4.2);
   FOR_ITERATOR(i, v6){
   set_div(v6, i, 2.0);
   }
   assert(get_nz_data(v6, 0) == 2.1);
   assert(get_val(v6, 3) == 2.1);
   assert(nnz(v6) == 1);

   mat mymat = init_mat("1 2 3; 3 2 1; 1 2 3", 3, 3);

   remove("stam");
   it_file saved("stam");
   saved << Name("vector");
   saved << v5;
   saved << Name("matrix");
   saved << mymat;
   saved.close();

   vec v5_saved;
   it_file loaded("stam");
   loaded >> Name("vector");
   loaded >> v5_saved;
   assert(v5_saved.size() == 4);
   assert(v5_saved[0] == 1);
   assert(v5_saved[1] == 2);
   assert(v5_saved[2] == 3);
   assert(v5_saved[3] == 4);
   loaded >> Name("matrix");
   mat mymat2;
   loaded >> mymat2;
   assert(mymat2.rows() == 3 && mymat2.cols() == 3);
   assert(get_val(mymat2, 0, 0) == 1.0);
   assert(get_val(mymat2, 2, 2) == 3.0);


   vec start = head(v5, 2);
   assert(start.size() == 2);
   assert(start[0] == 1);
   assert(start[1] == 2);

   ivec delvec(4);
   delvec[0] = 1;
   delvec[1] = 2;
   delvec[2] = 3;  
   delvec[3] = 4;
   del(delvec,3);
   assert(delvec.size() == 3);
   assert(delvec[0] == 1);
   assert(delvec[1] == 2);
   assert(delvec[2] == 3);
   del(delvec,1);
   assert(delvec.size() == 2);
   assert(delvec[0] == 1);
   assert(delvec[1] == 3);

   ivec delvec2(4);
   delvec2[0] = 1;
   delvec2[1] = 2;
   delvec2[2] = 3;
   delvec2[3] = 4;
   del(delvec2, 0);
   assert(delvec2.size() == 3);
   assert(delvec2[0] = 2);
}
