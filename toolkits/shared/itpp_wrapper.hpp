/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://graphlab.org
 *
 * Code by Danny Bickson, CMU
 *
 */

#ifndef ITPP_WRAPPER
#define ITPP_WRAPPER

/***
 *
 *  SET OF WRAPPER FUNCTIONS TO BE USED WITH IT++
 *
 * */

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
using namespace itpp;

inline void debug_print_vec(const char * name,const vec& _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec[i] == 0)
      printf("      0    ");
    else printf("%12.4g    ", _vec[i]);
  printf("\n");
}
inline void debug_print_vec(const char * name,const double* _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec[i] == 0)
      printf("      0    ");
    else printf("%12.4g    ", _vec[i]);
  printf("\n");
}
inline void compact(sparse_vec & a){
  //TODO
}
inline void set_val(mat& A, int row, int col, double val){
  A.set(row, col, val);
}
inline void set_val(imat& A, int row, int col, int val){
  A.set(row, col, val);
}
inline double get_val(const mat& A, int row, int col){
  return A.get(row, col);
}
inline int get_val(const imat& A, int row, int col){
  return A.get(row, col);
}
inline vec get_col(mat& A, int col){
  return A.get_col(col);
}
inline vec get_row(mat& A, int row){
  return A.get_row(row);
}
inline void set_col(mat& A, int col, const vec & val){
  A.set_col(col, val);
}
inline void set_row(mat& A, int row, const vec & val){
  A.set_row(row, val);
}
inline void set_diag(mat &A, vec &v){
  A = diag(v);
}
inline vec init_vec(const char * string, int size){
  return vec(string);
}
inline mat init_mat(const char * string, int row, int col){
  return mat(string);
}
inline vec init_dbl_vec(const char * string, int size){
  return vec(string);
}
inline vec head(const vec &v, int num){
  return v.mid(0,num);
}
inline vec mid(const vec&v, int start, int num){
  return v.mid(start, std::min(num, v.size() - start));
}
inline vec tail(const vec&v, int num){
  return v.mid(v.size()-num, num);
}
inline ivec head(const ivec &v, int num){
  return v.mid(0,num);
}
inline ivec sort(ivec &a){
   Sort<int> sorter;
   sorter.sort(0, a.size()-1, a);
   return a;
}
inline void del(ivec&a, int i){
  a.del(i);
}
inline ivec sort_index(vec& a){
  Sort<double> sorter;
  return sorter.sort_index(0, a.size()-1, a);
}
inline void set_val(vec & v, int pos, double val){
  v.set(pos, val);
}
inline mat get_cols(const mat &A, const ivec & ind){
  return A.get_cols(ind);
}
inline const double * data(const mat &A){
  return A._data();
}
inline const int * data(const imat &A){
  return A._data();
}
inline const double * data(const vec &v){
  return v._data();
}
inline void set_size(sparse_vec &v, int size){
  v.set_size(size);
}
inline void set_size(mat &a, int row, int col){
  a.set_size(row, col);
}
inline void set_new(sparse_vec&v, int ind, double val){
  v.set_new(ind, val);
} 
#define FOR_ITERATOR(i, v) \
    for (int i = 0; i < v.nnz(); i++)
inline int get_nz_index(sparse_vec &v, int i){
  return v.get_nz_index(i);
}
inline double get_nz_data(sparse_vec &v, int i){
  return v.get_nz_data(i);
}
inline int nnz(sparse_vec & v){
  return v.nnz();
}
inline double dot_prod(sparse_vec &v1, sparse_vec & v2){
  return v1*v2;
}
inline double dot_prod(vec &v1, vec & v2){
  return v1*v2;
}
inline double dot_prod(const sparse_vec &v1, const vec & v2){
  return v1*v2;
}
inline double dot_prod(vec &v1, sparse_vec & v2){
  return v1*v2;
}
inline double get_val(sparse_vec & v1, int i){
   FOR_ITERATOR(j, v1){
      if (v1.get_nz_index(j) == i)
         return v1.get_nz_data(j);
   }
   return 0;
}
inline double get_val(vec & v1, int i){
  return v1[i];
}
inline void set_div(sparse_vec&v, int i, double val){
  v.set(v.get_nz_index(i) ,v.get_nz_data(i) / val);
}
inline sparse_vec minus(sparse_vec &v1,sparse_vec &v2){
/*  sparse_vec ret; 
  for (int i=0; i< v1.nnz(); i++){
      ret.set_new(v1.get_nz_index(i), v1.get_nz_data(i) - get_val(v2, v1.get_nz_index(i)));
  }
  for (int i=0; i< v2.nnz(); i++){
      ret.set_new(v2.get_nz_index(i), get_val(v1, v2.get_nz_index(i)) - v2.get_nz_data(i));
  }
  return ret;*/
  return v1+(-v2);
}
inline vec minus( sparse_vec &v1,  vec &v2){
  vec ret = -v2;;
  FOR_ITERATOR(i, v1){  
    ret.set(v1.get_nz_index(i), ret.get(v1.get_nz_index(i)) + v1.get_nz_data(i));
  }
  return ret;
}
inline void plus( vec &v1,  sparse_vec &v2){
  FOR_ITERATOR(i, v2){ 
     v1[get_nz_index(v2, i)] += get_nz_data(v2, i);
  }
}
inline void minus( vec &v1, sparse_vec &v2){
  FOR_ITERATOR(i, v2){ 
     v1[get_nz_index(v2, i)] -= get_nz_data(v2, i);
  }
}
inline sparse_vec fabs( sparse_vec & dvec1){
   sparse_vec ret(dvec1.size(), dvec1.nnz());
   FOR_ITERATOR(i, dvec1){
       set_new(ret,get_nz_index(dvec1, i), fabs(get_nz_data(dvec1, i)));
   }
   return ret;
	
};

inline vec fabs(const vec & a){
   return abs(a);
}
inline double abs_sum(const mat& A){
   return sumsum(abs(A));
}
inline double abs_sum(const vec&v){
  return sum(fabs(v));
}
inline bool eig_sym(const mat & T, vec & eigenvalues, mat & eigenvectors){
  itpp::eig_sym(T,eigenvalues, eigenvectors);
  eigenvalues = reverse(eigenvalues);
  mat reverse_cols = zeros(eigenvectors.rows(), eigenvectors.cols());
  for (int i=0; i< eigenvectors.cols(); i++){
    reverse_cols.set_col(i, eigenvectors.get_col(eigenvectors.rows() - i - 1));
  }
  eigenvectors = reverse_cols;
  return true;
}

inline double sum_sqr(sparse_vec & v){
  double sum = 0;
  FOR_ITERATOR(i, v){
     sum+= powf(v.get_nz_data(i),2);
  }
  return sum;
}
inline void dot2(const vec&  x1, const vec& x3, mat & Q, int j, int len){
	for (int i=0; i< len; i++){
		Q.set(i,j,(x1[i] * x3[i]));
	}
}
inline void assign(vec & v1, sparse_vec & v2, int N){
  v1 = zeros(N);
  FOR_ITERATOR(i, v2){
     v1[get_nz_index(v2, i)] = get_nz_data(v2, i);
  }

}

inline double get(sparse_vec & v1, int pos){
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


inline double min( sparse_vec & dvec){
 
  double dmin = 1e100;
  FOR_ITERATOR(i, dvec){
     dmin = std::min(dmin, get_nz_data(dvec, i));
  }
  return dmin;
}

inline double max( sparse_vec & dvec){
 
  double dmax = -1e100;
  FOR_ITERATOR(i, dvec){
     dmax = std::max(dmax, get_nz_data(dvec, i));
  }
  return dmax;
}
inline void plus_mul( vec &v1,  sparse_vec &v2, double factor){
  FOR_ITERATOR(i, v2){  
    v1[get_nz_index(v2, i)] += factor*get_nz_data(v2, i);
  }
}

inline double sum( sparse_vec & dvec){
  double sum = 0;
  FOR_ITERATOR(i, dvec){
     sum += get_nz_data(dvec, i);
  }
  return sum;
}
inline vec dbl_fzeros(int size){ 
  return itpp::zeros(size);
}
inline mat dbl_fzeros(int rows, int cols){
  return itpp::zeros(rows, cols);
}
inline void print(sparse_vec & vec){
  int cnt = 0;
  FOR_ITERATOR(i, vec){
    std::cout<<get_nz_index(vec, i)<<":"<< get_nz_data(vec, i) << " ";
    cnt++;
    if (cnt >= 20)
       break;
  }
  std::cout<<std::endl;
}
inline vec init_vec(double * array, int size){
  return vec(array, size);
}
/**
 * It seems that it++ random number generator is not thread safe so
 * we are using graphlab's
 */
inline ivec randi(int size, int from, int to){
  ivec ret(size);
  for (int i=0; i<size; i++)
    ret[i]= graphlab::random::uniform<int>(from,to);
  return ret;
}
inline int randi(int from, int to){
  return graphlab::random::uniform<int>(from,to);
}
inline vec sqrt(const vec & v){
  return itpp::sqrt(v);
}
#endif //ITPP_WRAPPER
