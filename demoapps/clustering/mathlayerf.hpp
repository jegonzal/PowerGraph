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
 *      http://www.graphlab.ml.cmu.edu
 *  
 *  Code written by Danny Bickson, CMU
 */

#ifndef MATH_LAYER_FLOAT_GRAPHLAB
#define MATH_LAYER_FLOAT_GRAPHLAB




/**
 *
 * SET OF WRAPPER FUNCTIONS TO ALLOW USING EIGEN
 */
#ifdef HAS_EIGEN

#include <iostream>
#include <fstream>
#include <ostream>

#include "Eigen/Dense"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include "Eigen/Sparse"
#include "Eigen/Cholesky"
#include "Eigen/Eigenvalues"


using namespace Eigen;

typedef MatrixXf fmat;
typedef VectorXf fvec;
typedef SparseVector<float> sparse_fvec;

/*
inline mat eye(int size){
  return mat::Identity(size, size);
}
inline vec ones(int size){
  return vec::Ones(size);
}*/

inline void sort(fvec & a){
   std::sort(a.data(), a.data()+a.size());
}
inline ivec sort_index(const fvec &a){
  ivec ret(a.size()); 
  std::vector<std::pair<float,int> > D;
  // 	
  D.reserve(a.size());
  for (int i=0;i<a.size();i++)
    D.push_back(std::make_pair<float,int>(a.coeff(i),i));
  std::sort(D.begin(),D.end());
  for (int i=0;i<a.size();i++)
  { 
    ret[i]=D[i].second;
  } 
  return ret;
}


inline fvec vec2fvec(const vec & a){
  fvec ret(a.size());
  for (int i=0; i< a.size(); i++)
   ret[i] = (float)a[i];
  return ret;
}
inline vec fvec2vec(const fvec & a){
  vec ret(a.size());
  for (int i=0; i< a.size(); i++)
   ret[i] = (double)a[i];
  return ret;
}
inline mat fmat2mat(const fmat  & a){
  mat ret(a.rows(), a.cols());
  for (int i=0; i< a.rows(); i++)
    for (int j=0; j< a.cols(); j++)
      ret(i,j) = a(i,j);
  return ret;
}
inline fmat mat2fmat(const mat  & a){
  fmat ret(a.rows(), a.cols());
  for (int i=0; i< a.rows(); i++)
    for (int j=0; j< a.cols(); j++)
      ret(i,j) = a(i,j);
  return ret;
}
inline void compact(sparse_fvec & v){
   //v.conservativeResize(v.nonZeros());
   //TODO
}
inline fvec init_fvec(const float * array, int size){
  fvec ret(size);
  memcpy(ret.data(), array, size*sizeof(float));
  return ret;
}
inline fvec fzeros(int size){
  return fvec::Zero(size);
}
inline fmat fzeros(int rows, int cols){
  return fmat::Zero(rows, cols);
}
inline void debug_print_vec(const char * name,const fvec& _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec(i) == 0)
      printf("      0    ");
    else printf("%12.4f    ", _vec(i));
  printf("\n");
}

inline void dot2(const fvec&  x1, const fvec& x3, fmat & Q, int j, int len){
	for (int i=0; i< len; i++){
		Q(i,j) = (x1(i) * x3(i));
	}
}

inline bool ls_solve_chol(const fmat &A, const fvec &b, fvec &result){
    //result = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    result = A.ldlt().solve(b);
    return true;
}
inline bool ls_solve(const fmat &A, const fvec &b, fvec &result){
    //result = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    result = A.ldlt().solve(b);
    return true;
}
inline bool chol(fmat& sigma, fmat& out){
   out = sigma.llt().matrixLLT();
   return true;
}
inline bool backslash(const fmat& A, const fvec & b, fvec & x){
   x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
   return true;
} 
inline fmat transpose(fmat & A){
   return A.transpose();
}
inline void set_val(fmat &A, int row, int col, float val){
  A(row, col) = val;
}
inline float get_val(const fmat &A, int row, int col){
  return A(row, col);
}
inline fvec get_col(const fmat& A, int col){
  return A.col(col);
}
inline fvec get_row(const fmat& A, int row){
  return A.row(row);
}
inline void set_col(fmat& A, int col, const fvec & val){
  A.col(col) = val;
}
inline void set_row(fmat& A, int row, const fvec & val){
  A.row(row) = val;
}
inline void set_diag(fmat &A, fvec & v){
   A.diagonal()=v;
}
inline float sumsum(const fmat & A){
   return A.sum();
}
inline fmat init_fmat(const char * string, int row, int col){
  fmat out(row, col);
  char buf[2056];
  strcpy(buf, string);
  char *pch = strtok(buf," \r\n\t;");
  for (int i=0; i< row; i++){
    for (int j=0; j< col; j++){
     out(i,j) = atof(pch);
     pch = strtok (NULL, " \r\n\t;");
    }
  }
  return out;
}
inline fvec init_fvec(const char * string, int size){
  fvec out(size);
  char buf[2056];
  strcpy(buf, string);
  char *pch = strtok (buf," \r\n\t;");
  int i=0;
  while (pch != NULL)
  {
     out(i) =atof(pch);
     pch = strtok (NULL, " \r\n\t;");
     i++;
  }
  assert(i == size);
  return out;
}
inline float norm(const fmat &A, int pow=2){
     return A.squaredNorm();
}
inline fmat inv(const fmat&A){
   return A.inverse();
}
inline bool inv(const fmat&A, fmat &out){
   out = A.inverse();
   return true;
}
inline fmat outer_product(const fvec&a, const fvec&b){
   return a*b.transpose();
}
inline fvec head(const fvec& v, int num){
   return v.head(num);
}
inline fvec mid(const fvec&v, int start, int num){
   return v.segment(start, std::min(num, (int)(v.size()-start)));
}
inline fvec tail(const fvec&v,  int num){
   return v.segment(v.size() - num, num);
}
inline fvec elem_mult(const fvec&a, const fvec&b){
   fvec ret = a;
   for (int i=0; i<b.size(); i++)
      ret(i) *= b(i);
   return ret;
}
inline sparse_fvec elem_mult(const sparse_fvec&a, const sparse_fvec&b){
   return a.cwiseProduct(b);
}
inline float sum(const fvec & a){
  return a.sum();
}
template<>
inline double sum_sqr<fvec>(const fvec & a){
  fvec ret = a.array().pow(2);
  return ret.sum();
}
inline float trace(const fmat & a){
  return a.trace();
}
/*inline float min(const fvec &a){
  return a.minCoeff();
}
inline float max(const fvec & a){
  return a.maxCoeff();
}*/
inline fvec frandu(int size){
  return fvec::Random(size);
}
inline float frandu(){
  return fvec::Random(1)(0);
}
inline fmat get_cols(const fmat&A, ivec & cols){
  fmat a(A.rows(), cols.size());
  for (int i=0; i< cols.size(); i++)
    set_col(a, i, get_col(A, cols[i]));
  return a;
}
inline void set_val(fvec & v, int pos, float val){
  v(pos) = val;
}
inline float dot(const fvec&a, const fvec& b){
   return a.dot(b);
}
inline fvec reverse(fvec& a){
   return a.reverse();
}
inline const float * data(const fmat &A){
  return A.data();
}
inline const float * data(const fvec &v){
  return v.data();
}
inline void set_size(sparse_fvec &v, int size){
  //did not find a way to declare vector dimension, yet
}
inline void set_new(sparse_fvec&v, int ind, float val){
  v.insert(ind) = val;
} 
inline int nnz(sparse_fvec& v){
  return v.nonZeros();
}
inline int get_nz_index(sparse_fvec &v, sparse_fvec::InnerIterator& i){
  return i.index();
}
inline float get_nz_data(sparse_fvec &v, sparse_fvec::InnerIterator& i){
  return i.value();
}
#define FOR_ITERATOR2(i,v) \
  for (sparse_fvec::InnerIterator i(v); i; ++i)

template<>
inline double sum_sqr<sparse_fvec>(const sparse_fvec & a){
  double sum=0;
  FOR_ITERATOR2(i,a){
    sum+= powf(i.value(),2);
  }
  return sum;
}

inline float get_nz_data(sparse_fvec &v, int i){
  assert(nnz(v) > i);
  int cnt=0;
  FOR_ITERATOR2(j, v){
    if (cnt == i){
      return j.value();
    }
    cnt++;
  }
  return 0.0;
}



inline void print(sparse_fvec & vec){
  int cnt = 0;
  FOR_ITERATOR2(i, vec){
    std::cout<<get_nz_index(vec, i)<<":"<< get_nz_data(vec, i) << " ";
    cnt++;
    if (cnt >= 20)
       break;
  }
  std::cout<<std::endl;
}

inline fvec pow(const fvec&v, int exponent){
  fvec ret = fvec(v.size());
  for (int i=0; i< v.size(); i++)
    ret[i] = powf(v[i], exponent);
  return ret;
}
inline float dot_prod(sparse_fvec &v1, sparse_fvec & v2){
  return v1.dot(v2);
}
inline float dot_prod(const fvec &v1, const fvec & v2){
  return v1.dot(v2);
}
inline float dot_prod(sparse_fvec &v1, const fvec & v2){
  float sum = 0;
  for (int i=0; i< v2.size(); i++){
    sum+= v2[i] * v1.coeffRef(i);
  }
  return sum;
}
inline fvec cumsum(fvec& v){
  fvec ret = v;
  for (int i=1; i< v.size(); i++)
     for (int j=0; j< i; j++)
       ret(i) += v(j);
  return ret;
}
inline float get_val(sparse_fvec & v1, int i){ //TODO optimize performance
  for (sparse_fvec::InnerIterator it(v1); it; ++it)
    if (it.index() == i)
       return it.value();

  return 0;
} 
inline float get_val(fvec & v1, int i){
  return v1(i);
}
inline void set_div(sparse_fvec&v, sparse_fvec::InnerIterator i, float val){
   v.coeffRef(i.index()) /= val;
}
inline sparse_fvec minus(sparse_fvec &v1,sparse_fvec &v2){
   return v1-v2;
}
inline fvec minus( sparse_fvec &v1,  fvec &v2){
   return v1-sparse_fvec(v2);
}
inline void plus( fvec &v1,  sparse_fvec &v2){
   FOR_ITERATOR2(i, v2){
     v1[i.index()] += i.value();
   }
}
inline void plus_mul( fvec &v1,  sparse_fvec &v2, float factor){
  FOR_ITERATOR2(i, v2){  
    v1[get_nz_index(v2, i)] += factor*get_nz_data(v2, i);
  }
}


inline void minus( fvec &v1, sparse_fvec &v2){
   FOR_ITERATOR2(i, v2){
      v1[i.index()] -= i.value();
   }
}
inline sparse_fvec fabs( sparse_fvec & dvec1){
   sparse_fvec ret = dvec1;
   FOR_ITERATOR2(i, ret){
      ret.coeffRef(i.index()) = fabs(i.value()); 
   }	
   return ret;
};

inline fvec fabs( const fvec & dvec1){
   fvec ret(dvec1.size());
   for (int i=0; i< dvec1.size(); i++){
      ret(i) = fabs(dvec1(i));
   }	
   return ret;
};
inline float abs_sum(const fmat& A){
  float sum =0;
  for (int i=0; i< A.rows(); i++)
    for (int j=0; j< A.cols(); j++)
      sum += fabs(A(i,j));
  return sum;
}
inline float abs_sum(const fvec &v){
  float sum =0;
  for (int i=0; i< v.size(); i++)
      sum += fabs(v(i));
  return sum;
}
inline float sum(const sparse_fvec &v){
  float  sum =0;
  FOR_ITERATOR2(i, v){
      sum += i.value();
  }
  return sum;
}
inline fvec sqrt(fvec & v){
   fvec ret(v.size());
   for (int i=0; i< v.size(); i++){
      ret[i] = sqrt(v(i));
   }
   return ret;
}
inline void assign(fvec & v1, sparse_fvec & v2, int N){
  v1 = fzeros(N);
  FOR_ITERATOR2(i, v2){
     v1[get_nz_index(v2, i)] = get_nz_data(v2, i);
  }

}

#else //eigen is not found
/***
 *
 *  SET OF WRAPPER FUNCTIONS TO BE USED WITH IT++
 *
 * */
#if defined(HAS_ITPP)

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
using namespace itpp;

//#undef sparse_fvec
typedef Sparse_Vec<float> sparse_fvec;
typedef Vec<float> fvec;
typedef Mat<float> fmat;
inline void compact(sparse_fvec &v){
   v.compact();
}
inline fvec fzeros(int size){
  fvec ret(size);
  for (int i=0; i< size; i++)
    ret[i] = 0;
  return ret;
}
inline fmat fzeros(int rows, int cols){
  fmat ret(rows, cols);
  for (int i=0; i< rows; i++)
    for (int j=0; j< cols; j++)
       ret.set(i,j,0);

  return ret;
}
inline void set_val(fmat& A, int row, int col, float val){
  A.set(row, col, val);
}
inline float get_val(const fmat& A, int row, int col){
  return A.get(row, col);
}
inline float get_val(const fvec & v, int pos){
  return v[pos];
}
inline fvec get_col(fmat& A, int col){
  return A.get_col(col);
}
inline fvec get_row(fmat& A, int row){
  return A.get_row(row);
}
inline void set_col(fmat& A, int col, const fvec & val){
  A.set_col(col, val);
}
inline void set_row(fmat& A, int row, const fvec & val){
  A.set_row(row, val);
}
inline void set_diag(fmat &A, fvec &v){
  A = diag(v);
}
inline fvec init_fvec(const char * string, int size){
  fvec out(size);
  char buf[2056];
  strcpy(buf, string);
  char *pch = strtok (buf," \r\n\t;");
  int i=0;
  while (pch != NULL)
  {
     out[i] =(float)atof(pch);
     pch = strtok (NULL, " \r\n\t;");
     i++;
  }
  assert(i == size);
  return out;
}
inline fmat init_fmat(const char * string, int row, int col){
  fmat out(row, col);
  char buf[2056];
  strcpy(buf, string);
  char *pch = strtok(buf," \r\n\t;");
  for (int i=0; i< row; i++){
    for (int j=0; j< col; j++){
     out.set(i,j,(float) atof(pch));
     pch = strtok (NULL, " \r\n\t;");
    }
  }
  return out;
}
inline fvec head(const fvec &v, int num){
  return v.mid(0,num);
}
inline fvec mid(const fvec&v, int start, int num){
  return v.mid(start, std::min(v.size(), num));
}
inline fvec tail(const fvec&v, int num){
  return v.mid(v.size()-num, num);
}
inline void set_val(fvec & v, int pos, float val){
  v.set(pos, val);
}
inline fmat get_cols(const fmat &A, const ivec & ind){
  return A.get_cols(ind);
}
inline const float * data(const fmat &A){
  return A._data();
}
inline const float * data(const fvec &v){
  return v._data();
}
inline void set_size(sparse_fvec &v, int size){
  v.set_size(size);
}
inline void set_new(sparse_fvec&v, int ind, float val){
  v.set_new(ind, val);
} 
inline int get_nz_index(const sparse_fvec &v, int i){
  return v.get_nz_index(i);
}
inline float get_nz_data(sparse_fvec &v, int i){
  return v.get_nz_data(i);
}
inline int nnz(sparse_fvec & v){
  return v.nnz();
}
inline float dot_prod(sparse_fvec &v1, sparse_fvec & v2){
  return itpp::operator*(v1,v2);
}
inline float dot_prod(const sparse_fvec &v1, const fvec & v2){
  return itpp::operator*(v1,v2);
}
inline float dot_prod(fvec &v1, sparse_fvec & v2){
  return itpp::operator*(v1,v2);
}
inline float get_val(sparse_fvec & v1, int i){
   FOR_ITERATOR(j, v1){
      if (v1.get_nz_index(j) == i)
         return v1.get_nz_data(j);
   }
   return 0;
}

inline void set_div(sparse_fvec&v, int i, float val){
  v.set(v.get_nz_index(i) ,v.get_nz_data(i) / val);
}
inline sparse_fvec minus(sparse_fvec &v1,sparse_fvec &v2){
/*  sparse_fvec ret; 
  for (int i=0; i< v1.nnz(); i++){
      ret.set_new(v1.get_nz_index(i), v1.get_nz_data(i) - get_val(v2, v1.get_nz_index(i)));
  }
  for (int i=0; i< v2.nnz(); i++){
      ret.set_new(v2.get_nz_index(i), get_val(v1, v2.get_nz_index(i)) - v2.get_nz_data(i));
  }
  return ret;*/
  return v1+(-v2);
}
inline fvec minus( sparse_fvec &v1,  fvec &v2){
  fvec ret = -v2;;
  FOR_ITERATOR(i, v1){  
    ret.set(v1.get_nz_index(i), ret.get(v1.get_nz_index(i)) + v1.get_nz_data(i));
  }
  return ret;
}
inline void plus( fvec &v1,  sparse_fvec &v2){
  FOR_ITERATOR(i, v2){ 
     v1[get_nz_index(v2, i)] += get_nz_data(v2, i);
  }
}
inline void minus( fvec &v1, sparse_fvec &v2){
  FOR_ITERATOR(i, v2){ 
     v1[get_nz_index(v2, i)] -= get_nz_data(v2, i);
  }
}
inline sparse_fvec fabs( sparse_fvec & dvec1){
   sparse_fvec ret(dvec1.size(), dvec1.nnz());
   FOR_ITERATOR(i, dvec1){
       set_new(ret,get_nz_index(dvec1, i), fabs(get_nz_data(dvec1, i)));
   }
   return ret;
	
};

inline fvec fabs(const fvec & a){
   fvec ret = fzeros(a.size());
   for (int i=0; i< a.size(); i++)
     ret[i] = fabs(a[i]);
   return ret;
}
inline fmat abs(const fmat & a){
   fmat ret(a.rows(), a.cols());
   for (int i=0; i< a.rows(); i++)
     for (int j=0; j< a.cols(); j++)
       ret.set(i,j,fabs(a.get(i,j)));
   return ret;
}
inline float abs_isum(const fmat& A){
   return sumsum(abs(A));
}
inline float abs_sum(const fvec&v){
  return sum(fabs(v));
}


/*inline float sum_sqr(sparse_fvec & v){
  float sum = 0;
  FOR_ITERATOR(i, v){
     sum+= powf(v.get_nz_data(i),2);
  }
  return sum;
}*/
inline float sum_sqr(sparse_fvec & v){
  double sum = 0;
  FOR_ITERATOR(i, v){
     sum+= powf(v.get_nz_data(i),2);
  }
  return sum;
}

inline fvec pow(const fvec &v, float exponent){
  fvec ret(v.size());
  for (int i=0; i< v.size(); i++)
    ret[i] = powf(v[i], exponent);
  return ret;
}
inline fvec vec2fvec(const vec & a){
  fvec ret(a.size());
  for (int i=0; i< a.size(); i++)
   ret[i] = (float)a[i];
  return ret;
}
inline vec fvec2vec(const fvec & a){
  vec ret(a.size());
  for (int i=0; i< a.size(); i++)
   ret[i] = (double)a[i];
  return ret;
}
inline mat fmat2mat(const fmat  & a){
  mat ret(a.rows(), a.cols());
  for (int i=0; i< a.rows(); i++)
    for (int j=0; j< a.cols(); j++)
      ret.set(i,j,a.get(i,j));
  return ret;
}
inline fmat mat2fmat(const mat  & a){
  fmat ret(a.rows(), a.cols());
  for (int i=0; i< a.rows(); i++)
    for (int j=0; j< a.cols(); j++)
      ret.set(i,j,a.get(i,j));
  return ret;
}


inline fvec frandu(int size){
  return vec2fvec(randu(size));
}

inline void debug_print_vec(const char * name,const fvec& _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec[i] == 0)
      printf("      0    ");
    else printf("%12.4f    ", _vec[i]);
  printf("\n");
}
inline float get(sparse_fvec & v1, int pos){
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
inline float min( sparse_fvec & dvec){
 
  float dmin = 1e100;
  FOR_ITERATOR(i, dvec){
     dmin = std::min(dmin, get_nz_data(dvec, i));
  }
  return dmin;
}

inline float max( sparse_fvec & dvec){
 
  float dmax = -1e100;
  FOR_ITERATOR(i, dvec){
     dmax = std::max(dmax, get_nz_data(dvec, i));
  }
  return dmax;
}
inline void plus_mul( fvec &v1,  sparse_fvec &v2, float factor){
  FOR_ITERATOR(i, v2){  
    v1[get_nz_index(v2, i)] += factor*get_nz_data(v2, i);
  }
}
inline float sum( sparse_fvec & dvec){
  float sum = 0;
  FOR_ITERATOR(i, dvec){
     sum += get_nz_data(dvec, i);
  }
  return sum;
}
inline void assign(fvec & v1, sparse_fvec & v2, int N){
  v1 = fzeros(N);
  FOR_ITERATOR(i, v2){
     v1[get_nz_index(v2, i)] = get_nz_data(v2, i);
  }

}
inline fvec init_fvec(float * array, int size){
  return fvec(array, size);
}

inline void print(sparse_fvec & vec){
  int cnt = 0;
  FOR_ITERATOR(i, vec){
    std::cout<<get_nz_index(vec, i)<<":"<< get_nz_data(vec, i) << " ";
    cnt++;
    if (cnt >= 20)
       break;
  }
  std::cout<<std::endl;
}




#endif

#endif //eigen
#endif //MATH_LAYER_FLOAT_GRAPHLAB
