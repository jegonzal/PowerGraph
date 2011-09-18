#ifndef MATH_LAYER_GRAPHLAB
#define MATH_LAYER_GRAPHLAB

//#define HAS_EIGEN
#ifdef HAS_EIGEN


#include "../../../deps/eigen-eigen-3.0.2/Eigen/Dense"
#include "../../../deps/eigen-eigen-3.0.2/Eigen/Cholesky"
#include "../../../deps/eigen-eigen-3.0.2/Eigen/Eigenvalues"
using namespace Eigen;

typedef MatrixXd mat;
typedef VectorXd vec;
typedef VectorXi ivec;

mat randn1(int dx, int dy, int col);

inline mat eye(int size){
  return mat::Identity(size, size);
}
inline vec ones(int size){
  return vec::Ones(size);
}
inline vec zeros(int size){
  return vec::Zero(size);
}
inline mat zeros(int rows, int cols){
  return mat::Zero(rows, cols);
}
inline void debug_print_vec(const char * name,const vec& _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec(i) == 0)
      printf("      0    ");
    else printf("%12.4g    ", _vec(i));
  printf("\n");
}

inline void dot2(const vec&  x1, const vec& x3, mat & Q, int j, int len){
	for (int i=0; i< len; i++){
		Q(i,j) = (x1(i) * x3(i));
	}
}

inline bool ls_solve_chol(const mat &A, const vec &b, vec &result){
    //result = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    result = A.ldlt().solve(b);
    return true;
}
inline bool ls_solve(const mat &A, const vec &b, vec &result){
    //result = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    result = A.ldlt().solve(b);
    return true;
}
inline bool chol(mat& sigma, mat& out){
  //out = sigma.llt();
  //TODO!!
   return true;
}
inline mat transpose(mat & A){
   return A.transpose();
}
inline void set_val(mat &A, int row, int col, double val){
  A(row, col) = val;
}
inline double get_val(const mat &A, int row, int col){
  return A(row, col);
}
inline vec get_col(const mat& A, int col){
  return A.col(col);
}
inline vec get_row(const mat& A, int row){
  return A.row(row);
}
inline void set_col(mat& A, int col, const vec & val){
  A.col(col) = val;
}
inline void set_row(mat& A, int row, const vec & val){
  A.row(row) = val;
}
inline mat randn(int dx, int dy){
  return randn1(dx,dy,-1);
}
inline void set_diag(mat &A, vec & v){
   A.diagonal()=v;
}
inline double sumsum(const mat & A){
   return A.sum();
}
inline mat init_mat(const char * string, int row, int col){
  mat out(row, col);
  char buf[2056];
  strcpy(buf, string);
  char *pch = strtok(buf," \r\n\t;");
  while (pch != NULL)
  {
     out << atof(pch);
     pch = strtok (NULL, " \r\n\t;");
  }
  return out;
}
inline vec init_vec(const char * string, int size){
  vec out(size);
  char buf[2056];
  strcpy(buf, string);
  char *pch = strtok (buf," \r\n\t;");
  while (pch != NULL)
  {
     out << atof(pch);
     pch = strtok (NULL, " \r\n\t;");
  }
  return out;
}
inline double norm(const mat &A, int pow=2){
     return A.squaredNorm();
}
inline mat inv(const mat&A){
   return A.inverse();
}
inline bool inv(const mat&A, mat &out){
   out = A.inverse();
   return true;
}
inline vec outer_product(const vec &a, const vec &b){
  return a*b.transpose();
}
inline bool eig_sym(const mat & T, vec & eigenvalues, mat & eigenvectors){
   VectorXcd eigs = T.eigenvalues();
   eigenvalues = eigs.real(); //TODO - what happen with complex
   return true;
}
inline vec head(const vec& v, int num){
   return v.head(num);
}
inline ivec head(const ivec& v, int num){
   return v.head(num);
}
inline vec elem_mult(const vec&a, const vec&b){
   vec ret = a;
   for (int i=0; i<b.size(); i++)
      ret(i) *= b(i);
   return ret;
}
inline double sum(const vec & a){
  return a.sum();
}
inline double sum_sqr(const vec & a){
  vec ret = a.array().pow(2);
  return ret.sum();
}
inline double trace(const mat & a){
  return a.trace();
}
inline double min(const vec &a){
  return a.minCoeff();
}
inline double max(const vec & a){
  return a.maxCoeff();
}
inline vec randu(int size){
  return vec::Random(size);
}
inline ivec randi(int size, int from, int to){
  ivec ret(size);
  for (int i=0; i<size; i++)
    ret[i]= internal::random<int>(from,to);
  return ret;
}
inline ivec concat(const ivec&a, const ivec&b){ 
   ivec ret(a.size()+b.size());
   ret << a,b;
   return ret;
}
inline void sort(ivec &a){
   std::sort(a.data(), a.data()+a.size());
}
inline ivec sort_index(const vec&a){
  ivec ret(a.size()); 
  std::vector<std::pair<double,int> > D;
  // 	
  D.reserve(a.size());
  for (int i=0;i<a.size();i++)
    D.push_back(std::make_pair<double,int>(a.coeff(i),i));
  std::sort(D.begin(),D.end());
  for (int i=0;i<a.size();i++)
  { 
    ret[i]=D[i].second;
  } 
  return ret;
}
inline void del(ivec&a, int i){
   ivec ret(a.size() - 1);
   ret << a.head(i-1), a.tail(a.size() -i-1);
   a= ret;
}
inline mat get_cols(const mat&A, ivec & cols){
  mat a(A.rows(), cols.size());
  for (int i=0; i< cols.size(); i++)
    set_col(a, i, get_col(A, cols[i]));
  return a;
}
inline void set_val(vec & v, int pos, double val){
  v(pos) = val;
}
inline double dot(const vec&a, const vec& b){
   return a.dot(b);
}
inline vec reverse(vec& a){
   return a.reverse();
}
inline ivec reverse(ivec& a){
   return a.reverse();
}

#else //eigen is not found
#if defined(HAS_ITPP)

#include <itpp/itbase.h>
#include <itpp/itstat.h>
#include <itpp/stat/misc_stat.h>
#include "itppvecutils.hpp"
using namespace itpp;

inline void set_val(mat& A, int row, int col, double val){
  A.set(row, col, val);
}
inline double get_val(const mat& A, int row, int col){
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
inline mat init_mat(const char * string, int row, int col){
  return mat(string);
}
inline vec init_vec(const char * string, int size){
  return vec(string);
}
inline vec head(const vec &v, int num){
  return v.mid(0,num);
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
inline mat get_cols(mat &A, ivec & ind){
  return A.get_cols(ind);
}
#endif

#endif //eigen
#endif //MATH_LAYER_GRAPHLAB
