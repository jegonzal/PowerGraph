#ifndef _ITPP_VEC_UTILS_H
#define _ITPP_VEC_UTILS_H

#include <itpp/itbase.h>
/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU


*/
using namespace itpp;

namespace graphlab {

template<>
oarchive& operator<< <itpp::Vec<double> > (oarchive& arc, const itpp::Vec<double> &vec) {
  arc << vec.length();
  serialize(arc, vec._data(), sizeof(double)*vec.length());
  return arc;
}
template<>
oarchive& operator<< <itpp::Mat<double> > (oarchive& arc, const itpp::Mat<double> &mat) {
  arc << mat.rows() << mat.cols();
  serialize(arc, mat._data(), sizeof(double)*mat._datasize());  
  return arc;
}

template<>
iarchive& operator>> <itpp::Vec<double> > (iarchive& arc, itpp::Vec<double> &vec) {
  size_t vlength;
  arc >> vlength;
  vec.set_size(vlength);
  deserialize(arc, vec._data(), sizeof(double)*vec.length());
  return arc;
}

template<>
iarchive& operator>> <itpp::Mat<double> > (iarchive& arc, itpp::Mat<double> &mat) {
  size_t rows, cols;
  arc >> rows >> cols;
  mat.set_size(rows,cols);
  deserialize(arc, mat._data(), sizeof(double)*mat._datasize());   
  return arc;
}


};

void debug_print_vec(const char * name,const vec& _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec[i] == 0)
      printf("      0    ");
    else printf("%12.4g    ", _vec[i]);
  printf("\n");
}

inline void dot2(const vec&  x1, const vec& x3, mat & Q, int j, int len){
	for (int i=0; i< len; i++){
		Q.set(i,j,(x1[i] * x3[i]));
	}
}


mat GenDiffMat(int K){
    mat ret(K,K); 
    ret.zeros();
    for (int i=0; i<K; i++){
        ret(i,i) = 2;
    }
    for (int i=1; i<K; i++){
	ret(i-1,i) = -1;
    }
    for (int i=1; i<K; i++){
	ret(i,i-1) = -1;
    }
   return ret;
}

#endif //_ITPP_VEC_UTILS_H
