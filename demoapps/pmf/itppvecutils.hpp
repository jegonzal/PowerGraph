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
/*
template<>
graphlab::oarchive& operator<< <QQR> (graphlab::oarchive&, const QQR &data) {
  // TODO
}

template<>
graphlab::iarchive& operator>> <QQR> (graphlab::iarchive&, QQR &data) {
  // TODO
}
template<>
graphlab::oarchive& operator<< <mult_QQR> (graphlab::oarchive&, const mult_QQR &data) {
  // TODO
}

template<>
graphlab::iarchive& operator>> <mult_QQR> (graphlab::iarchive&, mult_QQR &data) {
  // TODO
}
template<>
graphlab::oarchive& operator<< <mult_vec> (graphlab::oarchive&, const mult_vec &data) {
  // TODO
}

template<>
graphlab::iarchive& operator>> <mult_vec> (graphlab::iarchive&, mult_vec &data) {
  // TODO
}
*/



};

void debug_print_vec(const char * name,const vec& _vec, int len){
  printf("%s ) ", name);
  for (int i=0; i< len; i++)
    if (_vec[i] == 0)
      printf("      0    ");
    else printf("%12.4g    ", _vec[i]);
  printf("\n");
}



 double * vec2vec(const vec * _vec){
	double * ret = new double[_vec->size()];
        for (int i=0; i< _vec->size(); i++)
            ret[i] = _vec->get(i);
        return ret;
}
itpp::vec * vec2vec(const double * _vec, int len){
	itpp::vec * ret = new vec(len);
        //for (int i=0; i< len; i++)
        //    ret->set(i, _vec[i]);
        memcpy(ret->_data(), _vec, len*sizeof(double));
        return ret;
}
void vec2vec(double * _vec, const itpp::vec & _vec2, int len){
        //for (int i=0; i< len; i++)
        //    _vec[i] = _vec2[i];
        memcpy(_vec, _vec2._data(), len*sizeof(double));
}
void vec2vec2(const double * _vec, itpp::vec & _vec2, int len){
      _vec2 = vec(len);  
      //for (int i=0; i< len; i++)
      //      _vec2[i] = _vec[i];
      memcpy(_vec2._data(), _vec, len*sizeof(double));
}


 vec dot(const double * x1, const double * x2, int len){
     vec ret(len);
     for (int i=0; i< len; i++)
        ret(i) = x1[i]*x2[i];
    return ret;
 }
 
inline void dot2(double * x1, const vec & x3, double * ret, int len){
	for (int i=0; i< len; i++){
		ret[i] = (x1[i] * x3[i]);
	}
}
inline void dot2(const vec& x1, const vec & x3, double * ret, int len){
	for (int i=0; i< len; i++){
		ret[i] = (x1[i] * x3[i]);
	}
}


inline void dot2(double * x1, const double * x3, mat & Q, int j, int len){
	for (int i=0; i< len; i++){
		Q.set(i,j,(x1[i] * x3[i]));
	}
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
