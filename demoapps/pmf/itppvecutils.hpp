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
 */


#ifndef _ITPP_VEC_UTILS_H
#define _ITPP_VEC_UTILS_H

#include <itpp/itbase.h>
/**

 Probabalistic matrix/tensor factorization written Danny Bickson, CMU


*/
using namespace itpp;

namespace graphlab {
/*
template<>
inline oarchive& operator<< <itpp::Vec<double> > (oarchive& arc, const itpp::Vec<double> &vec) {
  arc << vec.length();
  serialize(arc, vec._data(), sizeof(double)*vec.length());
  return arc;
}
template<>
inline oarchive& operator<< <itpp::Mat<double> > (oarchive& arc, const itpp::Mat<double> &mat) {
  arc << mat.rows() << mat.cols();
  serialize(arc, mat._data(), sizeof(double)*mat._datasize());  
  return arc;
}

template<>
inline iarchive& operator>> <itpp::Vec<double> > (iarchive& arc, itpp::Vec<double> &vec) {
  size_t vlength;
  arc >> vlength;
  vec.set_size(vlength);
  deserialize(arc, vec._data(), sizeof(double)*vec.length());
  return arc;
}

template<>
inline iarchive& operator>> <itpp::Mat<double> > (iarchive& arc, itpp::Mat<double> &mat) {
  size_t rows, cols;
  arc >> rows >> cols;
  mat.set_size(rows,cols);
  deserialize(arc, mat._data(), sizeof(double)*mat._datasize());   
  return arc;
}
*/

};

inline void debug_print_vec(const char * name,const vec& _vec, int len){
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


#endif //_ITPP_VEC_UTILS_H
