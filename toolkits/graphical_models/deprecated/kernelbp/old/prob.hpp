/*  
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


#ifndef NPROB_HPP
#define NPROB_HPP


#include <cmath>
//#include <itpp/itbase.h>
//#include <itpp/stat/misc_stat.h>
#define pi 3.14152965

#include "graphlab/util/random.hpp"

using namespace itpp;
using namespace std;


void randv(int n, vec & ret){
   assert(n>=1);
   for (int i=0; i< n; i++)
       //ret[i] = drand48();
       ret[i] = graphlab::random::rand01();
}
mat randn1(int Dx, int Dy){
  if (Dx == 0)
    Dx = 1;
  assert(Dy>=1);
  mat ret = zeros(Dx,Dy);
  vec us = zeros(ceil(Dx*Dy/2.0)*2); 
  randv(ceil(Dx*Dy/2.0)*2, us);
  int k=0;
  for (int i=0; i<Dx; i++){
     for (int j=0; j< Dy; j++){
         if (k % 2 == 0)
         	ret(i,j) = sqrt(-2.0*std::log(us[k/2]))*std::cos(2*pi*us[k/2+1]);
         else
         	ret(i,j) = sqrt(-2.0*std::log(us[k/2]))*std::sin(2*pi*us[k/2+1]);
         k++;
     }
  }
  assert(k == Dx*Dy);
  assert(ret.rows() == Dx && ret.cols() == Dy);
  return ret;
}



#endif //NPROB_HPP
