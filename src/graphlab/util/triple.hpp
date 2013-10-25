/**  
 * Copyright (c) 2013 Shanghai Jiao Tong University. 
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
 * author: rong chen (rongchen@sjtu.edu.cn) 2013.7
 *
 */

#ifndef GRAPHLAB_TRIPLE_HPP
#define GRAPHLAB_TRIPLE_HPP


#include <iostream>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

namespace graphlab {

  template <typename _T1, typename _T2, typename _T3>
  struct triple {
    typedef _T1 first_type;
    typedef _T2 second_type;
    typedef _T3 third_type;
  
    first_type first;
    second_type second;
    third_type third;
    
    triple() : first(_T1()), second(_T2()), third(_T3()) {}
    triple(const _T1& x, const _T2& y, const _T3& z) : 
            first(x), second(y), third(z) {}
  
    triple(const triple<_T1, _T2, _T3>& o) : 
            first(o.first), second(o.second), third(o.third){}

    void load(iarchive& iarc) {
      iarc >> first;
      iarc >> second;
      iarc >> third;
    }

    void save(oarchive& oarc) const {
      oarc << first;
      oarc << second;
      oarc << third;
    }
  };

  template <typename _T1, typename _T2, typename _T3>
  inline bool operator == (const triple<_T1, _T2, _T3>& x,
                            const triple<_T1, _T2, _T3>& y) { 
    return x.first == y.first && x.second == y.second && x.third == y.third; 
  }
  
  template <typename _T1, typename _T2, typename _T3>
  inline bool operator < (const triple<_T1, _T2, _T3>& l, 
                          const triple<_T1, _T2, _T3>& r) { 
    return (l.first < r.first) || 
           (!(r.first < l.first) 
              && (l.second < r.second)) || 
           (!(r.first < l.first) 
              && !(r.second < l.second) 
              && (l.third < r.third)); 

  }

  template <typename _T1, typename _T2, typename _T3>
  inline bool operator != (const triple<_T1, _T2, _T3>& l, 
                            const triple<_T1, _T2, _T3>& r) {
    return !(l == r);
  }

  template <typename _T1, typename _T2, typename _T3>
  inline bool operator > (const triple<_T1, _T2, _T3>& l, 
                          const triple<_T1, _T2, _T3>& r) {
    return r < l;
  }

  template <typename _T1, typename _T2, typename _T3>
  inline bool operator <= (const triple<_T1, _T2, _T3>& l, 
                            const triple<_T1, _T2, _T3>& r) {
    return !(r < l);
  }

  template <typename _T1, typename _T2, typename _T3>
  inline bool operator >= (const triple<_T1, _T2, _T3>& l, 
                            const triple<_T1, _T2, _T3>& r) {
    return !(l < r);
  }

  template <typename _T1, typename _T2, typename _T3>
  inline triple<_T1, _T2, _T3> make_triple(
        const _T1& x, const _T2& y, const _T3& z) {
    return triple<_T1, _T2, _T3>(x, y, z);
  }
}; // end of graphlab namespace

#endif

