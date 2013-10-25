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

#ifndef GRAPHLAB_TETRAD_HPP
#define GRAPHLAB_TETRAD_HPP


#include <iostream>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

namespace graphlab {

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  struct tetrad {
    typedef _T1 first_type;
    typedef _T2 second_type;
    typedef _T3 third_type;
    typedef _T4 fourth_type;
  
    first_type first;
    second_type second;
    third_type third;
    fourth_type fourth;
    
    tetrad() : first(_T1()), second(_T2()), third(_T3()), fourth(_T4()) {}
    tetrad(const _T1& x, const _T2& y, const _T3& z, const _T4& w) : 
            first(x), second(y), third(z), fourth(w) {}
  
    tetrad(const tetrad<_T1, _T2, _T3, _T4>& o) : 
            first(o.first), second(o.second), third(o.third), fourth(o.fourth){}

    void load(iarchive& iarc) {
      iarc >> first;
      iarc >> second;
      iarc >> third;
      iarc >> fourth;
    }

    void save(oarchive& oarc) const {
      oarc << first;
      oarc << second;
      oarc << third;
      oarc << fourth;
    }
  };

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline bool operator == (const tetrad<_T1, _T2, _T3, _T4>& x, 
                           const tetrad<_T1, _T2, _T3, _T4>& y) { 
    return x.first == y.first && x.second == y.second 
        && x.third == y.third && x.fourth == y.fourth; 
  }
  
  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline bool operator < (const tetrad<_T1, _T2, _T3, _T4>& l, 
                           const tetrad<_T1, _T2, _T3, _T4>& r) { 
    return (l.first < r.first) || 
           (!(r.first < l.first) 
              && (l.second < r.second)) || 
           (!(r.first < l.first) 
              && !(r.second < l.second) 
              && (l.third < r.third)) ||
           (!(r.first < l.first) 
              && !(r.second < l.second) 
              && !(r.third < l.third)
              && (l.fourth< r.fourth)); 

  }

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline bool operator != (const tetrad<_T1, _T2, _T3, _T4>& l, 
                            const tetrad<_T1, _T2, _T3, _T4>& r) {
    return !(l == r);
  }

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline bool operator > (const tetrad<_T1, _T2, _T3, _T4>& l, 
                          const tetrad<_T1, _T2, _T3, _T4>& r) {
    return r < l;
  }

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline bool operator <= (const tetrad<_T1, _T2, _T3, _T4>& l, 
                            const tetrad<_T1, _T2, _T3, _T4>& r) {
    return !(r < l);
  }

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline bool operator >= (const tetrad<_T1, _T2, _T3, _T4>& l, 
                            const tetrad<_T1, _T2, _T3, _T4>& r) {
    return !(l < r);
  }

  template <typename _T1, typename _T2, typename _T3, typename _T4>
  inline tetrad<_T1, _T2, _T3, _T4> make_tetrad(
        const _T1& x, const _T2& y, const _T3& z, const _T4& w) {
    return tetrad<_T1, _T2, _T3, _T4>(x, y, z, w);
  }
}; // end of graphlab namespace

#endif


