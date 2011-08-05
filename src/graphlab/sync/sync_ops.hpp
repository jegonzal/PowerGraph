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


#ifndef GRAPHLAB_SYNC_OPS
#define GRAPHLAB_SYNC_OPS

#include <graphlab/scope/iscope.hpp>


namespace graphlab {

  template<typename Graph>
  struct sync_ops {
    typedef iscope<Graph> iscope_type;

    template<typename Accum, Accum(*MapFun)(iscope_type& scope)>
    struct sum_group {
      Accum acc;
      sum_group() : acc() {}
      sum_group(const Accum& acc) : acc(acc) { }
      sum_group(iscope_type& scope) : acc(MapFun(scope)) { }
      operator Accum () const { return acc; }
      void operator+=(const sum_group& other) { acc += other.acc; }
    };      

    template<typename Accum, Accum(*MapFun)(iscope_type& scope)>
    struct max_group {
      Accum acc;
      max_group() : acc() {}
      max_group(const Accum& acc) : acc(acc) { }
      max_group(iscope_type& scope) : acc(MapFun(scope)) { }
      operator Accum() { return acc; }
      void operator+=(const max_group& other) { 
        acc = std::max(acc, other.acc); }
    };

    

    
    /**
     * Directly applies the new value in place of the old value:
     *  current value = new data;
     */
    template<typename T, typename Accum>
    static void identity(T& lvalue, const Accum& rvalue) {
      lvalue = rvalue;
    }

    /**
     * Same as the identity operation except that it prints the new
     * value before applying:
     *  print(new data);
     *  current value = new data;
     */
    template<typename T, typename Accum>
    static void identity_print(T& lvalue, const Accum& rvalue) {
      std::cout << rvalue << std::endl;
      identity<T, Accum>(lvalue, rvalue);
    }


    /**
     * Computes the square root:
     *  current value = sqrt(new_data);
     */
    template<typename T, typename Accum>
    static void sqrt(T& lvalue, const Accum& rvalue) {
      lvalue = sqrt(rvalue);
    }

    /**
     * Increments the current value by the new value
     *  current value += new value;
     */
    template<typename T, typename Accum>
    static void increment(T& lvalue, const Accum& rvalue) {
      lvalue += rvalue;
    }


    /**
     * Decrements the current value by the new value
     *  current value -= new value;
     */
    template<typename T, typename Accum>
    static void decrement(T& lvalue, const Accum& rvalue) {
      lvalue -= rvalue;
    }
  };

}; // end of namespace graphlab



#endif

