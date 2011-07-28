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


#ifndef GRAPHLAB_SHARED_DATA_OPS
#define GRAPHLAB_SHARED_DATA_OPS

#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/generics/any.hpp>


namespace graphlab {

  /**
   * Some basic apply operations:
   * Usage:
   *   apply_ops<double>::increment_apply;
   */
  struct glshared_apply_ops {
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


    
  /**
   * Some basic sync operations (behave like folds)
   */
  template<typename Graph>
  struct glshared_sync_ops {
    typedef typename Graph::vertex_data_type vertex_data_type;
    typedef iscope<Graph> iscope_type;


    /**
     * Compute the sum of the values
     *   result = sum x[i]
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void sum(iscope_type& scope,
                    any& accumulator) {
      accumulator.as<AccumType>() += Getter(scope.vertex_data());
    }


    /**
     * Compute the L1 sum of the values
     *   result = sum abs(x[i])
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void l1_sum(iscope_type& scope,
                       any& accumulator) {
      accumulator.as<AccumType>() += std::abs(Getter(scope.vertex_data()));
    }

    /**
     * Compute the L2 sum of the values:
     *   result = sum x[i]^2
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void l2_sum(iscope_type& scope,
                       any& accumulator) {
      AccumType res = Getter(scope.vertex_data());
      accumulator.as<AccumType>() += res*res;
    }

    /**
     * Compute the max of the values:
     *   result = max {x_1, ... x_n}
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void max(iscope_type& scope,
                    any& accumulator) {
      accumulator.as<AccumType>() =
        std::max(accumulator.as<AccumType>(), Getter(scope.vertex_data()));
    }
  };
  
  
  
  
  
    
  /**
   * Some basic sync operations (behave like folds)
   */
  struct glshared_merge_ops {
    /**
     * Compute the sum of the values
     *   result = sum x[i]
     */
    template<typename AccumType>  
    static void sum(any& dest,
                    const any& src) {
      dest.as<AccumType>() += src.as<AccumType>();
    }


    /**
     * Compute the L1 sum of the values
     *   result = sum abs(x[i])
     */
    template<typename AccumType>  
    static void l1_sum(any& dest,
                       const any& src) {
      dest.as<AccumType>() += std::abs(src.as<AccumType>());
    }

    /**
     * Compute the L2 sum of the values:
     *   result = sum x[i]^2
     */
    template<typename AccumType>  
    static void l2_sum(any& dest,
                       const any& src) {
      dest.as<AccumType>() += src.as<AccumType>() * src.as<AccumType>();
    }

    /**
     * Compute the max of the values:
     *   result = max {x_1, ... x_n}
     */
    template<typename AccumType>  
    static void max(any& dest,
                    const any& src) {
      dest.as<AccumType>() =
        std::max(src.as<AccumType>(), dest.as<AccumType>());
    }
  };
}; // end of namespace graphlab



#endif

