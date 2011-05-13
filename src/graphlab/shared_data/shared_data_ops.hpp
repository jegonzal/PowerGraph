/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_SHARED_DATA_OPS
#define GRAPHLAB_SHARED_DATA_OPS

#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/scope/iscope_factory.hpp>

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
    template<typename AccumType>
    static void identity(any& current_data,
                         const any& new_data) {
      current_data.as<AccumType>() = new_data.as<AccumType>();
    }

    /**
     * Same as the identity operation except that it prints the new
     * value before applying:
     *  print(new data);
     *  current value = new data;
     */
    template<typename AccumType>
    static void identity_print(any& current_data,
                               const any& new_data) {
      std::cout << new_data.as<AccumType>() << std::endl;
      current_data.as<AccumType>() = new_data.as<AccumType>();
    }


    /**
     * Computes the square root:
     *  current value = sqrt(new_data);
     */
    template<typename AccumType>
    static void sqrt(any& current_data,
                     const any& new_data) {
      current_data.as<AccumType>() =
        std::sqrt(new_data.as<AccumType>());
    }

    /**
     * Increments the current value by the new value
     *  current value += new value;
     */
    template<typename AccumType>
    static void increment(any& current_data,
                          const any& new_data) {
      current_data.as<AccumType>() += new_data.as<AccumType>();
    }


    /**
     * Decrements the current value by the new value
     *  current value -= new value;
     */
    template<typename AccumType>
    static void decrement(any& current_data,
                          const any& new_data) {
      current_data.as<AccumType>() -= new_data.as<AccumType>();
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

