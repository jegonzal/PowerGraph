#ifndef GRAPHLAB_SHARED_DATA_OPS
#define GRAPHLAB_SHARED_DATA_OPS

#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/scope/iscope_factory.hpp>
#include <graphlab/shared_data/ishared_data.hpp>

namespace graphlab {

  /**
   * \deprecated Use glshared and glshared_apply_ops
   * Some basic apply operations:
   * Usage:
   *   apply_ops<double>::increment_apply;
   */
  template<typename Graph>
  struct apply_ops {
    typedef ishared_data<Graph> ishared_data_type;
    /**
     * Directly applies the new value in place of the old value:
     *  current value = new data;
     */
    template<typename AccumType>
    static void identity(size_t index,
                         const ishared_data_type& shared_data,
                         any& current_data,
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
    static void identity_print(size_t index,
                               const ishared_data_type& shared_data,
                               any& current_data,
                               const any& new_data) {
      std::cout << new_data.as<AccumType>() << std::endl;
      current_data.as<AccumType>() = new_data.as<AccumType>();
    }


    /**
     * Computes the square root:
     *  current value = sqrt(new_data);
     */
    template<typename AccumType>
    static void sqrt(size_t index,
                     const ishared_data_type& shared_data,
                     any& current_data,
                     const any& new_data) {
      current_data.as<AccumType>() =
        std::sqrt(new_data.as<AccumType>());
    }

    /**
     * Increments the current value by the new value
     *  current value += new value;
     */
    template<typename AccumType>
    static void increment(size_t index,
                          const ishared_data_type& shared_data,
                          any& current_data,
                          const any& new_data) {
      current_data.as<AccumType>() += new_data.as<AccumType>();
    }


    /**
     * Decrements the current value by the new value
     *  current value -= new value;
     */
    template<typename AccumType>
    static void decrement(size_t index,
                          const ishared_data_type& shared_data,
                          any& current_data,
                          const any& new_data) {
      current_data.as<AccumType>() -= new_data.as<AccumType>();
    }      
  };

    
  /**
   * \deprecated Use glshared and glshared_sync_ops
   * Some basic sync operations (behave like folds)
   */
  template<typename Graph>
  struct sync_ops {
    typedef typename Graph::vertex_data_type vertex_data_type;
    typedef ishared_data<Graph> ishared_data_type;
    typedef iscope<Graph> iscope_type;


    /**
     * Compute the sum of the values
     *   result = sum x[i]
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void sum(size_t index,
                    const ishared_data_type& shared_data,
                    iscope_type& scope,
                    any& accumulator) {
      accumulator.as<AccumType>() += Getter(scope.vertex_data());
    }


    /**
     * Compute the L1 sum of the values
     *   result = sum abs(x[i])
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void l1_sum(size_t index,
                       const ishared_data_type& shared_data,
                       iscope_type& scope,
                       any& accumulator) {
      accumulator.as<AccumType>() += std::abs(Getter(scope.vertex_data()));
    }

    /**
     * Compute the L2 sum of the values:
     *   result = sum x[i]^2
     */
    template<typename AccumType, 
             AccumType (*Getter)(const vertex_data_type&) >  
    static void l2_sum(size_t index,
                       const ishared_data_type& shared_data,
                       iscope_type& scope,
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
    static void max(size_t index,
                    const ishared_data_type& shared_data,
                    iscope_type& scope,
                    any& accumulator) {
      accumulator.as<AccumType>() =
        std::max(accumulator.as<AccumType>(), Getter(scope.vertex_data()));
    }
  };






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
    typedef ishared_data<Graph> ishared_data_type;
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
}; // end of namespace graphlab



#endif
