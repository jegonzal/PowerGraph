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


#ifndef GRAPHLAB_ATOMIC_OPS_HPP
#define GRAPHLAB_ATOMIC_OPS_HPP

#include <stdint.h>


namespace graphlab {
  /**
   * \ingroup util
     atomic instruction that is equivalent to the following:
     \code
     if (a==oldval) {    
       a = newval;           
       return true;          
     }
     else {
       return false;
    }
    \endcode
  */
  template<typename T>
  bool atomic_compare_and_swap(T& a, T oldval, T newval) {
    return __sync_bool_compare_and_swap(&a, oldval, newval);
  };

  /**
   * \ingroup util
     atomic instruction that is equivalent to the following:
     \code
     if (a==oldval) {    
       a = newval;           
       return true;          
     }
     else {
       return false;
    }
    \endcode
  */
  template<typename T>
  bool atomic_compare_and_swap(volatile T& a, 
                               T oldval, 
                               T newval) {
    return __sync_bool_compare_and_swap(&a, oldval, newval);
  };

  /**
   * \ingroup util
     atomic instruction that is equivalent to the following:
     \code
     if (a==oldval) {    
       a = newval;           
       return true;          
     }
     else {
       return false;
    }
    \endcode
  */
  template <>
  inline bool atomic_compare_and_swap(volatile double& a, 
                                      double oldval, 
                                      double newval) {
    volatile uint64_t* a_ptr = reinterpret_cast<volatile uint64_t*>(&a);
    const uint64_t* oldval_ptr = reinterpret_cast<const uint64_t*>(&oldval);
    const uint64_t* newval_ptr = reinterpret_cast<const uint64_t*>(&newval);
    return __sync_bool_compare_and_swap(a_ptr, *oldval_ptr, *newval_ptr);
  };

  /**
   * \ingroup util
     atomic instruction that is equivalent to the following:
     \code
     if (a==oldval) {    
       a = newval;           
       return true;          
     }
     else {
       return false;
    }
    \endcode
  */
  template <>
  inline bool atomic_compare_and_swap(volatile float& a, 
                                      float oldval, 
                                      float newval) {
    volatile uint32_t* a_ptr = reinterpret_cast<volatile uint32_t*>(&a);
    const uint32_t* oldval_ptr = reinterpret_cast<const uint32_t*>(&oldval);
    const uint32_t* newval_ptr = reinterpret_cast<const uint32_t*>(&newval);
    return __sync_bool_compare_and_swap(a_ptr, *oldval_ptr, *newval_ptr);
  };

  /** 
    * \ingroup util
    * \brief Atomically exchanges the values of a and b.
    * \warning This is not a full atomic exchange. Read of a,
    * and the write of b into a is atomic. But the write into b is not.
    */
  template<typename T>
  void atomic_exchange(T& a, T& b) {
    b = __sync_lock_test_and_set(&a, b);
  };

  /** 
    * \ingroup util
    * \brief Atomically exchanges the values of a and b.
    * \warning This is not a full atomic exchange. Read of a,
    * and the write of b into a is atomic. But the write into b is not.
    */
  template<typename T>
  void atomic_exchange(volatile T& a, T& b) {
    b = __sync_lock_test_and_set(&a, b);
  };

  /** 
    * \ingroup util
    * \brief Atomically sets a to the newval, returning the old value
    */
  template<typename T>
  T fetch_and_store(T& a, const T& newval) {
    return __sync_lock_test_and_set(&a, newval);
  };

}
#endif

