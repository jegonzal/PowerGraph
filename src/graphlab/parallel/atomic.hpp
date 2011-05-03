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

#ifndef GRAPHLAB_ATOMIC_HPP
#define GRAPHLAB_ATOMIC_HPP

#include <stdint.h>

namespace graphlab {
  /**
   * \brief atomic object toolkit
   * \ingroup util
   * A templated class for creating atomic numbers.
   */
  template<typename T>
  class atomic{
  public:
    /// The current value of the atomic number
    volatile T value;
    
    /// Creates an atomic number with value "value"
    atomic(const T& value = 0) : value(value) { }
    
    /// Performs an atomic increment by 1, returning the new value
    T inc() { return __sync_add_and_fetch(&value, 1);  }
    /// Performs an atomic decrement by 1, returning the new value
    T dec() { return __sync_sub_and_fetch(&value, 1);  }
    /// Performs an atomic increment by 'val', returning the new value
    T inc(T val) { return __sync_add_and_fetch(&value, val);  }
    /// Performs an atomic decrement by 'val', returning the new value
    T dec(T val) { return __sync_sub_and_fetch(&value, val);  }
    
    /// Performs an atomic increment by 1, returning the old value
    T inc_ret_last() { return __sync_fetch_and_add(&value, 1);  }
    /// Performs an atomic decrement by 1, returning the old value
    T dec_ret_last() { return __sync_fetch_and_sub(&value, 1);  }
    /// Performs an atomic increment by 'val', returning the old value
    T inc_ret_last(T val) { return __sync_fetch_and_add(&value, val);  }
    /// Performs an atomic decrement by 'val', returning the new value
    T dec_ret_last(T val) { return __sync_fetch_and_sub(&value, val);  }

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
  bool atomic_compare_and_swap(T& a, const T &oldval, const T &newval) {
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
                               const T &oldval, 
                               const T &newval) {
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
  inline bool atomic_compare_and_swap(double& a, 
                                      const double &oldval, 
                                      const double &newval) {
    return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(&a), 
                                        *reinterpret_cast<const uint64_t*>(&oldval), 
                                        *reinterpret_cast<const uint64_t*>(&newval));
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
  inline bool atomic_compare_and_swap(float& a, const float &oldval, const float &newval) {
    return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(&a), 
                                        *reinterpret_cast<const uint32_t*>(&oldval), 
                                        *reinterpret_cast<const uint32_t*>(&newval));
  };

  /** 
    * \ingroup util
    * \brief Atomically exchanges the values of a and b.
    */
  template<typename T>
  void atomic_exchange(T& a, T& b) {
    b =__sync_lock_test_and_set(&a, b);
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
