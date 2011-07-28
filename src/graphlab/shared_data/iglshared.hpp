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


#ifndef GRAPHLAB_IGLSHARED_HPP
#define GRAPHLAB_IGLSHARED_HPP

// #include <boost/shared_ptr.hpp>
// #include <boost/function.hpp>
// #include <boost/type_traits/function_traits.hpp>
// #include <boost/type_traits/remove_reference.hpp>

// #include <graphlab/parallel/atomic.hpp>
// #include <graphlab/parallel/pthread_tools.hpp>
// #include <graphlab/logger/assertions.hpp>


#include <graphlab/util/generics/any.hpp>


namespace graphlab {

  /**
   * Common base class for all glshared<T> objects.  Exposes a common
   *  interface allowing all glshared<T> objects to be manipulated in
   *  the same way.
   */
  class iglshared {
  public:
    /**
     * The type of an apply function. The apply function performs an
     *  atomic operation on the contents of a shared object. The apply
     *  function takes a reference to the current value of the object
     *  (current_value: wrapped inside an any), an additional
     *  parameter (param), and makes modifications to the current
     *  value.
     */
    typedef void(*apply_function_type)(any& current_data, const any& param);

    /**
     * Gets the value of the shared variable wrapped in an any.
     */
    virtual any get_any() const = 0;
  
    /**
     * Sets the value of the shared variable using an any. The type of
     * the any must match the type of the shared object.
     */
    virtual void set_any(const any&) = 0;
  
    // /**
    //  * Performs an atomic modification on the value of the shared
    //  * object.  essentially calls fun(current_value, srcd) where
    //  * current_value is the value of this variable wrapped inside an
    //  * any.
    //  */
    // virtual void apply(apply_function_type fun,
    //                    const any& srcd) = 0;
                         
    /**
     * Returns true if there are no other active references to this
     * variable.
     */
    virtual bool is_unique() const = 0;


    // /**
    //  * Return a pointer to an aggergator type.  Note that aggregators
    //  * must be freed by the caller.
    //  */
    // virtual iaggregator* new_aggregator() = 0;

    
    // /**
    //  * Because it is inconvenient that the apply function specification
    //  * takes the current value as an "any" as opposed to using the true
    //  * type of the current value (T), this function adapts an apply
    //  * function written in the more intuitive form:
    //  *
    //  *  void applyfn(T&, const any&)
    //  *
    //  * to the regular apply function type.  apply_adapter<T, applyfn> is
    //  * a function which matches the regular apply function type and
    //  * calls applyfn.
    //  */
    // template<typename T, void (*applyfn)(T&, const any&) >  
    // void static apply_adapter(any& d, const any& param) {
    //   applyfn(d.as<T>(), param);
    // } // end of apply adapter

  };



} 
#endif

