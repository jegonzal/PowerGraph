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


#ifndef GRAPHLAB_GLSHARED_HPP
#define GRAPHLAB_GLSHARED_HPP
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <boost/type_traits/remove_reference.hpp>

#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {

  /**
   * \brief A shared data entry.  
   *
   * glshared<datatype> variable; will create a shared variable of the
   * defined datatype. 
   */
  template <typename T>
  class glshared : public iglshared_base{

  public:
    //! The type of the internal element
    typedef T contained_type;
    //! Type of the apply function inhereted from the gl_shared base
    typedef glshared_base::apply_function_type apply_function_type;

    //! Reference holder
    class const_ref {
      const T& data;
      rwlock& lock;
      const_ref(); // Not default constructable
      const_ref(const const_ref&); // Not copyable 
      operator=(const const_ref&); // not assignable
    public:
      const_ref(const glshared& shared) : 
        data(shared.data), lock(shared.lock) { 
        lock.read_lock();
      }
      ~const_ref() { lock.unlock(); }
      const T& val() { return data; }
    };

  private:
    //! The interal storage item
    contained_type data;
    rwlock lock;

  public:
    //! Construct initial shared pointers
    glshared() { }
  
    //! Add a delta function:
    void operator+=(const T& other) { 
      lock.writelock();
      data += other;
      lock.unlock();
    }


    /// Returns a copy of the data
    inline T val() const {
      T copy;
      lock.readlock();
      copy = data;
      lock.unlock();
      return copy;
    }

    /**
     * Gets the value of the shared variable wrapped in an any.
     */
    any get_any() const {
      return any<T>(val());
    }
  
    /**
     * Returns true if there are no other active references to this
     * variable. This should not be used to test for exclusive access,
     * and is meant for internal use.
     */
    bool is_unique() const {
      return true;
    }
  
    /**
     * Returns a shared_ptr to the data.  When the shared pointer goes
     * out of scope, its contained pointer becomes invalid. The user
     * should not request for the underlying pointer inside the
     * shared_ptr.
     *
     *  gl_shared<int> shared_x;
     * 
     *  Ok: 
     *    boost::shared_ptr<const int> var_p = shared_x.get_ptr();
     *    const int& x = *var_p;
     *
     *  Bad:
     *    const int& x = *shared_x.get_ptr();
     *
     *
     */
    inline const_ptr_type get_ptr() const{
      return boost::const_pointer_cast<const T, T>(*head);
    }

    /**
     * changes the data to 't'. This operation is atomic This
     * operation could stall forever if there are active shared
     * pointers to this variable which are never released.
     */
    void set(const T& t) {
      set_lock.lock();
      wait_for_buffer_release();
      *(*buffer) = t;
      exchange_buffer_and_head();
      set_lock.unlock();
    }
  
  
    /**
     * Sets the value of the shared variable using an any. The type of
     * the any must match the type of the shared object.
     */
    void set_any(const any &t) {
      set(t.as<T>());
    }
  
    /** 
     * Exchanges the data with 't'. This operation performs the
     * exchange atomically and could therefore stall
     * forever if there are active shared
     * pointers to this variable which are never released.
     */
    void exchange(T& t) {
      set_lock.lock();
      wait_for_buffer_release();
      T retval = get_val();
      *(*buffer) = t;
      t = retval;
      exchange_buffer_and_head();
      set_lock.unlock();
    }


    /**
     * apply's a function to this variable passing an additional
     * parameter. This operation tries to perform the modification atomically
     * and could stall forever if there are
     * active shared pointers to this variable which are never
     * released.
     */
    void apply(apply_function_type fun,
               const any& srcd) {
      set_lock.lock();
      wait_for_buffer_release();
      any temp = *(*head);
      fun(temp, srcd);
      *(*buffer) = temp.as<T>();
      exchange_buffer_and_head();
      set_lock.unlock();
    }
  };


  /**
   * Because it is inconvenient that the apply function specification
   * takes the current value as an "any" as opposed to using the true
   * type of the current value (T), this function adapts an apply
   * function written in the more intuitive form:
   *
   *  void applyfn(T&, const any&)
   *
   * to the regular apply function type.  apply_adapter<T, applyfn> is
   * a function which matches the regular apply function type and
   * calls applyfn.
   */
  template<typename T, void (*applyfn)(T&, const any&) >  
  void apply_adapter(any &d, const any& param) {
    applyfn(d.as<T>(), param);
  }

} 
#endif

