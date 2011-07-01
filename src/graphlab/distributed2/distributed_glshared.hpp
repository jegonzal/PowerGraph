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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_GLSHARED_HPP
#define GRAPHLAB_DISTRIBUTED_GLSHARED_HPP
#include <vector>
#include <boost/shared_ptr.hpp>
#include <graphlab/shared_data/glshared.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/distributed2/distributed_glshared_base.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
namespace graphlab {

class distributed_glshared_manager;



template <typename T>
class distributed_glshared: public distributed_glshared_base {

public:
  //! Type of the apply function inhereted from the gl_shared base
  typedef distributed_glshared_base::apply_function_type apply_function_type;
  //! Type of the boost shared pointer to a constant 
  typedef boost::shared_ptr<const T> const_ptr_type;
  //! Type of the boost shared pointer
  typedef boost::shared_ptr<T> ptr_type;


private:
  // two instances of the data
  T head;
  ptr_type ptr;
  // A lock used to sequentialize multiple writes
  rwlock lock;

  /**
    check if the backend storage has been modified
  */
  void check_invalidation(bool async = false) const;
  
  void issue_modification(bool async = false);
public:
  //! Construct initial shared pointers
  distributed_glshared() {
    any a(T()); // force instantiation of the serializer
    ptr.reset(&head, glshared_impl::empty_deleter<T>());
  }


  /// Returns a copy of the data
  inline T get_val() const{
    lock.readlock();
    T ret = head;
    lock.rdunlock();
    return ret;
  }

  /**
   * Gets the value of the shared variable wrapped in an any.
   */
  any get_any() const {
    return get_val();
  }

  /**
   * Returns true if there are no other active references to this
   * variable.
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
   * get_ptr not supported yet for distributed glshared
   */
  inline const_ptr_type get_ptr() __attribute__(( deprecated ));

  //// Unfortunate error is not supported on gcc 4.2 compilers 
  //   inline const_ptr_type get_ptr() 
  //     __attribute__(( error ("get_ptr not supported yet for distributed glshared") ));
  //  const{ return ptr; }

  /**
   * changes the data to 't'. This operation is atomic This
   * operation could stall forever if there are active shared
   * pointers to this variable which are never released.
   */
  void set(const T& t) {
    lock.writelock();
    head = t;
    lock.wrunlock();
    issue_modification(true);
  }


  /**
   * Sets the value of the shared variable using an any. The type of
   * the any must match the type of the shared object.
   */
  void set_any(const any &t) {
    set(t.as<T>());
  }

  /** 
   * Exchanges the data with 't'. This operation is atomic This
   * operation could stall forever if there are active shared
   * pointers to this variable which are never released.
   */
  void exchange(T& t);


  /**
   * apply's a function to this variable passing an additional
   * parameter.  This operation could stall forever if there are
   * active shared pointers to this variable which are never
   * released.
   */
  void apply(apply_function_type fun,
             const any& srcd);
  
  void save(oarchive &oarc) const{
    lock.readlock();
    oarc << head;
    lock.rdunlock();
  }
  void load(iarchive &iarc) {
    lock.writelock();
    iarc >> head;
    lock.wrunlock();
  }
  
  const char* type_name() const {
    return typeid(T).name();
  }

  procid_t preferred_machine() const;
  
};



} // namespace gaphlab


#include <graphlab/distributed2/distributed_glshared_manager.hpp>

namespace graphlab {
template <typename T>
void distributed_glshared<T>::check_invalidation(bool async) const{
  if (manager) {
    if (atomic_compare_and_swap(invalidated, true, false)) {
      manager->read_synchronize(id, async);
    }
  }
}
   
template<typename T>
void distributed_glshared<T>::issue_modification(bool async) {
  if (manager) {
    manager->write_synchronize(id, async);
  }
}


template<typename T>
void distributed_glshared<T>::exchange(T& t) {
    // If I have a manager, I need to do this through the manager to get
    // atomicity. If not I can just do it locally

  if (manager) {
    std::stringstream strm;
    oarchive oarc(strm);
    oarc << t;
    std::stringstream retstrm(manager->exchange(id, strm.str()));
    iarchive iarc(retstrm);
    iarc >> t;
  }
  else {
    lock.writelock();
    T retval = get_val();
    head = t;
    t = retval;
    lock.wrunlock();
  }
}

template<typename T>
void distributed_glshared<T>::apply(apply_function_type fun,
                                    const any& srcd) {
  if (manager) {
    // serialize the srcd
    manager->apply<T>(id, fun, srcd);
  }
  else {
    lock.writelock();
    any temp = head;
    fun(temp, srcd);
    head = temp.as<T>();
    lock.writelock();
  }
}

template <typename T>
procid_t distributed_glshared<T>::preferred_machine() const {
  if (manager) {
    return manager->preferred_machine(id);
  }
  else {
    return 0;
  }
}
}
#endif //GRAPHLAB_DISTRIBUTED_GLSHARED_HPP

