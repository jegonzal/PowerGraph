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
   *
   */
  inline __attribute((error("get_ptr not supported yet for distributed glshared"))) const_ptr_type get_ptr() const{
    return ptr;
  }

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

