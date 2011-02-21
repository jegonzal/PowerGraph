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
  T buffer_and_head[2];

  // shared pointers to the data
  boost::shared_ptr<T> buffer_and_head_ptr[2];

  // a pointer to the shared pointer which contains the write target
  boost::shared_ptr<T>* buffer;
  // a pointer to the shared pointer which contains the read target
  boost::shared_ptr<T>* head;

  // A lock used to sequentialize multiple writes
  mutex set_lock;

  // Waits until all references to the buffer are released
  inline void wait_for_buffer_release() {
    while(!(buffer->unique())) sched_yield();
  }

  // Performs an atomic exchange of the head and buffer pointers
  inline void exchange_buffer_and_head() {
    atomic_exchange(buffer, head);
  }
  /**
    check if the backend storage has been modified
  */
  void check_invalidation(bool async = false) const;
  
  void issue_modification(bool async = false);
public:
  //! Construct initial shared pointers
  distributed_glshared() {
    buffer_and_head_ptr[0].reset(&(buffer_and_head[0]), 
                                 glshared_impl::empty_deleter<T>());
    buffer_and_head_ptr[1].reset(&(buffer_and_head[1]), 
                                 glshared_impl::empty_deleter<T>());
    buffer = &(buffer_and_head_ptr[0]);
    head = &(buffer_and_head_ptr[1]);
  }


  /// Returns a copy of the data
  inline T get_val() const{
    check_invalidation(false);
    return *(*(head));
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
    return buffer->unique() && head->unique();
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
    check_invalidation(false);
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
    set_lock.lock();
    oarc << *(*head);
    set_lock.unlock();
  }
  void load(iarchive &iarc) {
    set_lock.lock();
    iarc >> *(*head);
    set_lock.unlock();
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
    set_lock.lock();
    wait_for_buffer_release();
    T retval = get_val();
    *(*buffer) = t;
    t = retval;
    exchange_buffer_and_head();
    set_lock.unlock();
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
    set_lock.lock();
    wait_for_buffer_release();
    any temp = *(*head);
    fun(temp, srcd);
    *(*buffer) = temp.as<T>();
    exchange_buffer_and_head();
    set_lock.unlock();
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
