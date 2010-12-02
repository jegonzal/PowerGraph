#ifndef GLSHARED_HPP
#define GLSHARED_HPP
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {






class glshared_base{
 public:

  typedef void(*apply_function_type)(any& current_data, const any& param);

  virtual any get_any() const = 0;
  virtual void set_any(const any&) = 0;
  virtual void apply(apply_function_type fun,
                         const any& srcd) = 0;
  virtual bool is_unique() const = 0;
};



namespace glshared_impl {
  template <typename T>
  struct empty_deleter {
    void operator()(T* d) { }
  };
};


/**
 * A shared data entry
 * glshared<datatype> variable;
 * 
 */
template <typename T>
class glshared:public glshared_base{
 private:
  T buffer_and_head[2];
  boost::shared_ptr<T> buffer_and_head_ptr[2];
  boost::shared_ptr<T>* buffer;
  boost::shared_ptr<T>* head;
  
  
  mutex set_lock;

  
  inline void wait_for_buffer_release() {
    while(!(buffer->unique())) sched_yield();
  }
  
  inline void exchange_buffer_and_head() {
    atomic_exchange(buffer, head);
  }
  
 public:
    typedef glshared_base::apply_function_type apply_function_type;
  

  glshared() {
    buffer_and_head_ptr[0].reset(&(buffer_and_head[0]), glshared_impl::empty_deleter<T>());
    buffer_and_head_ptr[1].reset(&(buffer_and_head[1]), glshared_impl::empty_deleter<T>());
    buffer = &(buffer_and_head_ptr[0]);
    head = &(buffer_and_head_ptr[1]);
  }
  

  /// Returns a copy of the data
  inline T get_val() const{
    return *(*(head));
  }

  any get_any() const {
    return *(*(head));
  }
  
  bool is_unique() const {
    return buffer->unique() && head->unique();
  }
  
  /**
   * Returns a shared_ptr to the data.
   * When the shared pointer goes out of scope, its contained
   * pointer becomes invalid. The user should not request
   * for the underlying pointer inside the shared_ptr.
   */
  inline boost::shared_ptr<const T> get_ptr() const{
    return boost::const_pointer_cast<const T, T>(*head);
  }

  /** changes the data to 't'. This operation is atomic
   * This operation could stall forever if there are active shared pointers
   * to this variable which are never released.
   */
  void set(const T& t) {
    set_lock.lock();
    wait_for_buffer_release();
    *(*buffer) = t;
    exchange_buffer_and_head();
    set_lock.unlock();
  }

  void set_any(const any &t) {
    set(t.as<T>());
  }
  
  /** Exchanges the data with 't'. This operation is atomic
   * This operation could stall forever if there are active shared pointers
   * to this variable which are never released.
   */
  void exchange(T& t) {
    set_lock.lock();
    wait_for_buffer_release();
    T retval = get_val();
    *(*buffer) = t;
    exchange_buffer_and_head();
    set_lock.unlock();
  }
  /**
   * apply's a function to this variable passing an additional parameter.
   * This operation could stall forever if there are active shared pointers
   * to this variable which are never released.
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


}
#endif