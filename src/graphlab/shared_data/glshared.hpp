#ifndef GLSHARED_HPP
#define GLSHARED_HPP
#include <boost/shared_ptr.hpp>
#include <graphlab/engine/iengine.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {


typedef void(*sync_function_type)(iscope_type& scope,
                                  any& accumulator);
  
typedef void(*apply_function_type)(T& current_data,
                                    const T& new_data);

typedef void(*merge_function_type)(any& merge_dest,
                                    const any& merge_src);


struct gl_shared_deleter{
  rwlock &lock;
  gl_shared_deleter(rwlock &lock_):lock(lock_) { }
  
  operator() {
    lock.unlock();
  }
};

/**
 * A shared data entry
 * glshared<datatype> variable;
 * 
 */
template <typename T>
class glshared{
 private:
  T buffer_and_head[2];
  shared_ptr buffer_and_head_ptr[2];
  shared_ptr* buffer;
  shared_ptr* head;
  
  iengine* engine;
  
  mutex set_lock;

  wait_for_buffer_release() {
    while(!(buffer->unique())) sched_yield;
  }
  
  exchange_buffer_and_head() {
    atomic_exchange(buffer, head);
  }
  
 public:
  gl_shared():engine(NULL) { 
    buffer_and_head_ptr[0].reset(&(buffer_and_head[0]));
    buffer_and_head_ptr[1].reset(&(buffer_and_head[1]));
    buffer = &(buffer_and_head_ptr[0]);
    head = &(buffer_and_head_ptr[1]);
  }
  
  /**
   * Associates this shared data variable with a sync and apply
   * operation in the engine.
   * The sync operation will be performed every 'sync_interval' updates.
   * and will perform a fold over all vertices within the range
   * [rangelow, rangehigh] inclusive
   */
  void associate_with(iengine* eng, 
                      sync_function_type sync_fun,
                      apply_function_type apply_fun,
                      const any& zero,
                      size_t sync_interval = -1,
                      size_t rangelow = 0,
                      size_t rangehigh = -1) {
    engine = eng;
  }
  
  /// Returns a copy of the data
  inline T get_val() const{
    return *(*(head));
  }
  
  /**
   * Returns a shared_ptr to the data.
   * When the shared pointer goes out of scope, its contained
   * pointer becomes invalid. The user should not request
   * for the underlying pointer inside the shared_ptr.
   */
  inline shared_ptr<const T> get_ptr() const{
    return (*head);
  }

  /// changes the data to 't'. This operation is atomic
  void set(const T& t) {
    set_lock.lock();
    wait_for_buffer_release();
    *(*buffer) = t;
    exchange_buffer_and_head();
    set_lock.unlock();
  }
  
  /// Exchanges the data with 't'. This operation is atomic
  void exchange(T& t) {
    set_lock.lock();
    wait_for_buffer_release();
    T retval = get_val();
    *(*buffer) = t;
    exchange_buffer_and_head();
    set_lock.unlock();
  }
  
  /** Applies a function to the data together with an additional
   * argument 'srcd'. This operation is atomic.
  */
  void apply(apply_function_type fun,
             const any& srcd) {
    set_lock.lock();
    wait_for_buffer_release();

    any destd = get_val();
    fun(destd, srcd);
    *(*buffer) = destd.as<T>();
    
    *(*buffer) = t;
    exchange_buffer_and_head();
    set_lock.unlock();
  }
  /** 
   * Performs a sync operation on this shared variable immediately.
   * Result will be available when the function completes.
   * Requires an engine to be registered.
   */
  void sync_now(){
    ASSERT_NE(engine, NULL);
    engine->sync_now();
  }

  /** Requests a defered sync operation on this shared variable. 
   *  Requires an engine to be registered.
   */
  void sync_soon(){
    ASSERT_NE(engine, NULL);
    engine->trigger_sync();
  }
};

}
#endif