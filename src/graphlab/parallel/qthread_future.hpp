#ifndef QTHREAD_EXTERNAL_FUTURE_HPP
#define QTHREAD_EXTERNAL_FUTURE_HPP
#include <boost/function.hpp>
#include <qthread.h>
namespace graphlab {

/**
 * Used by qthread functions to invoke an external operation
 * and to wait for a response.
 * The standard pattern is as follows
 *
 * \code
 * qthread_fn() {
 *   // create a future
 *   qthread_future<int> future;
 *   int& result = future.get();
 *   __call_external_function_to_do_stuff__(&result);
 *   // will wait until signal is raised on the future
 *   future.wait();
 * }
 *
 * void do_stuff(int* result) {
 *    .. doo stufff
 *    .. perhaps even give the handle and result pointers away and
 *    .. return.
 *    .. but eventually, someone somewhere must assign the result value
 *    .. and call qthread_future<int>::signal(result)
 * }
 *
 *
 * \endcode
 */
template <typename T>
class qthread_future {
 private:
  T _response;
 public:
  qthread_future() {
    // empty the var
    qthread_empty(reinterpret_cast<aligned_t*>(this));
  }

  /** Gets the address of the result of the future.
   * This should be passed together with the handle and should be
   * set before signal_future() is called.
   */
  T& get() {
    return _response;
  }

  /**
   * signals that a thread waiting on this future should wake up
   */
  static void signal(qthread_future<T>* handle) {
    qthread_fill(reinterpret_cast<aligned_t*>(handle));
  }

  /**
   * signals that a thread waiting on this future should wake up
   */
  static void signal(size_t handle) {
    qthread_fill(reinterpret_cast<aligned_t*>(handle));
  }


  bool is_ready() {
    return qthread_feb_status(reinterpret_cast<aligned_t*>(this));
  }


  /**
   * Waits for this future to be filled.
   */
  void wait() {
    qthread_readFF(NULL,reinterpret_cast<aligned_t*>(this));
  }
} __attribute__ ((aligned (8)));

} // namespace graphlab
#endif
