#ifndef OBJECT_REQUEST_FUTURE_HPP
#define OBJECT_REQUEST_FUTURE_HPP
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/function_ret_type.hpp>

namespace graphlab {


  /**
   * \ingroup rpc
   * The result of a future_remote_request call.
   * This class represents the outcome of a remote request sent to another
   * machine via the future_remote_request_call. The future_remote_request call
   * returns immediately with this object. Only when operator() is called on this
   * object, then it waits for a result from the remote machine.
   *
   * example:
   * \code
   * // this function returns immediately
   * graphlab::request_future<int> res = 
   *   rmi.future_remote_request(SOME_OTHER_MACHINE, 
   *                             function_which_returns_an_integer, ...);
   *
   * ... we can do other stuff ... 
   * // read the result, or wait for the result if it is not done yet.
   * int actual_result = res();
   * \endcode
   *
   * The future object holds a copy of the result of the request, and the
   * operator() call returns a reference to this result (once it is available).
   */
template <typename T>
struct request_future {
  typedef typename dc_impl::function_ret_type<T>::type result_type;
  mutable std::auto_ptr<dc_impl::reply_ret_type> reply;
  result_type result;
  bool hasval;

  /// default constructor
  request_future(): 
      reply(new dc_impl::reply_ret_type(REQUEST_WAIT_METHOD)),
      hasval(false) { }

  /** We can assign return values directly to the future in the
   * case where no remote calls are necessary. 
   * Thus allowing the following to be written easily:
   * \code
   * request_future<int> a_function(int arg) {
   *   if (arg == 0) return rmi.future_remote_request(... somewhere else ...) ;
   *   else return 10;
   * }
   * \endcode
   */
  request_future(const T& val): 
      reply(NULL),
      result(val), 
      hasval(true) { }

 
  /// copy constructor 
  request_future(const request_future<T>& val): 
      reply(val.reply),
      result(val.result), 
      hasval(val.hasval) { }

  /// operator=
  request_future& operator=(const request_future<T>& val) {
    reply = val.reply;
    result = val.result;
    hasval = val.hasval;
    return *this;
  }

  /// explicit call to wait(). Will wait only if the future has no value yet
  void wait() {
    if (!hasval) {
      reply->wait(); 
      iarchive iarc(reply->val.c, reply->val.len); 
      iarc >> result;  
      reply->val.free(); 
      hasval = true;
    }
  }

  /**
   * Waits for the request if it has not yet been received.
   * Otherwise, returns a reference to the received value.
   */
  result_type& operator()() {
    if (!hasval) wait();
    return result;
  }
};


template <>
struct request_future<void> {
  typedef typename dc_impl::function_ret_type<void>::type result_type;
  mutable std::auto_ptr<dc_impl::reply_ret_type> reply;
  bool hasval;

  request_future(): 
      reply(new dc_impl::reply_ret_type(REQUEST_WAIT_METHOD)),
      hasval(false) { }
  request_future(int val): 
      reply(NULL),
      hasval(true) { }
 
 
  request_future(const request_future<void>& val): 
      reply(val.reply),
      hasval(val.hasval) { }

  request_future& operator=(const request_future<void>& val) {
    reply = val.reply;
    hasval = val.hasval;
    return *this;
  }



  void wait() {
    if (!hasval) {
      result_type result;
      reply->wait(); 
      iarchive iarc(reply->val.c, reply->val.len); 
      iarc >> result;  
      reply->val.free(); 
      hasval = true;
    }
  }

  result_type operator()() {
    if (!hasval) wait();
    return 0;
  }
};


}
#endif
