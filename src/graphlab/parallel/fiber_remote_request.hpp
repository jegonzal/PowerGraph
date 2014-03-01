/*  
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


#ifndef GRAPHLAB_FIBER_RPC_FUTURE_HPP
#define GRAPHLAB_FIBER_RPC_FUTURE_HPP
#include <graphlab/rpc/request_future.hpp>
#include <graphlab/rpc/request_reply_handler.hpp>
#include <graphlab/parallel/fiber_control.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
namespace graphlab {

/**
 * A implementation of the ireply_container interface
 * that will wait for rpc requests, but if the request is issued from within 
 * a fiber, will deschedule the fiber.
 */
struct fiber_reply_container: public dc_impl::ireply_container {
  dc_impl::blob val;
  mutex lock;
  conditional cond;
  // if wait is in a fiber, this will contain the ID of the fiber to wake up
  // If 0, the wait is not in a fiber.
  size_t waiting_tid;
  // true when the blob is assigned
  bool valready;

  fiber_reply_container():waiting_tid(0),valready(false) { }

  ~fiber_reply_container() {
    val.free();
  }

  void wait() {
    if (fiber_control::in_fiber()) {    
      // if I am in a fiber, use the deschedule mechanism 
      lock.lock();
      waiting_tid = fiber_control::get_tid();
      while(!valready) {
        // set the waiting tid value
        // deschedule myself. This will deschedule the fiber
        // and unlock the lock atomically
        fiber_control::deschedule_self(&lock.m_mut);
        // unlock the condition variable, this does not re-lock the lock
        lock.lock();
      }
      lock.unlock();
    } else {
      // Otherwise use the condition variable
      waiting_tid = 0;
      lock.lock();
      while(!valready) cond.wait(lock);
      lock.unlock();
    }
  }

  void receive(procid_t source, dc_impl::blob b) {
    lock.lock();
    val = b;
    valready = true;
    if (waiting_tid) {
      // it is a fiber. wake it up.
      fiber_control::schedule_tid(waiting_tid);
    } else {
      // not in fiber. This is just a condition signal
      cond.signal();
    }
    lock.unlock();
  }
  bool ready() const {
    return valready;
  }

  dc_impl::blob& get_blob() {
    return val;
  }
};


#if DOXYGEN_DOCUMENTATION


/**
 * \brief Performs a nonblocking RPC call to the target machine
 * to run the provided function pointer which has an expected return value.
 *
 * fiber_remote_request() calls the function "fn" on a target remote machine.
 * Provided arguments are serialized and sent to the target.
 * Therefore, all arguments are necessarily transmitted by value.
 * If the target function has a return value, it is sent back to calling
 * machine.  fiber_remote_request() returns immediately a \ref
 * graphlab::request_future object which will allow you wait for the return
 * value.
 *
 * fiber_remote_request() has an identical interface to 
 * \ref graphlab::distributed_control::future_remote_request() , but has the 
 * additional capability that if a \ref graphlab::request_future::wait() is 
 * called on the request while within a fiber, it deschedules the fiber and
 * context switches, returning only when the future is ready. This allows
 * the future to be used from within a fiber.
 *
 * Since this function is not a member of the distributed_control class,
 * it uses the function \ref distributed_control::get_instance() to obtain
 * the last instance of the distribute_control class created. This should be
 * sufficient for most use cases.
 *
 * \ref graphlab::object_fiber_remote_request is the version of this function
 * for remotely calling class member functions.
 *
 * Example:
 * \code
 * // A print function is defined
 * int add_one(int i) {
 *   return i + 1;
 * }
 *
 * ... ...
 * // call the add_one function on machine 1
 * int i = 10;
 * graphlab::request_future<int> ret = fiber_remote_request(1, add_one, i);
 * // this is safe to do within a fiber as it will not halt other fibers.
 * int result = ret();
 * // result will be 11
 * \endcode
 *
 * \see graphlab::distributed_control::remote_request
 *      graphlab::distributed_control::future_remote_request
 *      graphlab::object_fiber_remote_request
 *
 * \param targetmachine The ID of the machine to run the function on
 * \param fn The function to run on the target machine. Must be a pointer to
 *            member function in the owning object.
 * \param ... The arguments to send to Fn. Arguments must be serializable.
 *            and must be castable to the target types.
 *
 * \returns Returns a future templated around the same type as the return 
 *          value of the called function
 */
  request_future<RetVal> fiber_remote_request(procid_t targetmachine, Fn fn, ...);



/**
 * \brief Performs a nonblocking RPC call to the target machine
 * to run the provided function pointer which has an expected return value.
 *
 * object_fiber_remote_request() calls the function "fn" on a target remote machine.
 * Provided arguments are serialized and sent to the target.
 * Therefore, all arguments are necessarily transmitted by value.
 * If the target function has a return value, it is sent back to calling
 * machine.  object_fiber_remote_request() returns immediately a \ref
 * graphlab::request_future object which will allow you wait for the return
 * value.
 *
 * object_fiber_remote_request() has an identical interface to 
 * \ref graphlab::dc_dist_object::future_remote_request() , but has the 
 * additional capability that if a \ref graphlab::request_future::wait() is 
 * called on the request while within a fiber, it deschedules the fiber and
 * context switches, returning only when the future is ready. This allows
 * the future to be used from within a fiber.
 *
 * Since this function is not a member of the \ref dc_dist_object class,
 * it needs to be provided a reference to the owning object's dc_dist_object.
 *
 * \ref graphlab::fiber_remote_request is the version of this function
 * for remotely calling global functions.
 *
 * Example:
 * \code
 * // A print function is defined in the distributed object
 * class distributed_obj_example {
 *  graphlab::dc_dist_object<distributed_obj_example> rmi;
 *   ... initialization and constructor ...
 *  private:
 *    int add_one(int i) {
 *      return i + 1;
 *    }
 *  public:
 *    int add_one_from_machine_1(int i) {
 *      // calls the add_one function on machine 1 with the argument i
 *      // this call returns immediately
 *      graphlab::request_future<int> future =
 *          object_future_remote_request(rmi, 1, &distributed_obj_example::add_one, i);
 *
 *      // ... we can do other stuff here
 *      // then when we want the answer
 *      // this is safe to do within a fiber as it will not halt other fibers.
 *      int result = future();
 *      return result;
 *    }
 * }
 * \endcode
 *
 * \see graphlab::dc_dist_object::remote_request
 *      graphlab::dc_dist_object::future_remote_request
 *      graphlab::fiber_remote_request
 *
 * \param rmiobj The dc_dist_object to use to send the request.
 * \param targetmachine The ID of the machine to run the function on
 * \param fn The function to run on the target machine. Must be a pointer to
 *            member function in the owning object.
 * \param ... The arguments to send to Fn. Arguments must be serializable.
 *            and must be castable to the target types.
 *
 * \returns Returns a future templated around the same type as the return 
 *          value of the called function
 */
  request_future<RetVal> object_fiber_remote_request(dc_dist_object<T> rmiobj,
                                                     procid_t targetmachine, 
                                                     Fn fn, ...);



#endif



#include <boost/preprocessor.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>

#define GENARGS(Z,N,_)  BOOST_PP_CAT(T, N) BOOST_PP_CAT(i, N)
#define GENI(Z,N,_) BOOST_PP_CAT(i, N)
#define GENT(Z,N,_) BOOST_PP_CAT(T, N)
#define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);

#define REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  BOOST_PP_TUPLE_ELEM(1,0,ARGS) (procid_t target, \
                                 F remote_function BOOST_PP_COMMA_IF(N) \
                                 BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
  request_future<__GLRPC_FRESULT> reply(new fiber_reply_container);      \
  distributed_control* dc = distributed_control::get_instance(); \
  ASSERT_NE(dc, NULL); \
  dc->custom_remote_request(target, reply.get_handle(), STANDARD_CALL, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  return reply; \
} 

BOOST_PP_REPEAT(7, REQUEST_INTERFACE_GENERATOR, (request_future<__GLRPC_FRESULT> fiber_remote_request) )

#include <graphlab/rpc/function_arg_types_undef.hpp>

#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#define OBJECT_REQUEST_INTERFACE_GENERATOR(Z,N,ARGS) \
template<typename RMI, typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
  BOOST_PP_TUPLE_ELEM(1,0,ARGS) (RMI& rmi, \
                                 procid_t target, \
                                 F remote_function BOOST_PP_COMMA_IF(N) \
                                 BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
  request_future<__GLRPC_FRESULT> reply(new fiber_reply_container);      \
  rmi.custom_remote_request(target, reply.get_handle(), STANDARD_CALL, remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENI ,_) ); \
  return reply; \
} 



  /*
  Generates the interface functions. 3rd argument is a tuple (interface name, issue name, flags)
  */
BOOST_PP_REPEAT(7, OBJECT_REQUEST_INTERFACE_GENERATOR, (request_future<__GLRPC_FRESULT> object_fiber_remote_request) )

#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

#undef OBJECT_REQUEST_INTERFACE_GENERATOR
#undef REQUEST_INTERFACE_GENERATOR
#undef GENARC
#undef GENT
#undef GENI
#undef GENARGS

} // namespace graphlab

#endif
