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


#ifndef REQUEST_ISSUE_HPP
#define REQUEST_ISSUE_HPP
#include <sstream>
#include <iostream>
#include <string>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/request_reply_handler.hpp>
#include <graphlab/rpc/request_future.hpp>
#include <graphlab/rpc/request_dispatch.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
#include <graphlab/rpc/dc_thread_get_send_buffer.hpp>
#include <graphlab/rpc/dc_compile_parameters.hpp>
#include <boost/preprocessor.hpp>

namespace graphlab {
namespace dc_impl {

/**

\internal
\ingroup rpc
\file request_issue.hpp

This is an internal function and should not be used directly.

This is an internal function and should not be used directly.
A request is an RPC which is performed "synchronously". The return value of the
function is returned.

The format of the RPC request is in the form of an archive and is as follows

The format of a "request" packet is in the form of an archive and is as follows

\li (dispatch_type*) -- pointer to target machine's dispatcher function
\li (void*)          -- pointer to target function
\li size_t           -- return ID
\li fn::arg1_type    -- target function's 1st argument
\li fn::arg2_type    -- target function's 2nd argument
\li  ...
\li fn::argN_type    -- target function's Nth argument


The ID here is a pointer to a ireply_container datastructure. When the remote machine completes
the function call, it will issue an RPC to the function reply_increment_counter on the originating machine.
The reply_increment_counter function  store the serialized return value in the ireply_container , as well
as perform an atomic increment on the ireply_container .

Here is an example of the marshall code for 1 argument
\code
namespace request_issue_detail {
  template < typename BoolType, typename F, typename T0 > 
  struct dispatch_selector1 {
    static dispatch_type dispatchfn () {
      return dc_impl::NONINTRUSIVE_REQUESTDISPATCH1 < distributed_control, F,
	T0 >;
    }
  };
  template < typename F, typename T0 > 
  struct dispatch_selector1 <boost::mpl::bool_ < true >, F, T0 > {
    static dispatch_type dispatchfn () {
      return dc_impl::REQUESTDISPATCH1 < distributed_control, F, T0 >;
    }
  };
}

template < typename F, typename T0 > 
class remote_request_issue1 {
 public:
  static void exec (dc_send * sender, size_t request_handle,
                    unsigned char flags, procid_t target, F remote_function,
                    const T0 & i0) {
    oarchive *ptr = get_thread_local_buffer (target);
    oarchive & arc = *ptr;
    size_t len =
      dc_send::write_packet_header (arc, _get_procid (), flags,
				    _get_sequentialization_key ());
    uint32_t beginoff = arc.off;
    dispatch_type d =
      request_issue_detail::dispatch_selector1 < typename is_rpc_call <
      F >::type, F, T0 >::dispatchfn ();
    arc << reinterpret_cast < size_t > (d);
    arc << reinterpret_cast < size_t > (remote_function);
    arc << request_handle;
    arc << i0;
    *(reinterpret_cast < uint32_t * >(arc.buf + len)) = arc.off - beginoff;
    release_thread_local_buffer (target, flags & CONTROL_PACKET);
    pull_flush_thread_local_buffer (target);
  }
};
\endcode

If the pointer to the dispatcher function is NULL, the next argument
will contain the name of the function. This is a "portable" call.
\see portable_issue.hpp
*/

#define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
#define GENT(Z,N,_) BOOST_PP_CAT(T, N)
#define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);


/**
The dispatch_selectorN structs are used to pick between the standard dispatcher and the nonintrusive dispatch
by checking if the function is a RPC style call or not.
*/
#define REMOTE_REQUEST_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
namespace request_issue_detail {      \
template <typename BoolType, typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
struct BOOST_PP_CAT(dispatch_selector, N){  \
  static dispatch_type dispatchfn() { return BOOST_PP_CAT(dc_impl::NONINTRUSIVE_REQUESTDISPATCH,N)<distributed_control,F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENT ,_) >; }  \
};\
template <typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
struct BOOST_PP_CAT(dispatch_selector, N)<boost::mpl::bool_<true>, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)>{  \
  static dispatch_type dispatchfn() { return BOOST_PP_CAT(dc_impl::REQUESTDISPATCH,N)<distributed_control,F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENT ,_) >; } \
}; \
}\
template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
class  BOOST_PP_CAT(FNAME_AND_CALL, N) { \
  public: \
  static void exec(dc_send* sender, size_t request_handle, unsigned char flags, procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    oarchive* ptr = get_thread_local_buffer(target);  \
    oarchive& arc = *ptr;                         \
    size_t len = dc_send::write_packet_header(arc, _get_procid(), flags, _get_sequentialization_key()); \
    uint32_t beginoff = arc.off; \
    dispatch_type d = BOOST_PP_CAT(request_issue_detail::dispatch_selector,N)<typename is_rpc_call<F>::type, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T) >::dispatchfn();   \
    arc << reinterpret_cast<size_t>(d);       \
    arc << reinterpret_cast<size_t>(remote_function); \
    arc << request_handle; \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    *(reinterpret_cast<uint32_t*>(arc.buf + len)) = arc.off - beginoff; \
    release_thread_local_buffer(target, flags & CONTROL_PACKET); \
    if (flags & FLUSH_PACKET) pull_flush_soon_thread_local_buffer(target); \
  }\
};


/**
Generates a function call issue. 3rd argument is the issue name
*/
BOOST_PP_REPEAT(7, REMOTE_REQUEST_ISSUE_GENERATOR,  remote_request_issue )



#undef GENARC
#undef GENT
#undef GENARGS
#undef REMOTE_REQUEST_ISSUE_GENERATOR


} // namespace dc_impl
} // namespace graphlab
#include <graphlab/rpc/function_arg_types_undef.hpp>

#endif

