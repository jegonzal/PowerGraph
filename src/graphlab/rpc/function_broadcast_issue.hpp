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


#ifndef FUNCTION_BROADCAST_ISSUE_HPP
#define FUNCTION_BROADCAST_ISSUE_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/dc_thread_get_send_buffer.hpp>
#include <graphlab/rpc/function_call_dispatch.hpp>
#include <graphlab/rpc/function_call_issue.hpp>
#include <graphlab/rpc/is_rpc_call.hpp>
#include <boost/preprocessor.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>

namespace graphlab{
namespace dc_impl {

/**

\ingroup rpc
\internal
\file function_broadcast_issue.hpp

See function_call_issue.hpp for details. This is equivalent to the macro
expansion in remote_call_issue with the difference that this takes an iterator 
sequence listing the machines to send to.

The code below generates the following for different number of arguments. Here, 
we demonstrate the 1 argument version.
\code
template < typename Iterator, typename F, typename T0 > 
class remote_broadcast_issue1 {
 public:
  static void exec (std::vector < dc_send * >&sender, unsigned char flags,
                    Iterator target_begin, Iterator target_end,
                    F remote_function, const T0 & i0) {
      oarchive arc;
      arc.buf = (char *) malloc (65536);
      arc.len = 65536;
      size_t len =
        dc_send::write_packet_header (arc, _get_procid (), flags,
              _get_sequentialization_key ());
      uint32_t beginoff = arc.off;
      dispatch_type d =
        function_call_issue_detail::dispatch_selector1 < typename is_rpc_call <
        F >::type, F, T0 >::dispatchfn ();
      arc << reinterpret_cast < size_t > (d);
      arc << reinterpret_cast < size_t > (remote_function);
      arc << i0;
      *(reinterpret_cast < uint32_t * >(arc.buf + len)) = arc.off - beginoff;
      Iterator iter = target_begin;
      while (iter != target_end) {
        oarchive *buf = get_thread_local_buffer (*iter);
        buf->write (arc.buf, arc.off);
        release_thread_local_buffer (*iter, flags & CONTROL_PACKET);
        ++iter;
      }
      free (arc.buf);
    }
};
\endcode
*/



#define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
#define GENI(Z,N,_) BOOST_PP_CAT(i, N)
#define GENT(Z,N,_) BOOST_PP_CAT(T, N)
#define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);


/**
The dispatch_selectorN structs are used to pick between the standard dispatcher and the nonintrusive dispatch
by checking if the function is a RPC style call or not.
*/
#define REMOTE_BROADCAST_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename Iterator, typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
class  BOOST_PP_CAT(FNAME_AND_CALL, N) { \
  public: \
  static void exec(std::vector<dc_send*>& sender, unsigned char flags, Iterator target_begin, Iterator target_end, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    oarchive arc;       \
    arc.buf = (char*)malloc(INITIAL_BUFFER_SIZE); \
    arc.len = INITIAL_BUFFER_SIZE; \
    size_t len = dc_send::write_packet_header(arc, _get_procid(), flags, _get_sequentialization_key()); \
    uint32_t beginoff = arc.off; \
    dispatch_type d = BOOST_PP_CAT(function_call_issue_detail::dispatch_selector,N)<typename is_rpc_call<F>::type, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T) >::dispatchfn();   \
    arc << reinterpret_cast<size_t>(d);       \
    arc << reinterpret_cast<size_t>(remote_function); \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    *(reinterpret_cast<uint32_t*>(arc.buf + len)) = arc.off - beginoff; \
    Iterator iter = target_begin; \
    while(iter != target_end) { \
      oarchive* buf = get_thread_local_buffer(*iter);  \
      buf->write(arc.buf, arc.off);  \
      release_thread_local_buffer(*iter, flags & CONTROL_PACKET); \
      ++iter;    \
    } \
    free(arc.buf); \
    if (flags & FLUSH_PACKET) pull_flush_soon_thread_local_buffer(); \
  }\
};



/**
Generates a function call issue. 3rd argument is the issue name
*/
BOOST_PP_REPEAT(6, REMOTE_BROADCAST_ISSUE_GENERATOR,  remote_broadcast_issue )



#undef GENARC
#undef GENT
#undef GENI
#undef GENARGS
#undef REMOTE_BROADCAST_ISSUE_GENERATOR

} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/function_arg_types_undef.hpp>

#endif

