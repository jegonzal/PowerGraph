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


#ifndef OBJECT_BROADCAST_ISSUE_HPP
#define OBJECT_BROADCAST_ISSUE_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/object_call_dispatch.hpp>
#include <graphlab/rpc/object_call_issue.hpp>
#include <graphlab/rpc/is_rpc_call.hpp>
#include <graphlab/rpc/dc_thread_get_send_buffer.hpp>
#include <boost/preprocessor.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>

namespace graphlab{
namespace dc_impl {

/**
\ingroup rpc
\internal
\file object_broadcast_issue.hpp
 This is an internal function and should not be used directly

See object_call_issue.hpp for details. This is equivalent to the macro
expansion in object_call_issue with the difference that this takes an iterator 
sequence listing the machines to send to.

The code below generates the following for different number of arguments. Here, 
we demonstrate the 1 argument version.

\code
template < typename Iterator, typename T, typename F, typename T0 > 
class object_broadcast_issue1 {
 public:
  static void exec (dc_dist_object_base * rmi,
                    std::vector < dc_send * >sender, unsigned char flags,
                    Iterator target_begin, Iterator target_end, size_t objid,
                    F remote_function, const T0 & i0) {
    oarchive arc;
    arc.buf = (char *) malloc (65536);
    arc.len = 65536;
    size_t len =
      dc_send::write_packet_header (arc, _get_procid (), flags,
				    _get_sequentialization_key ());
    uint32_t beginoff = arc.off;
    dispatch_type d =
      dc_impl::OBJECT_NONINTRUSIVE_DISPATCH1 < distributed_control, T, F,
      T0 >;
    arc << reinterpret_cast < size_t > (d);
    serialize (arc, (char *) (&remote_function), sizeof (F));
    arc << objid;
    arc << i0;
    uint32_t curlen = arc.off - beginoff;
    *(reinterpret_cast < uint32_t * >(arc.buf + len)) = curlen;
    Iterator iter = target_begin;
    while (iter != target_end) {
      oarchive *buf = get_thread_local_buffer (*iter);
      buf->write (arc.buf, arc.off);
      release_thread_local_buffer (*iter, flags & CONTROL_PACKET);
      if ((flags & CONTROL_PACKET) == 0) {
        rmi->inc_bytes_sent ((*iter), curlen);
      }
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


#define REMOTE_BROADCAST_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename Iterator, typename T, typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
class  BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,0,FNAME_AND_CALL), N) { \
  public: \
  static void exec(dc_dist_object_base* rmi, std::vector<dc_send*> sender, unsigned char flags, \
                    Iterator target_begin, Iterator target_end, size_t objid, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    oarchive arc;       \
    arc.buf = (char*)malloc(INITIAL_BUFFER_SIZE); \
    arc.len = INITIAL_BUFFER_SIZE; \
    size_t len = dc_send::write_packet_header(arc, _get_procid(), flags, _get_sequentialization_key()); \
    uint32_t beginoff = arc.off; \
    dispatch_type d = BOOST_PP_CAT(dc_impl::OBJECT_NONINTRUSIVE_DISPATCH,N)<distributed_control,T,F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENT ,_) >;   \
    arc << reinterpret_cast<size_t>(d);                                 \
    serialize(arc, (char*)(&remote_function), sizeof(F));               \
    arc << objid;                                                       \
    BOOST_PP_REPEAT(N, GENARC, _)                                       \
    uint32_t curlen = arc.off - beginoff;   \
    *(reinterpret_cast<uint32_t*>(arc.buf + len)) = curlen; \
    Iterator iter = target_begin;                                       \
    while(iter != target_end) { \
      oarchive* buf = get_thread_local_buffer(*iter);  \
      buf->write(arc.buf, arc.off);  \
      release_thread_local_buffer(*iter, flags & CONTROL_PACKET); \
      if ((flags & CONTROL_PACKET) == 0) {                                 \
        rmi->inc_bytes_sent((*iter), curlen); \
      } \
      ++iter; \
    } \
    free(arc.buf); \
    if (flags & FLUSH_PACKET) pull_flush_soon_thread_local_buffer(); \
  }  \
};



/**
Generates a function call issue. 3rd argument is a tuple (issue name, dispacther name)
*/
BOOST_PP_REPEAT(7, REMOTE_BROADCAST_ISSUE_GENERATOR,  (object_broadcast_issue, _) )



#undef GENARC
#undef GENT
#undef GENI
#undef GENARGS
#undef REMOTE_BROADCAST_ISSUE_GENERATOR

} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

#endif

