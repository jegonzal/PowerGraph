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


#ifndef OBJECT_CALL_ISSUE_HPP
#define OBJECT_CALL_ISSUE_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/rpc/object_call_dispatch.hpp>
#include <graphlab/rpc/is_rpc_call.hpp>
#include <graphlab/rpc/dc_thread_get_send_buffer.hpp>
#include <boost/preprocessor.hpp>
#include <graphlab/rpc/dc_compile_parameters.hpp>
#include <graphlab/util/generics/blob.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>

namespace graphlab{
namespace dc_impl {

/**
\ingroup rpc
\internal
\file object_call_issue.hpp
 This is an internal function and should not be used directly

Marshalls a object function call to a remote machine.
This is similar to the regular function call in function_call_issue.hpp
with the only difference that the object id needs to be transmitted as well.

\code
template < typename T, typename F, typename T0 > 
class object_call_issue1 {
 public:
  static void exec (dc_dist_object_base * rmi, dc_send * sender,
                    unsigned char flags, procid_t target, size_t objid,
                    F remote_function, const T0 & i0) {
    oarchive *ptr = get_thread_local_buffer (target);
    oarchive & arc = *ptr;
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
    release_thread_local_buffer (target, flags & CONTROL_PACKET);
    if ((flags & CONTROL_PACKET) == 0) {
      rmi->inc_bytes_sent (target, curlen);
    }
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
#define REMOTE_CALL_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename T, typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
class  BOOST_PP_CAT(BOOST_PP_TUPLE_ELEM(2,0,FNAME_AND_CALL), N) { \
  public: \
  static void exec(dc_dist_object_base* rmi, dc_send* sender, unsigned char flags, procid_t target, size_t objid, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    oarchive* ptr = get_thread_local_buffer(target);  \
    oarchive& arc = *ptr;                         \
    size_t len = dc_send::write_packet_header(arc, _get_procid(), flags, _get_sequentialization_key()); \
    uint32_t beginoff = arc.off; \
    dispatch_type d = BOOST_PP_CAT(dc_impl::OBJECT_NONINTRUSIVE_DISPATCH,N)<distributed_control,T,F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENT ,_) >;   \
    arc << reinterpret_cast<size_t>(d);       \
    serialize(arc, (char*)(&remote_function), sizeof(F)); \
    arc << objid;       \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    uint32_t curlen = arc.off - beginoff;   \
    *(reinterpret_cast<uint32_t*>(arc.buf + len)) = curlen; \
    release_thread_local_buffer(target, flags & CONTROL_PACKET); \
    if ((flags & CONTROL_PACKET) == 0) {                      \
      rmi->inc_bytes_sent(target, curlen);           \
    } \
    if (flags & FLUSH_PACKET) pull_flush_soon_thread_local_buffer(target); \
  } \
  \
};


/**
 * \ingroup rpc
 * \internal
 *
 * This generates a "split call". Where the header of the call message
 * is written to with split_call_begin, and the message actually sent with
 * split_call_end(). It is then up to the user to serialize the message arguments
 * into the oarchive returned. The split call can provide performance gains 
 * when the contents of the message are large, since this allows the user to
 * control the serialization process. For examples, see 
 * \ref dc_dist_object::split_call_begin
 */
template <typename T, typename F>
class object_split_call {
 public:
  static oarchive* split_call_begin(dc_dist_object_base* rmi, size_t objid, F remote_function) {
    oarchive* ptr = new oarchive;
    oarchive& arc = *ptr;
    arc.buf = (char*)malloc(INITIAL_BUFFER_SIZE); 
    arc.len = INITIAL_BUFFER_SIZE; 
    arc.advance(sizeof(packet_hdr));
    dispatch_type d = dc_impl::OBJECT_NONINTRUSIVE_DISPATCH2<distributed_control,T,F,size_t, wild_pointer>;
    arc << reinterpret_cast<size_t>(d);
    serialize(arc, (char*)(&remote_function), sizeof(F));
    arc << objid;
    // make a gap for the blob size argument
    // write the largest possible size_t. That will allow it to bypass
    // dynamic length encoding issues.
    // patch the header with the offset to this point. 
    (*reinterpret_cast<size_t*>(arc.buf)) = arc.off + 1;
    arc << (size_t)(-1);
    return ptr;
  }
  static void split_call_cancel(oarchive* oarc) {
    free(oarc->buf);
    delete oarc;
  }

/**
 * \ingroup rpc
 * \internal
 *
 * This sends a message first created with split_call_begin. The archive
 * pointer is consumed.
 */
  static void split_call_end(dc_dist_object_base* rmi,
                             oarchive* oarc, dc_send* sender, procid_t target, unsigned char flags) {
    // header points to the location of the blob size argument
    size_t blobsize_offset = *reinterpret_cast<size_t*>(oarc->buf);
    (*reinterpret_cast<size_t*>(oarc->buf + blobsize_offset)) = oarc->off - blobsize_offset - sizeof(size_t);
    // write the packet header
    packet_hdr* hdr = reinterpret_cast<packet_hdr*>(oarc->buf);
    hdr->len = oarc->off - sizeof(packet_hdr);
    hdr->src = _get_procid();
    hdr->packet_type_mask = flags;
    hdr->sequentialization_key = _get_sequentialization_key();
    size_t len = hdr->len;
    write_thread_local_buffer(target, oarc->buf, oarc->off, flags & CONTROL_PACKET);
    if ((flags & CONTROL_PACKET) == 0) {
      rmi->inc_bytes_sent(target, len);
    }
    if (flags & FLUSH_PACKET) pull_flush_soon_thread_local_buffer(target); 
    delete oarc;
  }
};

/**
Generates a function call issue. 3rd argument is a tuple (issue name, dispacther name)
*/

BOOST_PP_REPEAT(7, REMOTE_CALL_ISSUE_GENERATOR,  (object_call_issue, _) )


#undef GENARC
#undef GENT
#undef GENI
#undef GENARGS
#undef REMOTE_CALL_ISSUE_GENERATOR

} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

#endif

