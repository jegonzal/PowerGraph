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


#ifndef OBJECT_REQUEST_ISSUE_HPP
#define OBJECT_REQUEST_ISSUE_HPP
#include <sstream>
#include <iostream>
#include <string>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/object_request_dispatch.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <graphlab/rpc/request_future.hpp>
#include <graphlab/rpc/archive_memory_pool.hpp>
#include <graphlab/rpc/dc_compile_parameters.hpp>
#include <boost/preprocessor.hpp>

namespace graphlab {
namespace dc_impl {


#define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
#define GENT(Z,N,_) BOOST_PP_CAT(T, N)
#define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);


/**
\internal
\ingroup rpc
\file object_request_issue.hpp


This is an internal function and should not be used directly

This is the marshall function for the an object member function call.
This is very similar to the standard function request issue in request_issue.hpp
, with the only difference that an object id has to be transmitted

An annoyingly long sequence of type manipulations are needed to identify the return type.
\code
template<typename T,
        typename F ,
        typename T0> class object_request_issue1
{
    public:
    static request_future<typename function_ret_type<
              typename boost::remove_const<
              typename boost::remove_reference<
              typename boost::function<
              typename boost::remove_member_pointer<F>
                ::type>::result_type>::type>::type>::type (void)>
                exec(dc_send* sender, unsigned char flags, procid_t target,size_t objid, F remote_function , const T0 &i0 )
    {
        oarchive arc;
        arc.advance(sizeof(packet_hdr));
        reply_ret_type reply(1);
        dispatch_type d = dc_impl::OBJECT_NONINTRUSIVE_REQUESTDISPATCH1<distributed_control,T,F , T0 >;
        arc << reinterpret_cast<size_t>(d);
        serialize(arc, (char*)(&remote_function), sizeof(remote_function));
        arc << objid;
        arc << reinterpret_cast<size_t>(reply.reply.get());
        arc << i0;
        sender->send_data(target, flags, arc.buf, arc.off);
        reply.wait();
        iarchive iarc(reply.val.c, reply.val.len);
        typename function_ret_type<
            typename boost::remove_const<
            typename boost::remove_reference<
            typename boost::function<
            typename boost::remove_member_pointer<F>
                ::type>::result_type>::type>::type>::type result;
        iarc >> result;
        reply.val.free();
        return result;
    }
};
\endcode


*/
#define REMOTE_REQUEST_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename T,typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
class  BOOST_PP_CAT(FNAME_AND_CALL, N) { \
  public: \
  static request_future<__GLRPC_FRESULT> exec(dc_dist_object_base* rmi, dc_send* sender, unsigned char flags, procid_t target,size_t objid, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    oarchive* ptr = oarchive_from_pool();       \
    oarchive& arc = *ptr;                         \
    arc.advance(sizeof(size_t) + sizeof(packet_hdr));            \
    request_future<__GLRPC_FRESULT> reply;   \
    dispatch_type d = BOOST_PP_CAT(dc_impl::OBJECT_NONINTRUSIVE_REQUESTDISPATCH,N)<distributed_control,T,F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENT ,_) >;  \
    arc << reinterpret_cast<size_t>(d);       \
    serialize(arc, (char*)(&remote_function), sizeof(remote_function)); \
    arc << objid;       \
    arc << reinterpret_cast<size_t>(reply.reply.get());       \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    if (arc.off >= BUFFER_RELINQUISH_LIMIT) {  \
      sender->send_data(target,flags , arc.buf, arc.off);    \
      arc.buf = NULL; arc.len = 0;   \
    } else {        \
      char* newbuf = (char*)malloc(arc.off); memcpy(newbuf, arc.buf, arc.off); \
      sender->send_data(target,flags , newbuf, arc.off);    \
    }     \
    release_oarchive_to_pool(ptr); \
    if ((flags & CONTROL_PACKET) == 0)                       \
      rmi->inc_bytes_sent(target, arc.off);           \
    return reply;   \
  }\
};

BOOST_PP_REPEAT(6, REMOTE_REQUEST_ISSUE_GENERATOR,  object_request_issue )



#undef GENARC
#undef GENT
#undef GENARGS
#undef REMOTE_REQUEST_ISSUE_GENERATOR


} // namespace dc_impl
} // namespace graphlab
#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

#endif

