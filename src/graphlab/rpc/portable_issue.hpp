/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PORTABLE_ISSUE_HPP
#define PORTABLE_ISSUE_HPP
#include <graphlab/rpc/function_call_issue.hpp>
#include <graphlab/rpc/request_issue.hpp>
#include <graphlab/rpc/portable.hpp>
#include <graphlab/rpc/dc_send.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/is_rpc_call.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>
namespace graphlab{
namespace dc_impl {

/**
\ingroup rpc_internal
 \file
 
 This is an internal function and should not be used directly
 
This is a partial specialization of the remote_call_issue classes
described in function_call_issue.hpp and request_issue.hpp to handle the portable calls
A portable call permits function calls across different binaries, but required
registration.

Portable Call Formats \n
===================== \n
The format of a "portable call" packet is in the form of an archive and is as follows

\li size_t(NULL)     -- NULL. Corresponds to the dispatch_type* in the native call
\li  std::string      -- name of the function to call
\li  char(0)          -- flag that this is a call
\li  fn::arg1_type    -- target function's 1st argument
\li  fn::arg2_type    -- target function's 2nd argument
\li ...
\li  fn::argN_type    -- target function's Nth argument

Unlike the native call, argument casting is performed by the caller. The caller is
required to know the type of the function (At least through a function prototype)

---------\n
The format of a "portable request" packet is in the form of an archive and is as follows
\li  size_t(NULL)     -- NULL. Corresponds to the dispatch_type* in the native call
\li  std::string      -- name of the function to call
\li  char(1)          -- flag that this is a request
\li  fn::arg1_type    -- target function's 1st argument
\li  fn::arg2_type    -- target function's 2nd argument
\li   ...
\li  fn::argN_type    -- target function's Nth argument

Unlike the native call, argument casting is performed by the caller. The caller is
required to know the type of the function (At least through a function prototype)
At the end of the request, the dispatch will perform a fast call to the
 reply_increment_counter on the source machine passing the ID as an argument.
 The return data is passed using the actual result type of the function.

\see function_call_issue.hpp
*/

#define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
#define GENARC(Z,N,_) arc << ( typename BOOST_PP_CAT(get_cleaned_user_arg, N)<F>::arg_type ) BOOST_PP_CAT(i, N);

#define PORTABLE_CALL_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T) > \
  class BOOST_PP_CAT(FNAME_AND_CALL,N)<portable_call<F> BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> { \
   public: \
    static void exec(dc_send* sender, unsigned char flags, procid_t target, portable_call<F> remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)  ) {   \
    boost::iostreams::stream<resizing_array_sink_ref> &strm = get_thread_local_stream();    \
      oarchive arc(strm);                           \
      arc << 0;       \
      arc << remote_function.fname;      \
      arc << char(0);    \
      BOOST_PP_REPEAT(N, GENARC, _)                \
        strm.flush();                   \
      sender->send_data(target,  flags, strm->c_str(), strm->size());    \
    }  \
  };


BOOST_PP_REPEAT(6, PORTABLE_CALL_ISSUE_GENERATOR, remote_call_issue)
//BOOST_PP_REPEAT(6, PORTABLE_CALL_ISSUE_GENERATOR, basic_call_issue)


#define PORTABLE_REQUEST_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T) > \
  class BOOST_PP_CAT(FNAME_AND_CALL,N)<portable_call<F> BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> { \
   public: \
    static typename function_ret_type<__GLRPC_FRESULT>::type  exec(dc_send* sender, unsigned char flags, procid_t target, portable_call<F> remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)  ) {   \
    boost::iostreams::stream<resizing_array_sink_ref> &strm = get_thread_local_stream();    \
    oarchive arc(strm);                         \
    reply_ret_type reply(REQUEST_WAIT_METHOD);      \
    size_t fn = 0; \
    arc << fn;       \
    arc << remote_function.fname; \
    arc << char(1);       \
    arc << reinterpret_cast<size_t>(&reply);       \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    strm.flush();           \
    sender->send_data(target, flags, strm->c_str(), strm->size());    \
    reply.wait(); \
    boost::iostreams::stream<boost::iostreams::array_source> retstrm(reply.val.c, reply.val.len);    \
    iarchive iarc(retstrm);  \
    typename function_ret_type<__GLRPC_FRESULT>::type  result; \
    iarc >> result;  \
    reply.val.free(); \
    return result;  \
    }  \
  };

BOOST_PP_REPEAT(6, PORTABLE_REQUEST_ISSUE_GENERATOR, remote_request_issue)

#undef GENARC
#undef GENARGS
#undef PORTABLE_CALL_ISSUE_GENERATOR
#undef PORTABLE_REQUEST_ISSUE_GENERATOR

#include <graphlab/rpc/function_arg_types_undef.hpp>

} // namespace dc_impl
} // namespace graphlab

#endif

