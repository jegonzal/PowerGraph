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
This is a partial specialization of the remote_call_issue classes
described in function_call_issue.hpp and request_issue.hpp to handle the portable calls

\see function_call_issue.hpp
*/

#define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
#define GENARC(Z,N,_) arc << ( typename BOOST_PP_CAT(get_cleaned_user_arg, N)<F>::arg_type ) BOOST_PP_CAT(i, N);

#define PORTABLE_CALL_ISSUE_GENERATOR(Z,N,FNAME_AND_CALL) \
template<typename F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T) > \
  class BOOST_PP_CAT(FNAME_AND_CALL,N)<portable_call<F> BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> { \
   public: \
    static void exec(dc_send* sender, size_t flags, procid_t target, portable_call<F> remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)  ) {   \
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
    static typename function_ret_type<FRESULT>::type  exec(dc_send* sender, size_t flags, procid_t target, portable_call<F> remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)  ) {   \
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
    typename function_ret_type<FRESULT>::type  result; \
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
