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
#include <boost/preprocessor.hpp>

namespace graphlab {
namespace dc_impl {


#define GENARGS(Z,N,_)  BOOST_PP_CAT(const T, N) BOOST_PP_CAT(&i, N)
#define GENT(Z,N,_) BOOST_PP_CAT(T, N)
#define GENARC(Z,N,_) arc << BOOST_PP_CAT(i, N);

  
/**
\ingroup rpc_internal
\file

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
    static typename function_ret_type<
              typename boost::remove_const<
              typename boost::remove_reference<
              typename boost::function<
              typename boost::remove_member_pointer<F>
                ::type>::result_type>::type>::type>::type 
                exec(dc_send* sender, size_t flags, procid_t target,size_t objid, F remote_function , const T0 &i0 )
    {
        boost::iostreams::stream<resizing_array_sink> strm(128);
        oarchive arc(strm);
        reply_ret_type reply(1);
        dispatch_type d = dc_impl::OBJECT_NONINTRUSIVE_REQUESTDISPATCH1<distributed_control,T,F , T0 >;
        arc << reinterpret_cast<size_t>(d);
        serialize(arc, (char*)(&remote_function), sizeof(remote_function));
        arc << objid;
        arc << reinterpret_cast<size_t>(&reply);
        arc << i0;
        strm.flush();
        sender->send_data(target, flags, strm->str, strm->len);
        reply.wait();
        boost::iostreams::stream<boost::iostreams::array_source> retstrm(reply.val.c, reply.val.len);
        iarchive iarc(retstrm);
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
  static typename function_ret_type<__GLRPC_FRESULT>::type exec(dc_dist_object_base* rmi, dc_send* sender, size_t flags, procid_t target,size_t objid, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    boost::iostreams::stream<resizing_array_sink_ref> &strm = get_thread_local_stream();    \
    oarchive arc(strm);                         \
    reply_ret_type reply(REQUEST_WAIT_METHOD);      \
    dispatch_type d = BOOST_PP_CAT(dc_impl::OBJECT_NONINTRUSIVE_REQUESTDISPATCH,N)<distributed_control,T,F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENT ,_) >;  \
    arc << reinterpret_cast<size_t>(d);       \
    serialize(arc, (char*)(&remote_function), sizeof(remote_function)); \
    arc << objid;       \
    arc << reinterpret_cast<size_t>(&reply);       \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    strm.flush();           \
    sender->send_data(target, flags, strm->c_str(), strm->size());    \
    if ((flags & CONTROL_PACKET) == 0)                       \
      rmi->inc_bytes_sent(target, strm->size());           \
    reply.wait(); \
    boost::iostreams::stream<boost::iostreams::array_source> retstrm(reply.val.c, reply.val.len);  \
    iarchive iarc(retstrm);  \
    typename function_ret_type<__GLRPC_FRESULT>::type result; \
    iarc >> result;  \
    reply.val.free(); \
    return result;  \
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
