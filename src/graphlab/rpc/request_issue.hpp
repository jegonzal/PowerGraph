#ifndef REQUEST_ISSUE_HPP
#define REQUEST_ISSUE_HPP
#include <sstream>
#include <iostream>
#include <string>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/request_dispatch.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>

namespace graphlab {
namespace dc_impl {

/**
A request is an RPC which is performed "synchronously". The return value of the
function is returned.

The format of the RPC request is in the form of an archive and is as follows
arc << [pointer to target machine's dispatcher function]
arc << [pointer to the actual function to call]
arc << [An "ID" that uniquely identifies this request]
arc << arg1 << arg2 << arg3 ...

The ID here is a pointer to a reply_ret_type datastructure. When the remote machien completes
the function call, it will issue an RPC to the function reply_increment_counter on the originating machine.
The reply_increment_counter function  store the serialized return value in the reply_ret_type, as well
as perform an atomic increment on the reply_ret_type.

Here is an example of the marshall code for 1 argument

namespace request_issue_detail
{
    template <typename BoolType, typename F , typename T0> struct dispatch_selector1
    {
        static dispatch_type dispatchfn()
        {
            return dc_impl::NONINTRUSIVE_REQUESTDISPATCH1<distributed_control,F , T0 >;
        }
    };
    template <typename F , typename T0> struct dispatch_selector1<boost::mpl::bool_<true>, F , T0>
    {
        static dispatch_type dispatchfn()
        {
            return dc_impl::REQUESTDISPATCH1<distributed_control,F , T0 >;
        }
    };
}


template<typename F , typename T0> class remote_request_issue1
{
    public: static typename function_ret_type<
                  typename boost::remove_const<
                  typename boost::remove_reference<
                  typename boost::function<
                  typename boost::remove_pointer<F>::type>::result_type>
                  ::type>::type>::type 
      exec(dc_send* sender, size_t flags, procid_t target, F remote_function , const T0 &i0 )
    {
        boost::iostreams::stream<resizing_array_sink> strm(128);
        oarchive arc(strm);
        reply_ret_type reply(1);
        dispatch_type d = request_issue_detail::dispatch_selector1<typename is_rpc_call<F>::type, F , T0 >::dispatchfn();
        arc << reinterpret_cast<size_t>(d);
        arc << reinterpret_cast<size_t>(remote_function);
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
              typename boost::remove_pointer<F>::type>::result_type>
              ::type>::type>::type result;
        iarc >> result;
        reply.val.free();
        return result;
    }
};





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
  static typename function_ret_type<FRESULT>::type exec(dc_send* sender, size_t flags, procid_t target, F remote_function BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_) ) {  \
    boost::iostreams::stream<resizing_array_sink> strm(128);    \
    oarchive arc(strm);                         \
    reply_ret_type reply(REQUEST_WAIT_METHOD);      \
    dispatch_type d = BOOST_PP_CAT(request_issue_detail::dispatch_selector,N)<typename is_rpc_call<F>::type, F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T) >::dispatchfn();   \
    arc << reinterpret_cast<size_t>(d);       \
    arc << reinterpret_cast<size_t>(remote_function); \
    arc << reinterpret_cast<size_t>(&reply);       \
    BOOST_PP_REPEAT(N, GENARC, _)                \
    strm.flush();           \
    sender->send_data(target, flags, strm->str, strm->len);    \
    reply.wait(); \
    boost::iostreams::stream<boost::iostreams::array_source> retstrm(reply.val.c, reply.val.len);    \
    iarchive iarc(retstrm);  \
    typename function_ret_type<FRESULT>::type result; \
    iarc >> result;  \
    reply.val.free(); \
    return result;  \
  }\
}; 


/**
Generates a function call issue. 3rd argument is the issue name
*/
BOOST_PP_REPEAT(6, REMOTE_REQUEST_ISSUE_GENERATOR,  remote_request_issue )



#undef GENARC
#undef GENT
#undef GENARGS
#undef REMOTE_REQUEST_ISSUE_GENERATOR
  
  
} // namespace dc_impl
} // namespace graphlab
#include <graphlab/rpc/function_arg_types_undef.hpp>

#endif
