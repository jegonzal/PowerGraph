#ifndef PORTABLE_DISPATCH_HPP
#define PORTABLE_DISPATCH_HPP
#include <string>
#include <boost/unordered_map.hpp>
#include <boost/preprocessor.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/function_call_dispatch.hpp>
#include <graphlab/rpc/request_dispatch.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/rpc/portable.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>

/**
\ingroup rpc_internal
\file
This is an internal function and should not be used directly.

The portable calls work in a slightly different way as the regular RPC system.
The receiving function of any portable call must be registered through the REGISTER_RPC
macro. 

This macro does a few things. Firstly, it instantiates a find_dispatcher
class. The job of the find_dispatcher class is to expand the number of
arguments and the return type and figure out what is the right dispatch function 
type to use.

If the return type is void, the find_dispatcher is partially specialized,
and will instantiate a PORTABLE_DISPATCH function, which will be returned in the operator()
of the find_dispatch class.

Otherwise, it instantiates a PORTABLE_REQUESTDISPATCH function which sends back the 
return value of the function.

Finally, the macro inserts the pointer to the dispatch function into a hashmap in the 
distributed_control class.
*/
namespace graphlab {


namespace dc_impl {

namespace portable_detail {
/**
Returns the dispatch function for a variety of functions
template arguments are
\tparam F the function type
\tparam Fret the function return type
\tparam Nargs the number of arguments of F
\tparam f the function itself
\tparam IsRPCCall whether it is an rpccall
*/
template<typename F, typename Fret, size_t Nargs, F f, typename IsRPCCall>
struct find_dispatcher{
  static void* dispatch_call_fn() { return NULL; }
  static void* dispatch_request_fn() { return NULL; }
};
};

#define GENFN(N) BOOST_PP_CAT(__GLRPC_F, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENARGS(Z,N,_)  BOOST_PP_CAT(__GLRPC_F, N)(BOOST_PP_CAT(f, N))
#define GENPARAMS(Z,N,_)  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; iarc >> (BOOST_PP_CAT(f, N)) ;
#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));

// for the non-intrusive variety
#define GENNIARGS(Z,N,_)  BOOST_PP_CAT(__GLRPC_NIF, N)(BOOST_PP_CAT(f, N))


#define PORTABLE_DISPATCH_GENERATOR(Z,N,_) \
template<typename DcType, typename F, F f  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(PORTABLEDISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask, \
               std::istream &strm) { \
  iarchive iarc(strm); \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  f(dc, source BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)  ); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
} \
\
template<typename DcType, typename F, F f  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(PORTABLE_REQUESTDISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask, \
               std::istream &strm) { \
  iarchive iarc(strm); \
  size_t id; iarc >> id;                        \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  typename function_ret_type<__GLRPC_FRESULT>::type ret = function_ret_type<__GLRPC_FRESULT>::BOOST_PP_CAT(fcall, BOOST_PP_ADD(N, 2))   \
                                                  (f, dc, source BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
  boost::iostreams::stream<resizing_array_sink> retstrm(128);    \
  oarchive oarc(retstrm); \
  oarc << ret; \
  retstrm.flush(); \
  if (packet_type_mask & CONTROL_PACKET) { \
    dc.control_call(source, PORTABLE(reply_increment_counter), id, blob(retstrm->str, retstrm->len));\
  } \
  else {  \
    dc.fast_remote_call(source, PORTABLE(reply_increment_counter), id, blob(retstrm->str, retstrm->len));\
  } \
} \
\
template<typename DcType, typename F, F f  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(PORTABLE_NONINTRUSIVE_DISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask, \
               std::istream &strm) { \
  iarchive iarc(strm); \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  f(BOOST_PP_ENUM(N,GENNIARGS ,_)  ); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
} \
\
template<typename DcType, typename F, F f  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(PORTABLE_NONINTRUSIVE_REQUESTDISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask,  \
               std::istream &strm) { \
  iarchive iarc(strm); \
  size_t id; iarc >> id;                        \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  typename function_ret_type<__GLRPC_FRESULT>::type ret = function_ret_type<__GLRPC_FRESULT>::BOOST_PP_CAT(fcall, N) \
                                          (f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENNIARGS ,_)); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
  boost::iostreams::stream<resizing_array_sink> retstrm(128);    \
  oarchive oarc(retstrm); \
  oarc << ret; \
  retstrm.flush(); \
  if (packet_type_mask & CONTROL_PACKET) { \
    dc.control_call(source, PORTABLE(reply_increment_counter), id, blob(retstrm->str, retstrm->len));\
  } \
  else {  \
    dc.fast_remote_call(source, PORTABLE(reply_increment_counter), id, blob(retstrm->str, retstrm->len));\
  } \
  free(retstrm->str);                                                 \
} 




BOOST_PP_REPEAT(BOOST_PP_INC(6), PORTABLE_DISPATCH_GENERATOR, _)

#undef GENFN
#undef GENFN2
#undef GENARGS
#undef GENPARAMS
#undef GENNIARGS
#undef DISPATCH_GENERATOR
/*
template<typename CommType, typename F> 
struct find_dispatcher<CommType, F, 0>{
  static void* dispatchfn() {
    return DISPATCH0<CommType, F>;
  }
};


template<typename CommType, typename F> 
struct find_dispatcher<CommType, F, 1>{
  static typename distributed_control<CommType>::dispatch_type dispatchfn() {
    return (DISPATCH1<CommType, F, F0>);
  }
};*/

// c++ will select the most specific template first
#define GENFN(Z,N,_) BOOST_PP_EXPAND(BOOST_PP_CAT(__GLRPC_F, N))
#define GENNIFN(Z,N,_) BOOST_PP_EXPAND(BOOST_PP_CAT(__GLRPC_NIF, N))

namespace portable_detail {

#define PORTABLE_FIND_DISPATCH_GENERATOR(Z, N, _) \
template<typename F, typename Fret, F f, typename IsRPCCall> \
struct find_dispatcher<F, Fret, N, f, IsRPCCall>{  \
  static dispatch_type dispatch_call_fn() { \
    return BOOST_PP_CAT(PORTABLE_NONINTRUSIVE_DISPATCH, N)<distributed_control, F, f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENNIFN, _)>; \
  } \
 static dispatch_type dispatch_request_fn() { \
    return BOOST_PP_CAT(PORTABLE_NONINTRUSIVE_REQUESTDISPATCH, N)<distributed_control, F, f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENNIFN, _)>; \
  } \
};  \
\
template<typename F, typename Fret, F f> \
struct find_dispatcher<F, Fret, BOOST_PP_ADD(N, 2), f, boost::mpl::bool_<true> >{  \
  static dispatch_type dispatch_call_fn() { \
    return BOOST_PP_CAT(PORTABLEDISPATCH, N)<distributed_control, F, f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENFN, _)>; \
  }  \
  static dispatch_type dispatch_request_fn() { \
    return BOOST_PP_CAT(PORTABLE_REQUESTDISPATCH, N)<distributed_control, F, f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENFN, _)>; \
  } \
}; 
BOOST_PP_REPEAT(6, PORTABLE_FIND_DISPATCH_GENERATOR, _)
#undef PORTABLE_FIND_DISPATCH_GENERATOR
#undef GENFN

}

} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/function_arg_types_undef.hpp>
#endif
