#ifndef OBJECT_REQUEST_DISPATCH_HPP
#define OBJECT_REQUEST_DISPATCH_HPP
#include <sstream>
#include <iostream>
#include <string>
#include <functional>
#include <algorithm>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>

namespace graphlab {
namespace dc_impl{


/**
This is the dispatch function for the an object member function call.
This only does non-intrusive calls.
*/
#define GENFN(N) BOOST_PP_CAT(NIF, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENNIARGS(Z,N,_)  (BOOST_PP_CAT(f, N))
#define GENPARAMS(Z,N,_)  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; iarc >> (BOOST_PP_CAT(f, N)) ;
#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));

#define NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_) \
template<typename DcType,typename T, typename F  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(OBJECT_NONINTRUSIVE_REQUESTDISPATCH,N) (DcType& dc, procid_t source,  \
               std::istream &strm) { \
  iarchive iarc(strm); \
  F f; \
  deserialize(iarc, (char*)(&f), sizeof(F)); \
  size_t objid;   \
  iarc >> objid;  \
  T* obj = reinterpret_cast<T*>(dc.get_registered_object(objid)); \
  size_t id; iarc >> id;    \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  typename function_ret_type<FRESULT>::type ret = mem_function_ret_type<FRESULT>::BOOST_PP_CAT(fcall, N) \
                                              (f, obj BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENNIARGS ,_)); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
  boost::iostreams::stream<resizing_array_sink> retstrm(128);    \
  oarchive oarc(retstrm); \
  oarc << ret; \
  retstrm.flush(); \
  dc.fast_remote_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));\
} 

BOOST_PP_REPEAT(6, NONINTRUSIVE_DISPATCH_GENERATOR, _)


#undef GENFN
#undef GENFN2
#undef GENNIARGS
#undef GENPARAMS
#undef NONINTRUSIVE_DISPATCH_GENERATOR

} // namespace dc_impl
} // namespace graphlab
#include <graphlab/rpc/mem_function_arg_types_undef.hpp>
#endif
