#ifndef OBJECT_CALL_DISPATCH_HPP
#define OBJECT_CALL_DISPATCH_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>
namespace graphlab {
namespace dc_impl {



/**
Object calls.
This is similar to a regular function call with the only difference that
it needs to locate the object using dc.get_registered_object(...)
After the function call, it also needs to increment the call count for the object context.

template<typename DcType,
        typename T, 
        typename F , 
        typename T0> 
        void OBJECT_NONINTRUSIVE_DISPATCH1 (DcType& dc, 
                                          procid_t source, 
                                          unsigned char packet_type_mask, 
                                          std::istream &strm)
{
    iarchive iarc(strm);
    F f;
    deserialize(iarc, (char*)(&f), sizeof(F));
    size_t objid;
    iarc >> objid;
    T* obj = reinterpret_cast<T*>(dc.get_registered_object(objid));
    T0 (f0) ;
    iarc >> (f0) ;
    (obj->*f)( (f0) );
    charstring_free(f0);
    if ((packet_type_mask & CONTROL_PACKET) == 0) dc.get_rmi_instance(objid)->inc_calls_received(source);
}

} 

*/



#define GENFN(N) BOOST_PP_CAT(NIF, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENARGS(Z,N,_) (BOOST_PP_CAT(f, N))
#define GENPARAMS(Z,N,_)  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; iarc >> (BOOST_PP_CAT(f, N)) ;
#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));


#define OBJECT_NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_) \
template<typename DcType,typename T, typename F  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(OBJECT_NONINTRUSIVE_DISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask,  \
               std::istream &strm) { \
  iarchive iarc(strm); \
  F f; \
  deserialize(iarc, (char*)(&f), sizeof(F)); \
  size_t objid;   \
  iarc >> objid;  \
  T* obj = reinterpret_cast<T*>(dc.get_registered_object(objid)); \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  (obj->*f)(BOOST_PP_ENUM(N,GENARGS ,_)  ); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
  if ((packet_type_mask & CONTROL_PACKET) == 0) dc.get_rmi_instance(objid)->inc_calls_received(source); \
} 

BOOST_PP_REPEAT(6, OBJECT_NONINTRUSIVE_DISPATCH_GENERATOR, _)

#undef GENFN
#undef GENFN2
#undef GENARGS
#undef GENPARAMS
#undef NONINTRUSIVE_DISPATCH_GENERATOR

} // namespace dc_impl
} // namespace graphlab


#include <graphlab/rpc/mem_function_arg_types_undef.hpp>
#endif

