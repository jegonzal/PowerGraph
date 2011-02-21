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
This is the dispatch function for the an object request.
This is similar to the standard request dispatcher, except that the object
needs to be located using the object id.

template<typename DcType,
        typename T, 
        typename F , 
        typename T0> 
        void OBJECT_NONINTRUSIVE_REQUESTDISPATCH1 (DcType& dc, 
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
    size_t id;
    iarc >> id;
    T0 (f0) ;
    iarc >> (f0) ;
    typename function_ret_type<
          typename boost::remove_const<
          typename boost::remove_reference<
          typename boost::function<
          typename boost::remove_member_pointer<F>::type>::result_type>
          ::type>::type>::type  
             ret = mem_function_ret_type<typename boost::remove_const<
                            typename boost::remove_reference<
                            typename boost::function<
                            typename boost::remove_member_pointer<F>::type>::result_type>
                            ::type>::type>::fcall1 (f, obj , (f0));
    charstring_free(f0);
    boost::iostreams::stream<resizing_array_sink> retstrm(128);
    oarchive oarc(retstrm);
    oarc << ret;
    retstrm.flush();
    if (packet_type_mask & CONTROL_PACKET)
    {
        dc.control_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));
    }
    else
    {
        dc.fast_remote_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));
    } if ((packet_type_mask & CONTROL_PACKET) == 0) dc.get_rmi_instance(objid)->inc_calls_received(source);
}



*/
#define GENFN(N) BOOST_PP_CAT(NIF, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENNIARGS(Z,N,_) (BOOST_PP_CAT(f, N))

#define GENPARAMS(Z,N,_)                                                \
  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ;                             \
  iarc >> (BOOST_PP_CAT(f, N)) ;

#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));

#define NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_)                          \
  template<typename DcType, typename T,                                 \
           typename F  BOOST_PP_COMMA_IF(N)                             \
           BOOST_PP_ENUM_PARAMS(N, typename T) >                        \
  void BOOST_PP_CAT(OBJECT_NONINTRUSIVE_REQUESTDISPATCH,N) (DcType& dc, \
                                                            procid_t source, \
                                                            unsigned char packet_type_mask, \
                                                            std::istream &strm) { \
    iarchive iarc(strm);                                                \
    F f;                                                                \
    deserialize(iarc, (char*)(&f), sizeof(F));                          \
    size_t objid;                                                       \
    iarc >> objid;                                                      \
    T* obj = reinterpret_cast<T*>(dc.get_registered_object(objid));     \
    size_t id; iarc >> id;                                              \
    BOOST_PP_REPEAT(N, GENPARAMS, _);                                   \
    typename function_ret_type<FRESULT>::type ret =                     \
      mem_function_ret_type<FRESULT>::BOOST_PP_CAT(fcall, N)            \
      (f, obj BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENNIARGS ,_));      \
    BOOST_PP_REPEAT(N, CHARSTRINGFREE, _);                              \
    boost::iostreams::stream<resizing_array_sink> retstrm(128);         \
    oarchive oarc(retstrm);                                             \
    oarc << ret;                                                        \
    retstrm.flush();                                                    \
    if (packet_type_mask & CONTROL_PACKET) {                            \
      dc.control_call(source,                                           \
                      reply_increment_counter,                          \
                      id,                                               \
                      blob(retstrm->str, retstrm->len));                \
    } else {                                                            \
      dc.fast_remote_call(source,                                       \
                          reply_increment_counter,                      \
                          id,                                           \
                          blob(retstrm->str, retstrm->len));            \
    }                                                                   \
    if ((packet_type_mask & CONTROL_PACKET) == 0) {                     \
      dc.get_rmi_instance(objid)->inc_calls_received(source);           \
      dc.get_rmi_instance(objid)->inc_bytes_sent(source, retstrm->len); \
    }                                                                   \
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
