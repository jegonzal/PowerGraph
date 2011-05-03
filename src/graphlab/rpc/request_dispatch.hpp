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

#ifndef REQUEST_DISPATCH_HPP
#define REQUEST_DISPATCH_HPP
#include <sstream>
#include <iostream>
#include <string>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/reply_increment_counter.hpp>
#include <graphlab/rpc/function_ret_type.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>

namespace graphlab {
namespace dc_impl{

/**
\ingroup rpc_internal
\file

This is an internal function and should not be used directly.

Given  function F, as well as input types T1 ... Tn
it will construct an input archive and deserialize the types T1.... Tn,
and call the function f with it. The return value of the function
is then returned to the caller through the fast call to the 
source's reply_increment_counter. This code dispatches to the "intrusive" 
form of a function call (that is the function call must take a distributed_control
and a "procid_t source" as its first 2 arguments.

For instance, the 1 argument of this will be DISPATCH1:
\code
template<typename DcType, 
    typename F , 
    typename T0> void REQUESTDISPATCH1 (DcType& dc, 
                                        procid_t source, 
                                        unsigned char packet_type_mask, 
                                        std::istream &strm)
{
    iarchive iarc(strm);
    size_t s;
    iarc >> s;
    F f = reinterpret_cast<F>(s);
    size_t id;
    iarc >> id;
    T0 (f0) ;
    iarc >> (f0) ;
    typename function_ret_type<
        typename boost::remove_const<
        typename boost::remove_reference<
        typename boost::function<
        typename boost::remove_pointer<F>::type>
        ::result_type>::type>::type>::type 
        ret = function_ret_type<
                    typename boost::remove_const<
                    typename boost::remove_reference
                    <typename boost::function<
                    typename boost::remove_pointer<F>::type>
                    ::result_type>::type>::type>::fcall3 (f, dc, source , (f0));
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
    }
}
\endcode

charstring_free is a special template function which calls free(f1)
only if f1 is a character array (char*)

Note that the template around DcType is *deliberate*. This prevents this
function from instantiating the distributed_control until as late as possible, 
avoiding problems with circular references.

*/
#define GENFN(N) BOOST_PP_CAT(__GLRPC_F, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENARGS(Z,N,_)  (BOOST_PP_CAT(f, N))
#define GENPARAMS(Z,N,_)  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; iarc >> (BOOST_PP_CAT(f, N)) ;
#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));


#define DISPATCH_GENERATOR(Z,N,_) \
template<typename DcType, typename F  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(REQUESTDISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask, \
               std::istream &strm) { \
  iarchive iarc(strm); \
  size_t s; iarc >> s; F f = reinterpret_cast<F>(s); \
  size_t id; iarc >> id;    \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  typename function_ret_type<__GLRPC_FRESULT>::type ret = function_ret_type<__GLRPC_FRESULT>::BOOST_PP_CAT(fcall, BOOST_PP_ADD(N, 2))   \
                                                  (f, dc, source BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
  boost::iostreams::stream<resizing_array_sink> retstrm(128);    \
  oarchive oarc(retstrm); \
  oarc << ret; \
  retstrm.flush(); \
  if (packet_type_mask & CONTROL_PACKET) { \
    dc.control_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));\
  } \
  else {  \
    dc.fast_remote_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));\
  } \
  free(retstrm->str);                                                 \
} 

BOOST_PP_REPEAT(6, DISPATCH_GENERATOR, _)

#undef GENFN
#undef GENFN2
#undef GENARGS
#undef GENPARAMS
#undef DISPATCH_GENERATOR

/**
Same as above, but is the non-intrusive version.
*/
#define GENFN(N) BOOST_PP_CAT(NIF, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENNIARGS(Z,N,_)  (BOOST_PP_CAT(f, N))
#define GENPARAMS(Z,N,_)  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; iarc >> (BOOST_PP_CAT(f, N)) ;
#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));

#define NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_) \
template<typename DcType, typename F  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(NONINTRUSIVE_REQUESTDISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask, \
               std::istream &strm) { \
  iarchive iarc(strm); \
  size_t s; iarc >> s; F f = reinterpret_cast<F>(s); \
  size_t id; iarc >> id;    \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  typename function_ret_type<__GLRPC_FRESULT>::type ret = function_ret_type<__GLRPC_FRESULT>::BOOST_PP_CAT(fcall, N) \
                                          (f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENNIARGS ,_)); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
  boost::iostreams::stream<resizing_array_sink> retstrm(128);    \
  oarchive oarc(retstrm); \
  oarc << ret; \
  retstrm.flush(); \
  if (packet_type_mask & CONTROL_PACKET) { \
    dc.control_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));\
  } \
  else {  \
    dc.fast_remote_call(source, reply_increment_counter, id, blob(retstrm->str, retstrm->len));\
  } \
  free(retstrm->str);                                                 \
} 

BOOST_PP_REPEAT(6, NONINTRUSIVE_DISPATCH_GENERATOR, _)


#undef GENFN
#undef GENFN2
#undef GENNIARGS
#undef GENPARAMS
#undef NONINTRUSIVE_DISPATCH_GENERATOR

} // namespace dc_impl
} // namespace graphlab
#include <graphlab/rpc/function_arg_types_undef.hpp>
#endif
