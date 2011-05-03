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

#ifndef REPACK_DISPATCH_HPP
#define REPACK_DISPATCH_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>
namespace graphlab {
namespace dc_impl {

/**
\ingroup rpc_internal
\file
This is an internal function and should not be used directly

A "call" is an RPC which is performed asynchronously.
There are 2 types of calls. A "basic" call calls a standard C/C++ function
and does not require the function to be modified.
while the "regular" call requires the first 2 arguments of the function 
to be "distributed_control &dc, procid_t source".

A "dispatch" is a wrapper function on the receiving side of an RPC
which decodes the packet and performs the function call.

This scary looking piece of code is actually quite straightforward.
Given  function F, as well as input types T1 ... Tn
it will construct an input archive and deserialize the types T1.... Tn,
and call the function f with it. This code dispatches to the "intrusive" 
form of a function call (that is the function call must take a distributed_control
and a "procid_t source" as its first 2 arguments.

For instance, the 1 argument version of this is DISPATCH1:

\code
template<typename DcType, 
        typename F , 
        typename T0> void DISPATCH1 (DcType& dc, 
                                     procid_t source, 
                                     unsigned char packet_type_mask, 
                                     std::istream &strm)
{
    iarchive iarc(strm);
    size_t s;
    iarc >> s;
    F f = reinterpret_cast<F>(s);
    T0 (f0) ;
    iarc >> (f0) ;
    f(dc, source , (f0) );
    charstring_free(f0);
}
\endcode

charstring_free is a special template function which calls free(f1)
only if f1 is a character array (char*)

And similarly, the non-intrusive dispatch a little below

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
void BOOST_PP_CAT(DISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask, \
               std::istream &strm) { \
  iarchive iarc(strm); \
  size_t s; iarc >> s; F f = reinterpret_cast<F>(s); \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  f(dc, source BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N,GENARGS ,_)  ); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
} 

BOOST_PP_REPEAT(6, DISPATCH_GENERATOR, _)

#undef GENFN
#undef GENFN2
#undef GENARGS
#undef GENPARAMS
#undef DISPATCH_GENERATOR



/**
This is similar, but generates the non-intrusive version of a 
dispatcher. That is, the target function does not need to take
"distributed_control &dc, procid_t source" as its first 2 arguments.


template<typename DcType, 
        typename F , 
        typename T0> void NONINTRUSIVE_DISPATCH1 (DcType& dc, 
                                                procid_t source, 
                                                unsigned char packet_type_mask, 
                                                std::istream &strm)
{
    iarchive iarc(strm);
    size_t s;
    iarc >> s;
    F f = reinterpret_cast<F>(s);
    T0 (f0) ;
    iarc >> (f0) ;
    f( (f0) );
    charstring_free(f0);
}

*/



#define GENFN(N) BOOST_PP_CAT(__GLRPC_NIF, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENARGS(Z,N,_) (BOOST_PP_CAT(f, N))
#define GENPARAMS(Z,N,_)  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; iarc >> (BOOST_PP_CAT(f, N)) ;
#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));


#define NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_) \
template<typename DcType, typename F  BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, typename T)> \
void BOOST_PP_CAT(NONINTRUSIVE_DISPATCH,N) (DcType& dc, procid_t source, unsigned char packet_type_mask,  \
               std::istream &strm) { \
  iarchive iarc(strm); \
  size_t s; iarc >> s; F f = reinterpret_cast<F>(s); \
  BOOST_PP_REPEAT(N, GENPARAMS, _)                \
  f(BOOST_PP_ENUM(N,GENARGS ,_)  ); \
  BOOST_PP_REPEAT(N, CHARSTRINGFREE, _)                \
} 

BOOST_PP_REPEAT(6, NONINTRUSIVE_DISPATCH_GENERATOR, _)

#undef GENFN
#undef GENFN2
#undef GENARGS
#undef GENPARAMS
#undef NONINTRUSIVE_DISPATCH_GENERATOR

} // namespace dc_impl
} // namespace graphlab


#include <graphlab/rpc/function_arg_types_undef.hpp>
#endif

