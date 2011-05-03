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

#ifndef GRAPHLAB_OBJECT_CALL_DISPATCH_HPP
#define GRAPHLAB_OBJECT_CALL_DISPATCH_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <boost/preprocessor.hpp>
namespace graphlab {
namespace dc_impl {



/**
\ingroup rpc_internal
\file
This is an internal function and should not be used directly

This is similar to a regular function call in function_call_dispatch.hpp
with the only difference
that it needs to locate the object using dc.get_registered_object(...)
After the function call, it also needs to increment the call count for
the object context.
\code
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
    if ((packet_type_mask & CONTROL_PACKET) == 0) 
      dc.get_rmi_instance(objid)->inc_calls_received(source);
}

} 
\endcode
*/



#define GENFN(N) BOOST_PP_CAT(__GLRPC_NIF, N)
#define GENFN2(N) BOOST_PP_CAT(f, N)
#define GENARGS(Z,N,_) (BOOST_PP_CAT(f, N))

/**
 * This macro defines and deserializes each of the parameters to the
 * function.
 */
#define GENPARAMS(Z,N,_)  \
  BOOST_PP_CAT(T, N) (BOOST_PP_CAT(f, N)) ; \
  iarc >> (BOOST_PP_CAT(f, N)) ;

#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));


#define OBJECT_NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_)                   \
  template<typename DcType, typename T,                                 \
           typename F BOOST_PP_COMMA_IF(N)                              \
           BOOST_PP_ENUM_PARAMS(N, typename T) >                        \
  void BOOST_PP_CAT(OBJECT_NONINTRUSIVE_DISPATCH,N)(DcType& dc,         \
                                                    procid_t source,    \
                                                    unsigned char packet_type_mask, \
                                                    std::istream& strm){ \
    iarchive iarc(strm);                                                \
    F f;                                                                \
    deserialize(iarc, (char*)(&f), sizeof(F));                          \
    size_t objid;                                                       \
    iarc >> objid;                                                      \
    T* obj = reinterpret_cast<T*>(dc.get_registered_object(objid));     \
    /* Deserialize the arguments to f */                                \
    BOOST_PP_REPEAT(N, GENPARAMS, _);                                   \
    /* Invoke f */                                                      \
    (obj->*f)(BOOST_PP_ENUM(N,GENARGS ,_)  );                           \
    /* Free the buffers for the args */                                 \
    BOOST_PP_REPEAT(N, CHARSTRINGFREE, _) ;                             \
    /* Count the call if not a control call */                          \
    if ((packet_type_mask & CONTROL_PACKET) == 0)                       \
      dc.get_rmi_instance(objid)->inc_calls_received(source);           \
  } 



/**
 * This macro generates dispatch functions for functions for rpc calls
 * with up to 6 arguments.
 *
 * Remarks: If the compiler generates the following error "Too
 * few/many arguments to function" at this point is is due to the
 * caller not providing the correct number fo arguments in the RPC
 * call.  Note that default arguments are NOT supported in rpc calls
 * and so all arguments must be provided.
 *
 */
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

