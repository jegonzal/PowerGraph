/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef GRAPHLAB_OBJECT_PODCALL_DISPATCH_HPP
#define GRAPHLAB_OBJECT_PODCALL_DISPATCH_HPP
#include <iostream>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/rpc/dc_types.hpp>
#include <graphlab/rpc/dc_internal_types.hpp>
#include <graphlab/rpc/mem_function_arg_types_def.hpp>
#include <graphlab/rpc/pod_template_structs.hpp>
#include <boost/preprocessor.hpp>
namespace graphlab {
namespace dc_impl {




#define GENARGS(Z,N,_) (BOOST_PP_CAT(s->t, N))






#define CHARSTRINGFREE(Z,N,_)  charstring_free(BOOST_PP_CAT(f, N));


#define OBJECT_PODCALL_NONINTRUSIVE_DISPATCH_GENERATOR(Z,N,_)                   \
  template<typename DcType, typename T,                                 \
           typename F BOOST_PP_COMMA_IF(N)                              \
           BOOST_PP_ENUM_PARAMS(N, typename T) >                        \
  void BOOST_PP_CAT(OBJECT_PODCALL_NONINTRUSIVE_DISPATCH,N)(DcType& dc,         \
                                                    procid_t source,    \
                                                    unsigned char packet_type_mask, \
                                                    const char* data, \
                                                    size_t len){ \
    typedef BOOST_PP_CAT(pod_template_detail::pod_call_struct, N) <F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N, T)> pstruct; \
    const pstruct* s = reinterpret_cast<const pstruct*>(data);                 \
    T* obj = reinterpret_cast<T*>(dc.get_registered_object(s->objid));     \
    /* Deserialize the arguments to f */                                \
    /* Invoke f */                                                      \
    (obj->*(s->remote_function))(BOOST_PP_ENUM(N,GENARGS ,_)  );                           \
    /* Count the call if not a control call */                          \
    if ((packet_type_mask & CONTROL_PACKET) == 0)                       \
      dc.get_rmi_instance(s->objid)->inc_calls_received(source);           \
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
BOOST_PP_REPEAT(7, OBJECT_PODCALL_NONINTRUSIVE_DISPATCH_GENERATOR, _)



#undef GENARGS
#undef OBJECT_PODCALL_NONINTRUSIVE_DISPATCH_GENERATOR

} // namespace dc_impl
} // namespace graphlab


#include <graphlab/rpc/mem_function_arg_types_undef.hpp>
#endif

