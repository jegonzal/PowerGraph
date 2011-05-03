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

// check for multiple inclusions
#ifdef __GLRPC_F0
#error "multiple includes of function arg types"
#endif

#include <boost/type_traits/decay.hpp>
#include <graphlab/util/generics/remove_member_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <boost/function.hpp>
// This is the member function version of function_arg_types
/*
A huge collection of useful typedefs.
F0... F5: identifies the arguments for an RPC aware function F. Dropping const and dropping references 
          (therefore allowing you to use F0....F5 to do casting.

__GLRPC_NIF0... __GLRPC_NIF5: identifies the arguments for an RPC unaware function F

__GLRPC_R0.... __GLRPC_R7: Identifies the actual arguments of the function F, without de-consting and de-reffing

__GLRPC_FRESULT: de-const and de-refed type of F's return type

__GLRPC_FARITY: the number of arguments F takes
*/
#define REMOVE_CONST_REF(REF) typename boost::remove_const<typename boost::remove_reference<REF>::type>::type



//#define F0 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg1_type)
//#define __GLRPC_FRESULT REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::result_type)

#define __GLRPC_NIF0 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg1_type)
#define __GLRPC_NIF1 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg2_type)
#define __GLRPC_NIF2 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg3_type)
#define __GLRPC_NIF3 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg4_type)
#define __GLRPC_NIF4 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg5_type)
#define __GLRPC_NIF5 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg6_type)
#define __GLRPC_NIF6 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg7_type)
#define __GLRPC_NIF7 REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::arg8_type)



#define __GLRPC_R0 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg1_type
#define __GLRPC_R1 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg2_type
#define __GLRPC_R2 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg3_type
#define __GLRPC_R3 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg4_type
#define __GLRPC_R4 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg5_type
#define __GLRPC_R5 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg6_type
#define __GLRPC_R6 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg7_type
#define __GLRPC_R7 typename boost::function<typename boost::remove_member_pointer<F>::type>::arg8_type

#define __GLRPC_FRESULT REMOVE_CONST_REF(typename boost::function<typename boost::remove_member_pointer<F>::type>::result_type)

#define __GLRPC_FARITY boost::function<typename boost::remove_member_pointer<F>::type>::arity

