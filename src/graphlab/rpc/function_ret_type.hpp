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

#ifndef FUNCTION_RETURN_TYPE_HPP
#define FUNCTION_RETURN_TYPE_HPP
#include <boost/preprocessor.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
namespace graphlab {
namespace dc_impl {

  
/**
\ingroup rpc_internal
This struct performs two duties.
Firstly, it provides a consistent interface through a function called ::fcallN<F>
to complete a function call with a variable number of arguments.
Next, it provides the type of the return value of the function in ::type.
If the return type is void, it is promoted to an int. This makes the output
type of the function call be always serializable, simplifying the implementation
of "requests".
*/
template <typename RetType>
struct function_ret_type {
  typedef RetType type;
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(__GLRPC_R, N)  BOOST_PP_CAT(i, N)
 
  #define FCALL(Z, N, _) \
  template <typename F> \
  static RetType BOOST_PP_CAT(fcall, N)(F f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENARGS, _)){ \
    return f(BOOST_PP_ENUM_PARAMS(N, i)); \
  } 
    
  BOOST_PP_REPEAT(8, FCALL ,  _ )

  #undef FCALL
  #undef GENARGS

};

template <>
struct function_ret_type<void> {
  typedef size_t type;
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(__GLRPC_R, N) BOOST_PP_CAT(i, N)
 
  #define FCALL(Z, N, _) \
  template <typename F> \
  static size_t BOOST_PP_CAT(fcall, N)(F f BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENARGS, _)){ \
    f(BOOST_PP_ENUM_PARAMS(N, i)); \
    return 0; \
  } 
  
  BOOST_PP_REPEAT(8, FCALL ,  _ )

  #undef FCALL
  #undef GENARGS

};
#include <graphlab/rpc/function_arg_types_undef.hpp>


} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/mem_function_arg_types_def.hpp>

namespace graphlab {
namespace dc_impl {

/**
This struct performs two duties.
Firstly, it provides a consistent interface through a function called ::fcallN<F>
to complete a \b member function call with a variable number of arguments.
Next, it provides the type of the return value of the function in ::type.
If the return type is void, it is promoted to an int. This makes the output
type of the function call be always serializable, simplifying the implementation
of "requests".
*/
template <typename RetType>
struct mem_function_ret_type {
  typedef RetType type;
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(__GLRPC_R, N)  BOOST_PP_CAT(i, N)
 
  #define FCALL(Z, N, _) \
  template <typename F, typename T> \
  static RetType BOOST_PP_CAT(fcall, N)(F f , T t BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENARGS, _)){ \
    return (t->*f)(BOOST_PP_ENUM_PARAMS(N, i)); \
  }

  BOOST_PP_REPEAT(8, FCALL ,  _ )

  #undef FCALL
  #undef GENARGS

};

template <>
struct mem_function_ret_type<void> {
  typedef size_t type;
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(__GLRPC_R, N) BOOST_PP_CAT(i, N)
 
  #define FCALL(Z, N, _) \
  template <typename F, typename T> \
  static size_t BOOST_PP_CAT(fcall, N)(F f , T t BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM(N, GENARGS, _)){ \
     (t->*f)(BOOST_PP_ENUM_PARAMS(N, i)); \
     return 0; \
  }

  BOOST_PP_REPEAT(8, FCALL ,  _ )

  #undef FCALL
  #undef GENARGS

};




} // namespace dc_impl
} // namespace graphlab

#include <graphlab/rpc/mem_function_arg_types_undef.hpp>

#endif

