#ifndef FUNCTION_RETURN_TYPE_HPP
#define FUNCTION_RETURN_TYPE_HPP
#include <boost/preprocessor.hpp>
#include <graphlab/rpc/function_arg_types_def.hpp>
/**
Promotes a void to int so we can do request type calls on functions
with no return. 
function_ret_type<void>::type == int
and 
function_ret_type<T>::type == T for all other types T

mem_function_ret_type is similar but operates on member functions

*/
namespace graphlab {
namespace dc_impl {

template <typename RetType>
struct function_ret_type {
  typedef RetType type;
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(R, N)  BOOST_PP_CAT(i, N)
 
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
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(R, N) BOOST_PP_CAT(i, N)
 
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


template <typename RetType>
struct mem_function_ret_type {
  typedef RetType type;
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(R, N)  BOOST_PP_CAT(i, N)
 
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
  
  #define GENARGS(Z,N,_)  BOOST_PP_CAT(R, N) BOOST_PP_CAT(i, N)
 
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
