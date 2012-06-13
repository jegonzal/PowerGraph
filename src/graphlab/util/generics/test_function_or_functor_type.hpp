#ifndef GRAPHLAB_TEST_FUNCTION_OR_FUNCTOR_TYPE_HPP
#define GRAPHLAB_TEST_FUNCTION_OR_FUNCTOR_TYPE_HPP
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_pointer.hpp>

namespace graphlab {


  template <typename F,
            typename PreferredFunctionForm,
            typename RetType,
            typename Arg1>
  struct test_function_or_functor_1 {

    // test if the functor type matches
    template <typename T, RetType (T::*)(Arg1)>
    struct SFINAE1 {};

    template <typename T>
    static char test1(SFINAE1<T, &T::operator()>*){}

    template <typename T>
    static int test1(...){}

    static const bool value = ((sizeof(test1<F>(0)) == sizeof(char)) ||
                               boost::is_same<F, PreferredFunctionForm>::value ||
                                boost::is_same<typename boost::remove_pointer<F>::type, PreferredFunctionForm>::value);
  };


  
  template <typename F,
            typename PreferredFunctionForm,
            typename RetType,
            typename Arg1>
  struct test_function_or_const_functor_1 {

    // test if the functor type matches
    template <typename T, RetType (T::*)(Arg1) const>
    struct SFINAE1 {};

    template <typename T>
    static char test1(SFINAE1<T, &T::operator()>*){}

    template <typename T>
    static int test1(...){}
    
    static const bool value = ((sizeof(test1<F>(0)) == sizeof(char)) ||
                               boost::is_same<F, PreferredFunctionForm>::value ||
                                boost::is_same<typename boost::remove_pointer<F>::type, PreferredFunctionForm>::value);
  };




  template <typename F,
            typename PreferredFunctionForm,
            typename RetType,
            typename Arg1,
            typename Arg2>
  struct test_function_or_functor_2 {

    // test if the functor type matches
    template <typename T, RetType (T::*)(Arg1, Arg2)>
    struct SFINAE1 {};

    template <typename T>
    static char test1(SFINAE1<T, &T::operator()>*){}

    template <typename T>
    static int test1(...){}

    static const bool value = ((sizeof(test1<F>(0)) == sizeof(char)) ||
                               boost::is_same<F, PreferredFunctionForm>::value ||
                                boost::is_same<typename boost::remove_pointer<F>::type, PreferredFunctionForm>::value);
  };

  

  template <typename F,
            typename PreferredFunctionForm,
            typename RetType,
            typename Arg1,
            typename Arg2>
  struct test_function_or_const_functor_2 {

    // test if the functor type matches
    template <typename T, RetType (T::*)(Arg1, Arg2) const>
    struct SFINAE1 {};

    template <typename T>
    static char test1(SFINAE1<T, &T::operator()>*){}

    template <typename T>
    static int test1(...){}

    static const bool value = ((sizeof(test1<F>(0)) == sizeof(char)) ||
                               boost::is_same<F, PreferredFunctionForm>::value ||
                                boost::is_same<typename boost::remove_pointer<F>::type, PreferredFunctionForm>::value);
  };

  
  
} // namespace graphlab
#endif