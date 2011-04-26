#ifndef GRAPHLAB_IS_POD_HPP
#define GRAPHLAB_IS_POD_HPP
#include <boost/type_traits.hpp>

namespace graphlab {

  template <typename T>
  struct gl_is_pod{
    // it is a pod and is not an integer since we have special handlings for integers

    // (T is POD and  T is not an integer of size >= 2)
    BOOST_STATIC_CONSTANT(bool, value =
                          (
                           boost::type_traits::ice_and<
                             boost::is_pod<T>::value,
                             boost::type_traits::ice_not<
                               boost::type_traits::ice_and<
                                 boost::is_integral<T>::value,
                                 sizeof(T) >= 2
                                 >::value
                               >::value
                             >::value
                          ));

  };

}

#endif



