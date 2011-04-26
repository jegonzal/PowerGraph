#ifndef GRAPHLAB_ROBUST_CAST_HPP
#define GRAPHLAB_ROBUST_CAST_HPP

#include <boost/utility.hpp>
#include <boost/type_traits/is_convertible.hpp>
namespace graphlab {
  /** robust_cast performs a static cast from type A to type B
      if a cast can be done. Return B() otherwise */
  
  template <typename Target, typename Source>
  typename boost::disable_if_c<boost::is_convertible<Source, Target>::value, 
                               Target>::type
                               robust_cast(const Source &h) {
    return Target();
  }
  
  template <typename Target, typename Source>
  typename boost::enable_if_c<boost::is_convertible<Source, Target>::value, 
                              Target>::type
                              robust_cast(const Source &h) {
    return h;
  }
}

#endif


