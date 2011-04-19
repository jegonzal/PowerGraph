#ifndef GRAPHLAB_FLOAT_SELECTOR_HPP
#define GRAPHLAB_FLOAT_SELECTOR_HPP

namespace graphlab {
  
  template <int len>
  struct float_selector {
    // invalid
  };


  template <>
  struct float_selector<4> {
    typedef float float_type;
  };

  template <>
  struct float_selector<8> {
    typedef double float_type;
  };

  template <>
  struct float_selector<16> {
    typedef long double float_type;
  };

}
#endif


