#ifndef GRAPHLAB_MD5_HPP
#define GRAPHLAB_MD5_HPP

#include <graphlab/util/uint128.hpp>

namespace graphlab {
  /**
   * Computes the MD5 hash of a string
   */
  gl_uint128_t MD5(std::string str);
}

#endif