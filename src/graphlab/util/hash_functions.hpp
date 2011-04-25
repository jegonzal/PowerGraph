#ifndef GRAPHLAB_HASH_FUNCTIONS_HPP
#define GRAPHLAB_HASH_FUNCTIONS_HPP

namespace graphlab {
  /**
  \ingroup util_internal
  A hash function to maps a value to itself
  */
  class identity_hash{
  public:
    size_t operator()(const size_t &t) const{
      return t;
    }
  };
}
#endif
