#ifndef GRAPHLAB_SERIALIZE_VECTOR_HPP
#define GRAPHLAB_SERIALIZE_VECTOR_HPP

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

#include <vector>
namespace graphlab {
  /**
    Serializes a vector
    Returns true on success, false on failure  */
  template <typename T>
  oarchive& operator<<(oarchive& a, const std::vector<T>& vec){
    serialize_iterator(a,vec.begin(), vec.end());
    return a;
  }

  /**
    deserializes a vector
    Returns true on success, false on failure  */
  template <typename T>
  iarchive& operator>>(iarchive& a, std::vector<T>& vec){
    vec.clear();
    deserialize_iterator<T>(a, std::inserter(vec, vec.end()));
    return a;
  }
} // namespace prl

#endif //PRL_SERIALIZE_VECTOR_HPP
