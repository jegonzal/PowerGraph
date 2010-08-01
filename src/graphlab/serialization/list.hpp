#ifndef GRAPHLAB_SERIALIZE_LIST_HPP
#define GRAPHLAB_SERIALIZE_LIST_HPP

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

#include <list>

//TODO: Remove from the graphlab namespace

namespace graphlab {

  /** Serializes a list
      Returns true on success, false on failure  */
  template <typename T>
  oarchive& operator<<(oarchive& a, const std::list<T>& vec){
    serialize_iterator(a,vec.begin(), vec.end());
    return a;
  }

  /** deserializes a list
      Returns true on success, false on failure  */
  template <typename T>
  iarchive& operator>>(iarchive& a, std::list<T>& vec){
    vec.clear();
    deserialize_iterator<T>(a, std::inserter(vec, vec.end()));
    return a;
  }
} // namespace prl
#endif //PRL_SERIALIZE_LIST_HPP
