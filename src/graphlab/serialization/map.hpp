#ifndef GRAPHLAB_SERIALIZE_MAP_HPP
#define GRAPHLAB_SERIALIZE_MAP_HPP

#include <map>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

namespace graphlab {

  /** Serializes a map
      Returns true on success, false on failure  */
  template <typename T, typename U>
  oarchive& operator<<(oarchive& a, const std::map<T,U>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
    return a;
  }

  /** deserializes a map
      Returns true on success, false on failure  */
  template <typename T, typename U>
  iarchive& operator>>(iarchive& a, std::map<T,U>& vec){
    vec.clear();
    deserialize_iterator<std::pair<T,U> >(a, std::inserter(vec,vec.end()));
    return a;
  }
} // namespace prl

#endif //PRL_SERIALIZE_MAP_HPP
