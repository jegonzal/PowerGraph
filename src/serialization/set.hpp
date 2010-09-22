#ifndef GRAPHLAB_SERIALIZE_SET_HPP
#define GRAPHLAB_SERIALIZE_SET_HPP

#include <set>
#include <serialization/iarchive.hpp>
#include <serialization/oarchive.hpp>
#include <serialization/iterator.hpp>

namespace graphlab {
namespace archive_detail {
  /** serializes a set  */
  template <typename ArcType, typename T>
  struct serialize_impl<ArcType, std::set<T> > {
  static void exec(ArcType& a, const std::set<T>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a set  */
  template <typename ArcType, typename T>
  struct deserialize_impl<ArcType, std::set<T> > {
  static void exec(ArcType& a, std::set<T>& vec){
    vec.clear();
    deserialize_iterator<T>(a, std::inserter(vec,vec.end()));
  }
  };

} // archive_detail  
} // graphlab

#endif 
