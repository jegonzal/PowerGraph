#ifndef GRAPHLAB_SERIALIZE_UNORDERED_SET_HPP
#define GRAPHLAB_SERIALIZE_UNORDERED_SET_HPP

#include <boost/unordered_set.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

namespace graphlab {
namespace archive_detail {
  /** serializes a set  */
  template <typename ArcType, typename T>
  struct serialize_impl<ArcType, boost::unordered_set<T>, false > {
  static void exec(ArcType& a, const boost::unordered_set<T>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a set  */
  template <typename ArcType, typename T>
  struct deserialize_impl<ArcType, boost::unordered_set<T>, false > {
  static void exec(ArcType& a, boost::unordered_set<T>& vec){
    vec.clear();
    deserialize_iterator<ArcType, T>(a, std::inserter(vec,vec.end()));
  }
  };

} // archive_detail  
} // graphlab

#endif 
