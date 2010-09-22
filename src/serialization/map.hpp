#ifndef GRAPHLAB_SERIALIZE_MAP_HPP
#define GRAPHLAB_SERIALIZE_MAP_HPP

#include <map>
#include <serialization/iarchive.hpp>
#include <serialization/oarchive.hpp>
#include <serialization/iterator.hpp>

namespace graphlab {

namespace archive_detail {
  /** Serializes a map */
  template <typename ArcType, typename T, typename U>
  struct serialize_impl<ArcType, std::map<T,U> > {
  static void exec(ArcType& a, const std::map<T,U>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a map  */
      
  template <typename ArcType, typename T, typename U>
  struct deserialize_impl<ArcType, std::map<T,U> > {
  static void exec(ArcType& a, std::map<T,U>& vec){
    vec.clear();
    deserialize_iterator<std::pair<T,U> >(a, std::inserter(vec,vec.end()));
  }
  };

} // archive_detail  
} // graphlab
#endif 
