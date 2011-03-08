#ifndef GRAPHLAB_SERIALIZE_UNORDERED_MAP_HPP
#define GRAPHLAB_SERIALIZE_UNORDERED_MAP_HPP

#include <boost/unordered_map.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

namespace graphlab {

namespace archive_detail {
  /** Serializes a map */
  template <typename ArcType, typename T, typename U>
  struct serialize_impl<ArcType, boost::unordered_map<T,U> > {
  static void exec(ArcType& a, const boost::unordered_map<T,U>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a map  */
      
  template <typename ArcType, typename T, typename U>
  struct deserialize_impl<ArcType, boost::unordered_map<T,U> > {
  static void exec(ArcType& a, boost::unordered_map<T,U>& vec){
    vec.clear();
    deserialize_iterator<ArcType, std::pair<T,U> >(a, std::inserter(vec,vec.end()));
  }
  };

} // archive_detail  
} // graphlab
#endif 
