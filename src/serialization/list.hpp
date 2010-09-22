#ifndef GRAPHLAB_SERIALIZE_LIST_HPP
#define GRAPHLAB_SERIALIZE_LIST_HPP

#include <list>

#include <serialization/iarchive.hpp>
#include <serialization/oarchive.hpp>
#include <serialization/iterator.hpp>


namespace graphlab {
namespace archive_detail {
  /** serializes a list  */
  template <typename ArcType, typename T>
  struct serialize_impl<ArcType, std::list<T> > {
  static void exec(ArcType& a, const std::list<T>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a list  */
  template <typename ArcType, typename T>
  struct deserialize_impl<ArcType, std::list<T> > {
  static void exec(ArcType& a, std::list<T>& vec){
    vec.clear();
    deserialize_iterator<T>(a, std::inserter(vec,vec.end()));
  }
  };
} // archive_detail  
} // graphlab
#endif 
