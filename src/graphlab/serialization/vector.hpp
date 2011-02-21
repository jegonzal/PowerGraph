#ifndef GRAPHLAB_SERIALIZE_VECTOR_HPP
#define GRAPHLAB_SERIALIZE_VECTOR_HPP
#include <vector>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>


namespace graphlab {
  namespace archive_detail {
    /**
       Serializes a vector
       Returns true on success, false on failure  */
    template <typename ArcType, typename ValueType>
    struct serialize_impl<ArcType, std::vector<ValueType> > {
      static void exec(ArcType &a, const std::vector<ValueType>& vec) {
        serialize_impl<ArcType, size_t>::exec(a, vec.size());
        serialize_iterator(a,vec.begin(), vec.end());
      }
    };
    /**
       deserializes a vector
       Returns true on success, false on failure  */
    template <typename ArcType, typename ValueType>
    struct deserialize_impl<ArcType, std::vector<ValueType> > {
      static void exec(ArcType& a, std::vector<ValueType>& vec){
        size_t len;
        deserialize_impl<ArcType, size_t>::exec(a, len);
        vec.clear(); vec.reserve(len);
        deserialize_iterator<ArcType, ValueType>(a, std::inserter(vec, vec.end()));
      }
    };
                                               
    /** fast serializer for other scalar types
        can't do it for integers as our integer serializer pushes it to 64-bit 
        before writing to get a limited amount of cross-platform support
        This should not be used directly. */ 
#define __GEN_FAST_VECTOR_SERIALIZE__(tname)                            \
    template <typename ArcType>                                         \
    struct serialize_impl<ArcType, std::vector<tname> > {               \
      static void exec(ArcType &a, const std::vector<tname>& vec) {     \
        serialize_impl<ArcType, size_t>::exec(a, vec.size());           \
        serialize(a, &(vec[0]),sizeof(tname)*vec.size());               \
      }                                                                 \
    };                                                                  \
    template <typename ArcType>                                         \
    struct deserialize_impl<ArcType, std::vector<tname> > {             \
      static void exec(ArcType &a, std::vector<tname>& vec) {           \
        size_t len;                                                     \
        deserialize_impl<ArcType, size_t>::exec(a, len);                \
        vec.clear(); vec.resize(len);                                   \
        deserialize(a, &(vec[0]), sizeof(tname)*vec.size());            \
      }                                                                 \
    };
 
  } // archive_detail
} // namespace graphlab

#endif 
