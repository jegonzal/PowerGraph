#ifndef SERIALIZABLE_POD_HPP
#define SERIALIZABLE_POD_HPP
#include <graphlab/serialization/vector.hpp>

// if serializable POD is set on a struct / class, we generate a fast serializer / deserializer
// that will lose some cross platform support
// in particular, all platforms will need sizeof(struct) to be identical

#define SERIALIZABLE_POD(tname) \
  namespace graphlab { namespace archive_detail {             \
  template <typename ArcType> struct serialize_impl<ArcType, tname>{                        \
    static void exec(ArcType &a, const tname &i) {              \
    a.o->write(reinterpret_cast<const char*>(&i), sizeof(tname));\
    }                                                         \
  };                                                          \
  template <typename ArcType> struct deserialize_impl<ArcType, tname>{                      \
    static void exec(ArcType &a, tname &t) {                    \
    a.i->read(reinterpret_cast<char*>(&t), sizeof(tname)); \
    }                                                         \
  };                                                          \
  __GEN_FAST_VECTOR_SERIALIZE__(tname)                       \
  }}

#endif

