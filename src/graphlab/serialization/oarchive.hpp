#ifndef GRAPHLAB_OARCHIVE_HPP
#define GRAPHLAB_OARCHIVE_HPP

#include <iostream>
#include <string>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/assert.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/has_save.hpp>
namespace graphlab {


/**
The output archive object.
It is just a simple wrapper around a ostream
*/
class oarchive{
 public:
  std::ostream* o;

  oarchive(std::ostream& os)
    : o(&os) {}

  ~oarchive() { }
};


class oarchive_soft_fail{
 public:
  std::ostream* o;

  oarchive_soft_fail(std::ostream& os)
    : o(&os) {}

  oarchive_soft_fail(oarchive &oarc):o(oarc.o) {}
  
  ~oarchive_soft_fail() { }
};

namespace archive_detail {


template <typename ArcType, typename T>
struct serialize_hard_or_soft_fail {
  static void exec(ArcType &o, const T& t) {
    t.save(o);
  }
};


template <typename T>
struct serialize_hard_or_soft_fail<oarchive_soft_fail, T> {
  static void exec(oarchive_soft_fail &o, const T& t) {
    oarchive oarc(*(o.o));
    save_or_fail(oarc, t);
  }
};


/**
Implementation of the serializer for different types.
This is the catch-all and is used to call the .save function
if T is a class. Fails at runtime otherwise.
*/
template <typename ArcType, typename T>
struct serialize_impl {
  static void exec(ArcType &o, const T& t) {
    serialize_hard_or_soft_fail<ArcType, T>::exec(o, t);
  }
};


/**
Re-dispatch if for some reasons T already has a const
*/
template <typename ArcType, typename T>
struct serialize_impl<ArcType, const T> {
  static void exec(ArcType &o, const T& t) {
    serialize_impl<ArcType, T>::exec(o, t);
  }
};
}// archive_detail


/**
Allows Use of the "stream" syntax for serialization 
*/
template <typename T>
oarchive& operator<<(oarchive& a, const T& i) {
  archive_detail::serialize_impl<oarchive, T>::exec(a, i);
  return a;
}

template <typename T>
oarchive_soft_fail& operator<<(oarchive_soft_fail& a, const T& i) {
  archive_detail::serialize_impl<oarchive_soft_fail, T>::exec(a, i);
  return a;
}


/**
Serializes an arbitrary pointer + length to an archive 
*/
inline oarchive& serialize(oarchive& a, const void* i,const size_t length) {
  // save the length
  operator<<(a,length);
  a.o->write(reinterpret_cast<const char*>(i), length);
  assert(!a.o->fail());
  return a;
}


/**
Serializes an arbitrary pointer + length to an archive 
*/
inline oarchive_soft_fail& serialize(oarchive_soft_fail& a, const void* i,const size_t length) {
  // save the length
  operator<<(a,length);
  a.o->write(reinterpret_cast<const char*>(i), length);
  assert(!a.o->fail());
  return a;
}


}

/**
Macro to make it easy to define out-of-place saves (and loads)
to define an "out of place" save
OUT_OF_PLACE_SAVE(arc, typename, tval) 
  arc << tval;    // do whatever serialization stuff you need here
END_OUT_OF_PLACE_SAVE()

\note important! this must be defined in the global namespace!
See unsupported_serialize for an example
*/
#define BEGIN_OUT_OF_PLACE_SAVE(arc, tname, tval) \
  namespace graphlab{ namespace archive_detail {  \
  template <typename ArcType> struct serialize_impl<ArcType, tname> {     \
    static void exec(ArcType& arc, const tname & tval) {

#define END_OUT_OF_PLACE_SAVE() } }; } }

#include <graphlab/serialization/basic_types.hpp>
#endif  
