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


namespace archive_detail {

/**
Implementation of the serializer for different types.
This is the catch-all and is used to call the .save function
if T is a class. Fails at runtime otherwise.
*/
template <typename ArcType, typename T>
struct serialize_impl {
  static void exec(ArcType &o, const T& t) {
    save_or_fail(o, t);
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
