#ifndef GRAPHLAB_IARCHIVE_HPP
#define GRAPHLAB_IARCHIVE_HPP

#include <iostream>
#include <logger/assertions.hpp>
#include <serialization/has_load.hpp>
namespace graphlab {
/**
The input archive object.
It is just a simple wrapper around a istream
*/
class iarchive {
 public:
  std::istream* i;

  iarchive(std::istream& is)
    :i(&is) { }
};

namespace archive_detail {
/**
Implementation of the serializer for different types.
This is the catch-all and is used to call the .load function
if T is a class. Fails at runtime otherwise.
*/
template <typename ArcType, typename T>
struct deserialize_impl {
  static void exec(ArcType &i, T& t) {
    load_or_fail(i, t);
  }
};

} //namespace archive_detail

/**
Allows Use of the "stream" syntax for serialization 
*/
template <typename T>
iarchive& operator>>(iarchive& a, T &i) {
  archive_detail::deserialize_impl<iarchive, T>::exec(a, i);
  return a;
}

/**
deserializes an arbitrary pointer + length from an archive 
*/
inline iarchive& deserialize(iarchive& a, void* const i,const size_t length) {
  // Save the length and check if lengths match
  size_t length2;
  operator>>(a, length2);
  ASSERT_EQ(length, length2);

  //operator>> the rest
  a.i->read(reinterpret_cast<char*>(i), length);
  assert(!a.i->fail());
  return a;
}


} // namespace graphlab

#include <serialization/basic_types.hpp>
#endif 
