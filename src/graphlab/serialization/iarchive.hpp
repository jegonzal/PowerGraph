#ifndef GRAPHLAB_SERIALIZE_HPP
#include <graphlab/serialization/serialize.hpp>

#else


#ifndef GRAPHLAB_IARCHIVE_HPP
#define GRAPHLAB_IARCHIVE_HPP

#include <iostream>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/is_pod.hpp>
#include <graphlab/serialization/has_load.hpp>
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
      
  ~iarchive() {}
};



class iarchive_soft_fail{
 public:
  std::istream *i;

  iarchive_soft_fail(std::istream &is)
    : i(&is) {}

  iarchive_soft_fail(iarchive &iarc):i(iarc.i){}
  
  ~iarchive_soft_fail() { }
};


namespace archive_detail {

/// called by the regular archive The regular archive will do a hard fail
template <typename ArcType, typename T>
struct deserialize_hard_or_soft_fail {
  static void exec(ArcType &i, T& t) {
    t.load(i);
  }
};

/// called by the soft fail archive 
template <typename T>
struct deserialize_hard_or_soft_fail<iarchive_soft_fail, T> {
  static void exec(iarchive_soft_fail &i, T& t) {
    iarchive iarc(*(i.i));
    load_or_fail(iarc, t);
  }
};


/**
Implementation of the deserializer for different types.
This is the catch-all. If it gets here, it must be a non-POD and is a class.
We therefore call the .save function.
Here we pick between the archive types using serialize_hard_or_soft_fail
*/
template <typename ArcType, typename T, bool IsPOD>
struct deserialize_impl {
  static void exec(ArcType &i, T& t) {
    deserialize_hard_or_soft_fail<ArcType, T>::exec(i, t);
  }
};

// catch if type is a POD
template <typename ArcType, typename T>
struct deserialize_impl<ArcType, T, true>{
  static void exec(ArcType &a, T &t) {
    a.i->read(reinterpret_cast<char*>(&t), sizeof(T));
  }
};

} //namespace archive_detail

/**
Allows Use of the "stream" syntax for serialization 
*/
template <typename T>
iarchive& operator>>(iarchive& a, T &i) {
  archive_detail::deserialize_impl<iarchive, T, gl_is_pod<T>::value >::exec(a, i);
  return a;
}



/**
Allows Use of the "stream" syntax for serialization 
*/
template <typename T>
iarchive_soft_fail& operator>>(iarchive_soft_fail& a, T &i) {
  archive_detail::deserialize_impl<iarchive_soft_fail, T, gl_is_pod<T>::value >::exec(a, i);
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



/**
deserializes an arbitrary pointer + length from an archive 
*/
inline iarchive_soft_fail& deserialize(iarchive_soft_fail& a, void* const i,const size_t length) {
  // Save the length and check if lengths match
  size_t length2;
  operator>>(a, length2);
  ASSERT_EQ(length, length2);

  //operator>> the rest
  a.i->read(reinterpret_cast<char*>(i), length);
  assert(!a.i->fail());
  return a;
}

/**
Macro to make it easy to define out-of-place saves (and loads)
to define an "out of place" load
OUT_OF_PLACE_LOAD(arc, typename, tval) 
  arc >> tval;    // do whatever deserialization stuff you need here
END_OUT_OF_PLACE_LOAD()

\note important! this must be defined in the global namespace!
See unsupported_serialize for an example
*/
#define BEGIN_OUT_OF_PLACE_LOAD(arc, tname, tval) \
  namespace graphlab{ namespace archive_detail {  \
  template <typename ArcType> struct deserialize_impl<ArcType, tname, false>{     \
    static void exec(ArcType& arc, tname & tval) {             

#define END_OUT_OF_PLACE_LOAD() } }; } }




} // namespace graphlab


#endif 

#endif