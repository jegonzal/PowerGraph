/**
  This files defines the serializer/deserializer for all basic types
  (as well as string and pair)  
*/
#ifndef ARCHIVE_BASIC_TYPES_HPP
#define ARCHIVE_BASIC_TYPES_HPP

#include <string>
#include <graphlab/serialization/integer.hpp>
#include <graphlab/serialization/serializable_pod.hpp>
#include <graphlab/logger/assertions.hpp>

/***************************************************************************
 *                        Basic Serializers                                * 
 ***************************************************************************/
// generate the operator<< call for a whole bunch of integer types


#define INT_SERIALIZE(tname)                                            \
  template <typename ArcType> struct serialize_impl<ArcType, tname>{    \
    static void exec(ArcType &a, const tname &i_) {                     \
      int64_t i = i_ ;                                                  \
      char c[10];                                                       \
      unsigned char len = compress_int(i, c);                           \
      a.o->write(c + 10 - len, len);                                    \
    }                                                                   \
  };                                                                    \
  template <typename ArcType> struct deserialize_impl<ArcType, tname>{  \
    static void exec(ArcType &a, tname &t_) {                           \
      decompress_int<tname>(*(a.i), t_);                                \
    }                                                                   \
  };




namespace graphlab {
class oarchive;
class iarchive;
}


SERIALIZABLE_POD(char);
SERIALIZABLE_POD(bool);
SERIALIZABLE_POD(unsigned char);
SERIALIZABLE_POD(double);
SERIALIZABLE_POD(float);


namespace graphlab {
namespace archive_detail {

INT_SERIALIZE(short);
INT_SERIALIZE(unsigned short);
INT_SERIALIZE(int);
INT_SERIALIZE(long);
INT_SERIALIZE(long long);
INT_SERIALIZE(unsigned long);
INT_SERIALIZE(unsigned int);
INT_SERIALIZE(unsigned long long);


/********Serialization and deserialiation of char* **************/
template <typename ArcType>
struct serialize_impl<ArcType, const char*> {
  static void exec(ArcType &a, const char* const &s) {
    // save the length
    // ++ for the \0
    size_t length = strlen(s); length++;
    serialize_impl<ArcType, size_t>::exec(a, length);
    a.o->write(reinterpret_cast<const char*>(s), length);
    DASSERT_FALSE(a.o->fail());
  }
};


template <typename ArcType, size_t len>
struct serialize_impl<ArcType, char [len]> {
  static void exec(ArcType& a, const char s[len] ) { 
    size_t length = len;
    serialize_impl<ArcType, size_t>::exec(a, length);
    a.o->write(reinterpret_cast<const char*>(s), length);
    DASSERT_FALSE(a.o->fail());
  }
};

template <typename ArcType>
struct serialize_impl<ArcType, char*> {
  static void exec(ArcType &a, char* const &s) {
    // save the length
    // ++ for the \0
    size_t length = strlen(s); length++;
    serialize_impl<ArcType, size_t>::exec(a, length);
    a.o->write(reinterpret_cast<const char*>(s), length);
    DASSERT_FALSE(a.o->fail());
  }
};

template <typename ArcType>
struct deserialize_impl<ArcType, char*> {
  static void exec(ArcType& a, char*& s) {
    // Save the length and check if lengths match
    size_t length;
    deserialize_impl<ArcType, size_t>::exec(a, length);
    s = new char[length];
    //operator>> the rest
    a.i->read(reinterpret_cast<char*>(s), length);
    DASSERT_FALSE(a.i->fail());
  }
};
  
template <typename ArcType, size_t len>
struct deserialize_impl<ArcType, char [len]> {
  static void exec(ArcType& a, char s[len]) { 
    size_t length;
    deserialize_impl<ArcType, size_t>::exec(a, length);
    ASSERT_LE(length, len);
    a.i->read(reinterpret_cast<char*>(s), length);
    DASSERT_FALSE(a.i->fail());
  }
};



/********Serialization and deserialiation of strings **************/
template <typename ArcType>
struct serialize_impl<ArcType, std::string> {
  static void exec(ArcType &a, const std::string& s) {
    size_t length = s.length();
    serialize_impl<ArcType, size_t>::exec(a, length);
    a.o->write(reinterpret_cast<const char*>(s.c_str()), length);
    DASSERT_FALSE(a.o->fail());
  }
};


template <typename ArcType>
struct deserialize_impl<ArcType, std::string> {
  static void exec(ArcType &a, std::string &s) {
      //read the length
    size_t length;
    deserialize_impl<ArcType, size_t>::exec(a, length);
    //resize the string and read the characters
    s.resize(length);
    a.i->read(const_cast<char*>(s.c_str()), length);
    DASSERT_FALSE(a.i->fail());
  }
};

/******** Serialization and deserialization of pairs *************/


/********Serialization and deserialiation of strings **************/
template <typename ArcType, typename T, typename U>
struct serialize_impl<ArcType, std::pair<T, U> > {
  static void exec(ArcType &a, const std::pair<T, U> &s) {
    serialize_impl<ArcType, T>::exec(a, s.first);
    serialize_impl<ArcType, U>::exec(a, s.second);
  }
};



template <typename ArcType, typename T, typename U>
struct deserialize_impl<ArcType, std::pair<T, U> > {
  static void exec(ArcType &a, std::pair<T, U> &s) {
    deserialize_impl<ArcType, T>::exec(a, s.first);
    deserialize_impl<ArcType, U>::exec(a, s.second);
  }
};

} // namespace archive_detail
} // namespace graphlab
 
#undef INT_SERIALIZE
#endif
