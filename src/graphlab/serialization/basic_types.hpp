/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

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
  template <typename ArcType> struct serialize_impl<ArcType, tname, false>{ \
    static void exec(ArcType &a, const tname &i_) {                     \
      int64_t i = i_ ;                                                  \
      char c[10];                                                       \
      unsigned char len = compress_int(i, c);                           \
      a.o->write(c + 10 - len, len);                                    \
      /* a.o->write(c, len);     */                                     \
    }                                                                   \
  };                                                                    \
  template <typename ArcType> struct deserialize_impl<ArcType, tname, false>{ \
    static void exec(ArcType &a, tname &t_) {                           \
      decompress_int<tname>(*(a.i), t_);                                \
    }                                                                   \
  };




namespace graphlab {
  class oarchive;
  class iarchive;
}


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
    struct serialize_impl<ArcType, const char*, false> {
      static void exec(ArcType &a, const char* const &s) {
        // save the length
        // ++ for the \0
        size_t length = strlen(s); length++;
        serialize_impl<ArcType, size_t, false>::exec(a, length);
        a.o->write(reinterpret_cast<const char*>(s), length);
        DASSERT_FALSE(a.o->fail());
      }
    };


    template <typename ArcType, size_t len>
    struct serialize_impl<ArcType, char [len], false> {
      static void exec(ArcType& a, const char s[len] ) { 
        size_t length = len;
        serialize_impl<ArcType, size_t, false>::exec(a, length);
        a.o->write(reinterpret_cast<const char*>(s), length);
        DASSERT_FALSE(a.o->fail());
      }
    };

    template <typename ArcType>
    struct serialize_impl<ArcType, char*, false> {
      static void exec(ArcType &a, char* const &s) {
        // save the length
        // ++ for the \0
        size_t length = strlen(s); length++;
        serialize_impl<ArcType, size_t, false>::exec(a, length);
        a.o->write(reinterpret_cast<const char*>(s), length);
        DASSERT_FALSE(a.o->fail());
      }
    };

    template <typename ArcType>
    struct deserialize_impl<ArcType, char*, false> {
      static void exec(ArcType& a, char*& s) {
        // Save the length and check if lengths match
        size_t length;
        deserialize_impl<ArcType, size_t, false>::exec(a, length);
        s = new char[length];
        //operator>> the rest
        a.i->read(reinterpret_cast<char*>(s), length);
        DASSERT_FALSE(a.i->fail());
      }
    };
  
    template <typename ArcType, size_t len>
    struct deserialize_impl<ArcType, char [len], false> {
      static void exec(ArcType& a, char s[len]) { 
        size_t length;
        deserialize_impl<ArcType, size_t, false>::exec(a, length);
        ASSERT_LE(length, len);
        a.i->read(reinterpret_cast<char*>(s), length);
        DASSERT_FALSE(a.i->fail());
      }
    };



    /********Serialization and deserialiation of strings **************/
    template <typename ArcType>
    struct serialize_impl<ArcType, std::string, false> {
      static void exec(ArcType &a, const std::string& s) {
        size_t length = s.length();
        serialize_impl<ArcType, size_t, false>::exec(a, length);
        a.o->write(reinterpret_cast<const char*>(s.c_str()), length);
        DASSERT_FALSE(a.o->fail());
      }
    };


    template <typename ArcType>
    struct deserialize_impl<ArcType, std::string, false> {
      static void exec(ArcType &a, std::string &s) {
        //read the length
        size_t length;
        deserialize_impl<ArcType, size_t, false>::exec(a, length);
        //resize the string and read the characters
        s.resize(length);
        a.i->read(const_cast<char*>(s.c_str()), length);
        DASSERT_FALSE(a.i->fail());
      }
    };

    /******** Serialization and deserialization of pairs *************/


    /********Serialization and deserialiation of strings **************/
    template <typename ArcType, typename T, typename U>
    struct serialize_impl<ArcType, std::pair<T, U>, false > {
      static void exec(ArcType &a, const std::pair<T, U> &s) {
        serialize_impl<ArcType, T, gl_is_pod<T>::value >::exec(a, s.first);
        serialize_impl<ArcType, U, gl_is_pod<U>::value >::exec(a, s.second);
      }
    };



    template <typename ArcType, typename T, typename U>
    struct deserialize_impl<ArcType, std::pair<T, U>, false > {
      static void exec(ArcType &a, std::pair<T, U> &s) {
        deserialize_impl<ArcType, T, gl_is_pod<T>::value >::exec(a, s.first);
        deserialize_impl<ArcType, U, gl_is_pod<U>::value >::exec(a, s.second);
      }
    };

  } // namespace archive_detail
} // namespace graphlab
 
#undef INT_SERIALIZE
#endif
