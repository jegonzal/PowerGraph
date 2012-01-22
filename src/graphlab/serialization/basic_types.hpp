/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
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
      uint64_t i = (uint64_t)(i_) ;                                     \
      char c[10];                                                       \
      unsigned char len = compress_int(i, c);                           \
      a.o->write(c + 10 - len, (std::streamsize)len);                   \
    }                                                                   \
  };                                                                    \
  template <typename ArcType> struct deserialize_impl<ArcType, tname, false>{ \
    static void exec(ArcType &a, tname &t_) {                           \
      decompress_int<ArcType, tname>(a, t_);                                \
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
        a.read(reinterpret_cast<char*>(s), length);
        DASSERT_FALSE(a.i->fail());
      }
    };
  
    template <typename ArcType, size_t len>
    struct deserialize_impl<ArcType, char [len], false> {
      static void exec(ArcType& a, char s[len]) { 
        size_t length;
        deserialize_impl<ArcType, size_t, false>::exec(a, length);
        ASSERT_LE(length, len);
        a.read(reinterpret_cast<char*>(s), length);
        DASSERT_FALSE(a.i->fail());
      }
    };



    /********Serialization and deserialiation of strings **************/
    template <typename ArcType>
    struct serialize_impl<ArcType, std::string, false> {
      static void exec(ArcType &a, const std::string& s) {
        size_t length = s.length();
        serialize_impl<ArcType, size_t, false>::exec(a, length);
        a.o->write(reinterpret_cast<const char*>(s.c_str()), (std::streamsize)length);
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
        a.read(const_cast<char*>(s.c_str()), (std::streamsize)length);
        DASSERT_FALSE(a.fail());
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

