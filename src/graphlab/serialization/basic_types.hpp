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


/*
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
  template <typename OutArcType>                                           \
  struct serialize_impl<OutArcType, tname, false>{                         \
    static void exec(OutArcType& oarc, const tname& t) {                   \
      oarc.write(reinterpret_cast<const char*>(&t), sizeof(tname));     \
    }                                                                   \
  };                                                                    \
  template <typename InArcType>                                           \
  struct deserialize_impl<InArcType, tname, false>{                       \
    static void exec(InArcType& iarc, tname& t) {                         \
      decompress_int<InArcType, tname>(iarc, t);                          \
    }                                                                   \
  };




namespace graphlab {
  class oarchive;
  class iarchive;
}


namespace graphlab {
  namespace archive_detail {
    /*
     * Generate serializers and deserializers for all integer types
     */
    INT_SERIALIZE(bool);
    INT_SERIALIZE(char);
    INT_SERIALIZE(unsigned char);
    INT_SERIALIZE(short);
    INT_SERIALIZE(unsigned short);
    INT_SERIALIZE(int);
    INT_SERIALIZE(long);
    INT_SERIALIZE(long long);
    INT_SERIALIZE(unsigned long);
    INT_SERIALIZE(unsigned int);
    INT_SERIALIZE(unsigned long long);


    /** Serialization of null terminated const char* strings.
     * This is necessary to serialize constant strings like
     * \code 
     * oarc << "hello world";
     * \endcode
     */
    template <typename OutArcType>
    struct serialize_impl<OutArcType, const char*, false> {
      static void exec(OutArcType& oarc, const char* const& s) {
        // save the length
        // ++ for the \0
        size_t length = strlen(s); length++;
        serialize_impl<OutArcType, size_t, false>::exec(oarc, length);
        oarc.write(reinterpret_cast<const char*>(s), length);
        DASSERT_FALSE(oarc.fail());
      }
    };


    /// Serialization of fixed length char arrays
    template <typename OutArcType, size_t len>
    struct serialize_impl<OutArcType, char [len], false> {
      static void exec(OutArcType& oarc, const char s[len] ) { 
        size_t length = len;
        serialize_impl<OutArcType, size_t, false>::exec(oarc, length);
        oarc.write(reinterpret_cast<const char*>(s), length);
        DASSERT_FALSE(oarc.fail());
      }
    };


    /// Serialization of null terminated char* strings
    template <typename OutArcType>
    struct serialize_impl<OutArcType, char*, false> {
      static void exec(OutArcType& oarc, char* const& s) {
        // save the length
        // ++ for the \0
        size_t length = strlen(s); length++;
        serialize_impl<OutArcType, size_t, false>::exec(oarc, length);
        oarc.write(reinterpret_cast<const char*>(s), length);
        DASSERT_FALSE(oarc.fail());
      }
    };

    /// Deserialization of null terminated char* strings
    template <typename InArcType>
    struct deserialize_impl<InArcType, char*, false> {
      static void exec(InArcType& iarc, char*& s) {
        // Save the length and check if lengths match
        size_t length;
        deserialize_impl<InArcType, size_t, false>::exec(iarc, length);
        s = new char[length];
        //operator>> the rest
        iarc.read(reinterpret_cast<char*>(s), length);
        DASSERT_FALSE(iarc.fail());
      }
    };
  
    /// Deserialization of fixed length char arrays 
    template <typename InArcType, size_t len>
    struct deserialize_impl<InArcType, char [len], false> {
      static void exec(InArcType& iarc, char s[len]) { 
        size_t length;
        deserialize_impl<InArcType, size_t, false>::exec(iarc, length);
        ASSERT_LE(length, len);
        iarc.read(reinterpret_cast<char*>(s), length);
        DASSERT_FALSE(iarc.fail());
      }
    };



    /// Serialization of std::string
    template <typename OutArcType>
    struct serialize_impl<OutArcType, std::string, false> {
      static void exec(OutArcType& oarc, const std::string& s) {
        size_t length = s.length();
        serialize_impl<OutArcType, size_t, false>::exec(oarc, length);
        oarc.write(reinterpret_cast<const char*>(s.c_str()), 
                   (std::streamsize)length);
        DASSERT_FALSE(oarc.fail());
      }
    };


    /// Deserialization of std::string
    template <typename InArcType>
    struct deserialize_impl<InArcType, std::string, false> {
      static void exec(InArcType& iarc, std::string& s) {
        //read the length
        size_t length;
        deserialize_impl<InArcType, size_t, false>::exec(iarc, length);
        //resize the string and read the characters
        s.resize(length);
        iarc.read(const_cast<char*>(s.c_str()), (std::streamsize)length);
        DASSERT_FALSE(iarc.fail());
      }
    };

    /// Serialization of std::pair
    template <typename OutArcType, typename T, typename U>
    struct serialize_impl<OutArcType, std::pair<T, U>, false > {
      static void exec(OutArcType& oarc, const std::pair<T, U>& s) {
        serialize_impl<OutArcType, 
                       T, gl_is_pod<T>::value >::exec(oarc, s.first);
        serialize_impl<OutArcType, 
                       U, gl_is_pod<U>::value >::exec(oarc, s.second);
      }
    };


    /// Deserialization of std::pair
    template <typename InArcType, typename T, typename U>
    struct deserialize_impl<InArcType, std::pair<T, U>, false > {
      static void exec(InArcType& iarc, std::pair<T, U>& s) {
        deserialize_impl<InArcType, 
                         T, gl_is_pod<T>::value >::exec(iarc, s.first);
        deserialize_impl<InArcType, 
                         U, gl_is_pod<U>::value >::exec(iarc, s.second);
      }
    };

  } // namespace archive_detail
} // namespace graphlab
 
#undef INT_SERIALIZE
#endif

