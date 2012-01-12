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
   * The input archive object.
   * It is just a simple wrapper around a istream 
   */
  class iarchive {
  public:
    
    bool directbuffer; // if false, uses the istream, otherwise uses buf
    
    char* buf; size_t buflen; size_t bufread; bool buffail;
    
    std::istream* i;
    
    inline char read_char() {
      if (directbuffer == false) {
        char c;
        i->get(c);
        return c;
      }
      else {
        if (bufread < buflen) {
          ++bufread;
          return *buf++;
        }
        else {
          buffail = true;
          return 0;
        }
      }
    }
    
    inline bool has_directbuffer() {
      return directbuffer;
    }
    
    inline char* get_direct_buffer(size_t len) {
      if (len + bufread <= buflen) {
        char* ret = buf;
        buf += len;
        bufread += len;
        return ret;
      }
      else {
        buffail = true;
        return NULL;
      }
    }
    
    inline void read(char* c, size_t len) {
      if (directbuffer == false) i->read(c, len);
      else {
        if (len + bufread <= buflen) {
          memcpy(c, buf, len);
          buf += len;
          bufread += len;
        }
        else {
          buffail = true;
        }
      }
    }
    
    inline bool fail() {
      return (directbuffer == false)?i->fail():buffail;
    }
    /// constructor. Takes a generic std::istream object
    inline iarchive(std::istream& is)
      :directbuffer(false),i(&is) { }

    /// constructor. Takes a generic std::istream object
    inline iarchive(char* buf, size_t buflen)
          :directbuffer(true),
           buf(buf), buflen(buflen), bufread(0), buffail(false),
           i(NULL) { }

    ~iarchive() {}
  };


  /**
   * An alternate input archive object. 
   * When this object is used to deserialize an object,
   * and the object does not support serialization,
   * failure will only occur at runtime.
   * \see iarchive
   */
  class iarchive_soft_fail{
  public:
    bool directbuffer; // if false, uses the istream, otherwise uses buf
    
    char* buf; size_t buflen; size_t bufread; bool buffail;
    
    std::istream* i;
    
    inline char read_char() {
      if (directbuffer == false) {
        char c;
        i->get(c);
        return c;
      }
      else {
        if (bufread < buflen) {
          ++bufread;
          return *buf++;
        }
        else {
          buffail = true;
          return 0;
        }
      }
    }
    
    inline bool has_directbuffer() {
      return directbuffer;
    }
    
    inline char* get_direct_buffer(size_t len) {
      if (len + bufread <= buflen) {
        char* ret = buf;
        buf += len;
        bufread += len;
        return ret;
      }
      else {
        buffail = true;
        return NULL;
      }
    }
    
    inline void read(char* c, size_t len) {
      if (directbuffer == false) i->read(c, len);
      else {
        if (len + bufread <= buflen) {
          memcpy(c, buf, len);
          buf += len;
          bufread += len;
        }
        else {
          buffail = true;
        }
      }
    }
    
    inline bool fail() {
      return (directbuffer == false)?i->fail():buffail;
    }
    
    inline iarchive_soft_fail(std::istream &is)
      : i(&is) {}

    /// constructor. Takes a generic std::istream object
    inline iarchive_soft_fail(char* buf, size_t buflen)
          :directbuffer(true),
           buf(buf), buflen(buflen), bufread(0), buffail(false),
           i(NULL) { }


    inline iarchive_soft_fail(iarchive &iarc)
      :directbuffer(iarc.directbuffer), buf(iarc.buf), buflen(iarc.buflen),
       bufread(iarc.bufread), buffail(iarc.buffail), i(iarc.i){}
  
    ~iarchive_soft_fail() { }
  };


  namespace archive_detail {

    /// called by the regular archive The regular archive will do a hard fail
    template <typename ArcType, typename T>
    struct deserialize_hard_or_soft_fail {
      inline static void exec(ArcType &i, T& t) {
        t.load(i);
      }
    };

    /// called by the soft fail archive 
    template <typename T>
    struct deserialize_hard_or_soft_fail<iarchive_soft_fail, T> {
      inline static void exec(iarchive_soft_fail &i, T& t) {
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
      inline static void exec(ArcType &i, T& t) {
        deserialize_hard_or_soft_fail<ArcType, T>::exec(i, t);
      }
    };

    // catch if type is a POD
    template <typename ArcType, typename T>
    struct deserialize_impl<ArcType, T, true>{
      inline static void exec(ArcType &a, T &t) {
        a.read(reinterpret_cast<char*>(&t), sizeof(T));
      }
    };

  } //namespace archive_detail

  /**
     Allows Use of the "stream" syntax for serialization 
  */
  template <typename T>
  inline iarchive& operator>>(iarchive& a, T &i) {
    archive_detail::deserialize_impl<iarchive, T, gl_is_pod<T>::value >::exec(a, i);
    return a;
  }



  /**
     Allows Use of the "stream" syntax for serialization 
  */
  template <typename T>
  inline iarchive_soft_fail& operator>>(iarchive_soft_fail& a, T &i) {
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
    a.read(reinterpret_cast<char*>(i), (std::streamsize)length);
    assert(!a.fail());
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
    a.read(reinterpret_cast<char*>(i), (std::streamsize)length);
    assert(!a.fail());
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
#define BEGIN_OUT_OF_PLACE_LOAD(arc, tname, tval)                       \
  namespace graphlab{ namespace archive_detail {                        \
  template <typename ArcType> struct deserialize_impl<ArcType, tname, false>{ \
  static void exec(ArcType& arc, tname & tval) {             

#define END_OUT_OF_PLACE_LOAD() } }; } }




} // namespace graphlab


#endif 

#endif








