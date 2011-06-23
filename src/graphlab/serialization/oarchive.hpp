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


// This file should not be included directly. use serialize.hpp
#ifndef GRAPHLAB_SERIALIZE_HPP
#include <graphlab/serialization/serialize.hpp>

#else

#ifndef GRAPHLAB_OARCHIVE_HPP
#define GRAPHLAB_OARCHIVE_HPP

#include <iostream>
#include <string>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/is_pod.hpp>
#include <graphlab/serialization/has_save.hpp>

namespace graphlab {


  /**
   *  The output archive object.
   *  It is just a simple wrapper around a ostream. 
   */
  class oarchive{
  public:
    std::ostream* o;
    /// constructor. Takes a generic std::ostream object
    oarchive(std::ostream& os)
      : o(&os) {}

    ~oarchive() { }
  };

/**
   * An alternate output archive object. 
   * When this object is used to serialize an object,
   * and the object does not support serialization,
   * failure will only occur at runtime.
   * \see oarchive
   */
  class oarchive_soft_fail{
  public:
    std::ostream* o;

    oarchive_soft_fail(std::ostream& os)
      : o(&os) {}

    oarchive_soft_fail(oarchive &oarc):o(oarc.o) {}
  
    ~oarchive_soft_fail() { }
  };

  namespace archive_detail {

    /// called by the regular archive The regular archive will do a hard fail
    template <typename ArcType, typename T>
    struct serialize_hard_or_soft_fail {
      static void exec(ArcType &o, const T& t) {
        t.save(o);
      }
    };

    /// called by the soft fail archive 
    template <typename T>
    struct serialize_hard_or_soft_fail<oarchive_soft_fail, T> {
      static void exec(oarchive_soft_fail &o, const T& t) {
        oarchive oarc(*(o.o));
        save_or_fail(oarc, t);
      }
    };


    /**
       Implementation of the serializer for different types.
       This is the catch-all. If it gets here, it must be a non-POD and is a class.
       We therefore call the .save function.
       Here we pick between the archive types using serialize_hard_or_soft_fail
    */
    template <typename ArcType, typename T, bool IsPOD>
    struct serialize_impl {
      static void exec(ArcType &o, const T& t) {
        serialize_hard_or_soft_fail<ArcType, T>::exec(o, t);
      }
    };

    /** Catch if type is a POD */
    template <typename ArcType, typename T>
    struct serialize_impl<ArcType, T, true> {
      static void exec(ArcType &a, const T& t) {
        a.o->write(reinterpret_cast<const char*>(&t), sizeof(T));
      }
    };

    /**
       Re-dispatch if for some reasons T already has a const
    */
    template <typename ArcType, typename T, bool IsPOD>
    struct serialize_impl<ArcType, const T, IsPOD> {
      static void exec(ArcType &o, const T& t) {
        serialize_impl<ArcType, T, IsPOD>::exec(o, t);
      }
    };
  }// archive_detail


  /**
     Allows Use of the "stream" syntax for serialization 
  */
  template <typename T>
  oarchive& operator<<(oarchive& a, const T& i) {
    archive_detail::serialize_impl<oarchive, T, gl_is_pod<T>::value >::exec(a, i);
    return a;
  }

  template <typename T>
  oarchive_soft_fail& operator<<(oarchive_soft_fail& a, const T& i) {
    archive_detail::serialize_impl<oarchive_soft_fail, T, gl_is_pod<T>::value >::exec(a, i);
    return a;
  }


  /**
     Serializes an arbitrary pointer + length to an archive 
  */
  inline oarchive& serialize(oarchive& a, const void* i,const size_t length) {
    // save the length
    operator<<(a,length);
    a.o->write(reinterpret_cast<const char*>(i), (std::streamsize)length);
    assert(!a.o->fail());
    return a;
  }


  /**
     Serializes an arbitrary pointer + length to an archive 
  */
  inline oarchive_soft_fail& serialize(oarchive_soft_fail& a, const void* i,const size_t length) {
    // save the length
    operator<<(a,length);
    a.o->write(reinterpret_cast<const char*>(i), (std::streamsize)length);
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
#define BEGIN_OUT_OF_PLACE_SAVE(arc, tname, tval)                       \
  namespace graphlab{ namespace archive_detail {                        \
  template <typename ArcType> struct serialize_impl<ArcType, tname, false> { \
  static void exec(ArcType& arc, const tname & tval) {

#define END_OUT_OF_PLACE_SAVE() } }; } }


#endif  

#endif










