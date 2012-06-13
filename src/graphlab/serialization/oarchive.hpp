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
    std::ostream* out;
    /// constructor. Takes a generic std::ostream object
    inline oarchive(std::ostream& outstream)
      : out(&outstream) {}

    inline void write(const char* c, std::streamsize s) {
      out->write(c, s);
    }

    inline bool fail() {
      return out->fail();
    }
    
    inline ~oarchive() { }
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
    std::ostream* out;

    inline oarchive_soft_fail(std::ostream& outstream)
      : out(&outstream) {}

    inline oarchive_soft_fail(oarchive& oarc):out(oarc.out) {}
    
    inline void write(const char* c, std::streamsize s) {
      out->write(c, s);
    }
 
    inline bool fail() {
      return out->fail();
    }
    
    inline ~oarchive_soft_fail() { }
  };

  namespace archive_detail {

    /// called by the regular archive The regular archive will do a hard fail
    template <typename OutArcType, typename T>
    struct serialize_hard_or_soft_fail {
      inline static void exec(OutArcType& oarc, const T& t) {
        t.save(oarc);
      }
    };

    /// called by the soft fail archive 
    template <typename T>
    struct serialize_hard_or_soft_fail<oarchive_soft_fail, T> {
      inline static void exec(oarchive_soft_fail& oarc, const T& t) {
        // create a regular oarchive and
        // use the save_or_fail function which will
        // perform a soft fail
        oarchive regular_oarc(*(oarc.out));
        save_or_fail(regular_oarc, t);
      }
    };


    /**
       Implementation of the serializer for different types.
       This is the catch-all. If it gets here, it must be a non-POD and is a class.
       We therefore call the .save function.
       Here we pick between the archive types using serialize_hard_or_soft_fail
    */
    template <typename OutArcType, typename T, bool IsPOD>
    struct serialize_impl {
      static void exec(OutArcType& oarc, const T& t) {
        serialize_hard_or_soft_fail<OutArcType, T>::exec(oarc, t);
      }
    };

    /** Catch if type is a POD */
    template <typename OutArcType, typename T>
    struct serialize_impl<OutArcType, T, true> {
      inline static void exec(OutArcType& oarc, const T& t) {
        oarc.write(reinterpret_cast<const char*>(&t), sizeof(T));
      }
    };

    /**
       Re-dispatch if for some reasons T already has a const
    */
    template <typename OutArcType, typename T>
    struct serialize_impl<OutArcType, const T, true> {
      inline static void exec(OutArcType& oarc, const T& t) {
        serialize_impl<OutArcType, T, true>::exec(oarc, t);
      }
    };
    
    /**
       Re-dispatch if for some reasons T already has a const
    */
    template <typename OutArcType, typename T>
    struct serialize_impl<OutArcType, const T, false> {
      inline static void exec(OutArcType& oarc, const T& t) {
        serialize_impl<OutArcType, T, false>::exec(oarc, t);
      }
    };
  }// archive_detail


  /**
     Overloads the operator<< in the oarchive to
     allow the use of the stream syntax for serialization.
     It simply re-dispatches into the serialize_impl classes 
  */
  template <typename T>
  inline oarchive& operator<<(oarchive& oarc, const T& t) {
    archive_detail::serialize_impl<oarchive, 
                                   T, 
                                   gl_is_pod<T>::value >::exec(oarc, t);
    return oarc;
  }

  /**
     Overloads the operator<< in the oarchive_soft_fail to
     allow the use of the stream syntax for serialization.
     It simply re-dispatches into the serialize_impl classes 
  */
  template <typename T>
  inline oarchive_soft_fail& operator<<(oarchive_soft_fail& oarc, 
                                        const T& t) {
    archive_detail::serialize_impl<oarchive_soft_fail, 
                                  T, 
                                  gl_is_pod<T>::value >::exec(oarc, t);
    return oarc;
  }


  /**
     Serializes an arbitrary pointer + length to an archive 
  */
  inline oarchive& serialize(oarchive& oarc, 
                             const void* str,
                             const size_t length) {
    // save the length
    operator<<(oarc,length);
    oarc.write(reinterpret_cast<const char*>(str), 
                    (std::streamsize)length);
    assert(!oarc.fail());
    return oarc;
  }


  /**
     Serializes an arbitrary pointer + length to an archive 
  */
  inline oarchive_soft_fail& serialize(oarchive_soft_fail& oarc, 
                                       const void* str,
                                       const size_t length) {
    // save the length
    operator<<(oarc,length);
    oarc.write(reinterpret_cast<const char*>(str), 
                    (std::streamsize)length);
    assert(!oarc.fail());
    return oarc;
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
  template <typename OutArcType> struct serialize_impl<OutArcType, tname, false> { \
  static void exec(OutArcType& arc, const tname & tval) {

#define END_OUT_OF_PLACE_SAVE() } }; } }


#endif  

#endif










