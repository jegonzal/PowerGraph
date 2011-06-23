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


#ifndef GRAPHLAB_HAS_LOAD_HPP
#define GRAPHLAB_HAS_LOAD_HPP

#include <typeinfo>

/**
   Detects if a class has a load function implemented
*/
namespace graphlab {
  namespace archive_detail {

    // SFINAE method derived from
    // http://stackoverflow.com/questions/87372/is-there-a-technique-in-c-to-know-if-a-class-has-a-member-function-of-a-given-s/87846#87846
    template<typename ArcType, typename T>
    struct has_load_method
    {
      template<typename U, void (U::*)(ArcType&)> struct SFINAE {};
      template<typename U> static char Test(SFINAE<U, &U::load>*);
      template<typename U> static int Test(...);
      static const bool value = sizeof(Test<T>(0)) == sizeof(char);
    };


    template <typename ArcType, typename ValueType>
    typename boost::enable_if_c<has_load_method<ArcType, ValueType>::value, void>::type 
    load_or_fail(ArcType& o, ValueType &t) { 
      t.load(o);
    }
  
    template <typename ArcType, typename ValueType>
    typename boost::disable_if_c<has_load_method<ArcType, ValueType>::value, void>::type 
    load_or_fail(ArcType& o, ValueType &t) { 
      ASSERT_MSG(false, "Trying to deserializable type %s without valid load method.", typeid(ValueType).name()); 
    }
  
  }  // archive_detail
}  // graphlab

#endif

