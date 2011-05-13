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

