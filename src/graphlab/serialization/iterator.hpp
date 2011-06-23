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


#ifndef GRAPHLAB_SERIALIZE_ITERATOR_HPP
#define GRAPHLAB_SERIALIZE_ITERATOR_HPP

#include <iterator>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iarchive.hpp>

namespace graphlab {

  /**
    Serializes the contents between the iterators begin and end.
    This version prefers the availability of RandomAccessIterator since it needs
    a distance between the begin and end iterator.
    This function as implemented will work for other input iterators
    but is extremely inefficient.
    Returns true on success, false on failure  */
  template <typename ArcType, typename RandomAccessIterator>
  void serialize_iterator(ArcType& a, RandomAccessIterator begin,
                                      RandomAccessIterator end){
    size_t vsize = std::distance(begin, end);
    a << vsize;
    //store each element
    for(; begin != end; ++begin) {
      a << *begin;
    }
  }


  /**
    Serializes the contents between the iterators begin and end.
    This version takes all InputIterator types, but takes a "count" for
    efficiency. This count is checked and will return failure if the number
    of elements serialized does not match the count
    Returns true on success, false on failure  */
  template <typename ArcType, typename InputIterator>
  void serialize_iterator(ArcType& a, InputIterator begin,
                                      InputIterator end, size_t vsize){
    a << vsize;
    //store each element
    size_t count = 0;
    for(; begin != end; ++begin) {
      ++count;
      a << *begin;
    }
    // fail if count does not match
    assert(count == vsize);
  }

  /**
    The accompanying function to serialize_iterator()
    Reads elements from the stream and send it to the output iterator.
    Note that this requires an additional template parameter T which is the
    "type of object to deserialize"
    This is necessary for instance for the map type. The map<T,U>::value_type
    is pair<const T,U> which is not useful since I cannot assign to it.
    In this case, T=pair<T,U>

    Returns true on success, false on failure  */
  template <typename ArcType, typename T, typename OutputIterator>
  void deserialize_iterator(ArcType& a, OutputIterator result) {
    // get the number of elements to deserialize
    size_t length = 0;
    a >> length;
    
    // iterate through and send to the output iterator
    for (size_t x = 0; x < length ; ++x){
      T v;
      a >> v;
      (*result) = v;
      result++;
    }
  }
  
 
} // namespace prl
#endif //PRL_SERIALIZE_ITERATOR_HPP

