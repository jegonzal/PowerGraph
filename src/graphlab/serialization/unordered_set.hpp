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

#ifndef GRAPHLAB_SERIALIZE_UNORDERED_SET_HPP
#define GRAPHLAB_SERIALIZE_UNORDERED_SET_HPP

#include <boost/unordered_set.hpp>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

namespace graphlab {
namespace archive_detail {
  /** serializes a set  */
  template <typename ArcType, typename T>
  struct serialize_impl<ArcType, boost::unordered_set<T>, false > {
  static void exec(ArcType& a, const boost::unordered_set<T>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a set  */
  template <typename ArcType, typename T>
  struct deserialize_impl<ArcType, boost::unordered_set<T>, false > {
  static void exec(ArcType& a, boost::unordered_set<T>& vec){
    vec.clear();
    deserialize_iterator<ArcType, T>(a, std::inserter(vec,vec.end()));
  }
  };

} // archive_detail  
} // graphlab

#endif 
