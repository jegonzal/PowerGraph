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

#ifndef GRAPHLAB_SERIALIZE_LIST_HPP
#define GRAPHLAB_SERIALIZE_LIST_HPP

#include <list>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>


namespace graphlab {
namespace archive_detail {
  /** serializes a list  */
  template <typename ArcType, typename T>
  struct serialize_impl<ArcType, std::list<T>, false > {
  static void exec(ArcType& a, const std::list<T>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a list  */
  template <typename ArcType, typename T>
  struct deserialize_impl<ArcType, std::list<T>, false > {
  static void exec(ArcType& a, std::list<T>& vec){
    vec.clear();
    deserialize_iterator<T>(a, std::inserter(vec,vec.end()));
  }
  };
} // archive_detail  
} // graphlab
#endif 
