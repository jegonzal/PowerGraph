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

#ifndef GRAPHLAB_SERIALIZE_MAP_HPP
#define GRAPHLAB_SERIALIZE_MAP_HPP

#include <map>
#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/serialization/iterator.hpp>

namespace graphlab {

namespace archive_detail {
  /** Serializes a map */
  template <typename ArcType, typename T, typename U>
  struct serialize_impl<ArcType, std::map<T,U>, false > {
  static void exec(ArcType& a, const std::map<T,U>& vec){
    serialize_iterator(a,vec.begin(),vec.end(), vec.size());
  }
  };

  /** deserializes a map  */
      
  template <typename ArcType, typename T, typename U>
  struct deserialize_impl<ArcType, std::map<T,U>, false > {
  static void exec(ArcType& a, std::map<T,U>& vec){
    vec.clear();
    deserialize_iterator<ArcType, std::pair<T,U> >(a, std::inserter(vec,vec.end()));
  }
  };

} // archive_detail  
} // graphlab
#endif 
