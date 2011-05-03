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

#ifndef GRAPHLAB_ROBUST_CAST_HPP
#define GRAPHLAB_ROBUST_CAST_HPP

#include <boost/utility.hpp>
#include <boost/type_traits/is_convertible.hpp>
namespace graphlab {
  /** robust_cast performs a static cast from type A to type B
      if a cast can be done. Return B() otherwise */
  
  template <typename Target, typename Source>
  typename boost::disable_if_c<boost::is_convertible<Source, Target>::value, 
                               Target>::type
                               robust_cast(const Source &h) {
    return Target();
  }
  
  template <typename Target, typename Source>
  typename boost::enable_if_c<boost::is_convertible<Source, Target>::value, 
                              Target>::type
                              robust_cast(const Source &h) {
    return (Target)h;
  }
}

#endif


