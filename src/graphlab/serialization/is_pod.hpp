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

#ifndef GRAPHLAB_IS_POD_HPP
#define GRAPHLAB_IS_POD_HPP
#include <boost/type_traits.hpp>

namespace graphlab {

  template <typename T>
  struct gl_is_pod{
    // it is a pod and is not an integer since we have special handlings for integers

    // (T is POD and  T is not an integer of size >= 2)
    BOOST_STATIC_CONSTANT(bool, value =
                          (
                           boost::type_traits::ice_and<
                             boost::is_pod<T>::value,
                             boost::type_traits::ice_not<
                               boost::type_traits::ice_and<
                                 boost::is_integral<T>::value,
                                 sizeof(T) >= 2
                                 >::value
                               >::value
                             >::value
                          ));

  };

}

#endif



