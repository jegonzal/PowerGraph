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

#ifndef GRAPHLAB_FLOAT_SELECTOR_HPP
#define GRAPHLAB_FLOAT_SELECTOR_HPP

namespace graphlab {
  
  template <int len>
  struct float_selector {
    // invalid
  };


  template <>
  struct float_selector<4> {
    typedef float float_type;
  };

  template <>
  struct float_selector<8> {
    typedef double float_type;
  };

  template <>
  struct float_selector<16> {
    typedef long double float_type;
  };

}
#endif


