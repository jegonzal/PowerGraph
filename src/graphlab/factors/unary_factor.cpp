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

#include <graphlab/factors/unary_factor.hpp>

std::ostream& operator<<(std::ostream& out, 
                         const graphlab::unary_factor& fact) {
  out << "Unary Factor(" << fact.arity()  << ")"
      << std::endl;
  for(size_t i = 0; i < fact.arity(); ++i) {
    out << fact.logP(i) << " ";
  }
  out << std::endl;
  return out;
} // end of operator<<

