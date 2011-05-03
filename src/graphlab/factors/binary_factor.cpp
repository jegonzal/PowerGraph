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

#include <graphlab/factors/binary_factor.hpp>

std::ostream& operator<<(std::ostream& out, 
                         const graphlab::binary_factor& fact) {
  out << "Binary Factor(v_" << fact.var1() << " in {1..."
      << fact.arity1() << "}, " 
      << ", v_ " << fact.var2() << " in {1..." 
      << fact.arity2() << "})" << std::endl;
  for(size_t i = 0; i < fact.arity1(); ++i) {
    for(size_t j = 0; j < fact.arity2(); ++j) {
      out << fact.logP(i,j) << " ";
    }
    out << std::endl;
  }
  return out;
} // end of operator<<
