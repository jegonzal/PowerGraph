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

#include <vector>
#include <graphlab/distributed2/distributed_glshared.hpp>

namespace graphlab {
namespace distgl_impl {
	

std::vector<distributed_glshared_base*>& get_global_dist_glshared_registry() {
  static std::vector<distributed_glshared_base*> global_registry;
  return global_registry;
}	

void register_dist_glshared(distributed_glshared_base* glsharedobj) {
  get_global_dist_glshared_registry().push_back(glsharedobj);
}



} // namespace distgl_impl
} // namespace graphlab
