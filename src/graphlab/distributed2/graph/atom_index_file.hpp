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

#ifndef GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE
#define GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE

#include <string>
#include <vector>
#include <graphlab/graph/graph.hpp>

namespace graphlab {

std::vector<std::vector<size_t> > 
partition_atoms_sliced(const std::vector<std::set<uint16_t> > & atomadj, size_t nparts);

std::vector<std::vector<size_t> >
partition_atoms(const std::vector<std::set<uint16_t> > & atomadj, size_t nparts);


} // end namespace graphlab

#endif
