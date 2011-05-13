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
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/graph/graph.hpp>

namespace graphlab {


  struct atom_file_descriptor{
    std::string protocol;
    std::string file;
    vertex_id_t nverts;
    edge_id_t nedges;
    std::vector<vertex_id_t> adjatoms;
    std::vector<size_t> optional_weight_to_adjatoms;
    void save(oarchive &oarc) const{
      oarc << protocol
           << file
           << nverts
           << nedges
           << adjatoms
           << optional_weight_to_adjatoms;
    }
    
    void load(iarchive &iarc) {
      iarc >> protocol
           >> file
           >> nverts
           >> nedges
           >> adjatoms
           >> optional_weight_to_adjatoms;
    }
  };

  /**
  The atom index file is a file describing the locations and adjacency
  structure of a collection of atom files. It can be generated from
  a list of atom files using the disk_graph. The atom index file is stored
  as text for human readability.
  */
  struct atom_index_file {

    size_t nverts, nedges, ncolors, natoms;
    std::vector<atom_file_descriptor> atoms;

    void read_from_file(std::string indexfile);
    void write_to_file(std::string outfilename);
  };

  std::vector<std::vector<size_t> > 
  partition_atoms_sliced(const atom_index_file& atomindex, size_t nparts);

  std::vector<std::vector<size_t> >
  partition_atoms(const atom_index_file& atomindex, size_t nparts);


} // end namespace graphlab

#endif

