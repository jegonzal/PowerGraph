/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef GRAPHLAB_MR_DISK_GRAPH_CONSTRUCTION_IMPL_HPP
#define GRAPHLAB_MR_DISK_GRAPH_CONSTRUCTION_IMPL_HPP
#include <map>
#include <omp.h>
#include <graphlab/graph/disk_graph.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
  namespace mr_disk_graph_construction_impl {
    typedef disk_graph<int,int>::vertex_id_type vertex_id_type;
   

    vertex_id_type vertex_key_to_id(std::string s) {
      vertex_id_type vid;
      decompress_int<vertex_id_type>(s.c_str() + 1, vid);
      return vid;
    }

    std::pair<vertex_id_type, vertex_id_type> edge_key_to_id(std::string s) {
      ASSERT_EQ(s[0], 'e');
  
      const char *c = s.c_str() + 1;
      std::pair<vertex_id_type, vertex_id_type> edge;
      decompress_int_from_ref<vertex_id_type>(c, edge.first);
      ++c;
      decompress_int_from_ref<vertex_id_type>(c, edge.second);
      return edge;
    }

    struct atom_properties{
      size_t num_local_vertices;
      size_t num_local_edges;
      size_t max_color;
      std::string filename;
      std::map<uint16_t, uint32_t>  adjacent_atoms;
      void save(oarchive& oarc) const {
        oarc << num_local_vertices << num_local_edges
             << max_color << filename << adjacent_atoms;
      }
      void load(iarchive& iarc) {
        iarc >> num_local_vertices >> num_local_edges
             >> max_color >> filename >> adjacent_atoms;
      }
    };


    template <typename VertexData, typename EdgeData>
    atom_properties merge_parallel_disk_atom(std::vector<std::string> disk_atom_files, 
                                             std::string output_disk_atom,
                                             size_t idx,
                                             disk_graph_atom_type::atom_type atomtype) {
      std::vector<write_only_disk_atom*> atoms;
      atoms.resize(disk_atom_files.size());
  
      // open the atoms 
      for (size_t i = 0;i < disk_atom_files.size(); ++i) {
        atoms[i] = new write_only_disk_atom(disk_atom_files[i], idx);
      }

      // create the output store
      graph_atom* atomout = NULL;
      if (atomtype == disk_graph_atom_type::MEMORY_ATOM) {
        output_disk_atom += ".fast";
        atomout = new memory_atom(output_disk_atom, idx);
      }
      else if (atomtype == disk_graph_atom_type::WRITE_ONLY_ATOM) {
        output_disk_atom += ".dump";
        atomout = new write_only_disk_atom(output_disk_atom, idx);
      }
      else if (atomtype == disk_graph_atom_type::DISK_ATOM) {
        atomout = new disk_atom(output_disk_atom, idx);
      }

      atomout->clear();

      volatile uint32_t max_color = 0;
      // iterate through each database, joining the keys as we see it
#pragma omp parallel for
      for (int i = 0;i < (int)atoms.size(); ++i) {
        atoms[i]->play_back(atomout);
      }
      atomout->synchronize();
      atom_properties ret;
      ret.adjacent_atoms = atomout->enumerate_adjacent_atoms();
      ret.num_local_vertices = atomout->num_local_vertices();
      ret.num_local_edges = atomout->num_local_edges();
      ret.max_color = max_color;
      ret.filename = output_disk_atom;
      //std::cout << idx << " " << ret.num_local_vertices << " " << ret.num_local_edges << "\n";

      for (size_t i = 0;i < disk_atom_files.size(); ++i) {
        delete atoms[i];
      }
      delete atomout;
      return ret;
    }

    atom_index_file atom_index_from_properties(const std::map<size_t, atom_properties> &atomprops) {
      atom_index_file idx;
      std::map<size_t, atom_properties>::const_iterator iter = atomprops.begin();
      // compute the graph parameters
      idx.nverts = 0, idx.nedges = 0, idx.ncolors = 0;
      idx.natoms = atomprops.size();
      idx.atoms.resize(atomprops.size());
      while (iter != atomprops.end()) {
        idx.nverts += iter->second.num_local_vertices;
        idx.nedges += iter->second.num_local_edges;
        idx.ncolors = std::max<size_t>(idx.ncolors, iter->second.max_color);
        // i is the current atom index
        size_t i = iter->first;
        idx.atoms[i].protocol = "file";
        idx.atoms[i].file = iter->second.filename;
        idx.atoms[i].nverts = iter->second.num_local_vertices;
        idx.atoms[i].nedges = iter->second.num_local_edges;
        std::map<uint16_t, uint32_t>::const_iterator iteradj = iter->second.adjacent_atoms.begin();
        while (iteradj != iter->second.adjacent_atoms.end()) {
          idx.atoms[i].adjatoms.push_back(iteradj->first);
          idx.atoms[i].optional_weight_to_adjatoms.push_back(iteradj->second);
          ++iteradj;
        }
        ++iter;
      }
      return idx;
  
    }

  }   // namespace 

} // namespace graphlab
#endif

