#ifndef GRAPHLAB_MR_DISK_GRAPH_CONSTRUCTION_IMPL_HPP
#define GRAPHLAB_MR_DISK_GRAPH_CONSTRUCTION_IMPL_HPP
#include <map>
#include <omp.h>
#include <graphlab/graph/disk_graph.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {
  
namespace mr_disk_graph_construction_impl {

vertex_id_t vertex_key_to_id(std::string s) {
  vertex_id_t vid;
  decompress_int<vertex_id_t>(s.c_str() + 1, vid);
  return vid;
}

std::pair<vertex_id_t, vertex_id_t> edge_key_to_id(std::string s) {
  ASSERT_EQ(s[0], 'e');
  
  const char *c = s.c_str() + 1;
  std::pair<vertex_id_t, vertex_id_t> edge;
  decompress_int_from_ref<vertex_id_t>(c, edge.first);
  ++c;
  decompress_int_from_ref<vertex_id_t>(c, edge.second);
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
                              size_t idx) {
  std::vector<disk_atom*> atoms;
  atoms.resize(disk_atom_files.size());
  
  // open the atoms 
  for (size_t i = 0;i < disk_atom_files.size(); ++i) {
    atoms[i] = new disk_atom(disk_atom_files[i], idx);
  }

  // create the output store
  disk_atom atomout(output_disk_atom, idx);
  atomout.clear();

  volatile uint32_t max_color = 0;
  // iterate through each database, joining the keys as we see it
  #pragma omp parallel for
  for (size_t i = 0;i < atoms.size(); ++i) {
    // open a cursor
    std::string key, val;
    kyotocabinet::HashDB::Cursor* cur = atoms[i]->get_db().cursor();
    cur->jump();
    
    // begin iteration
    while (cur->get(&key, &val, true)) {
      ASSERT_GT(key.length(), 0);
      char c = key[0];
      
      // we only need to track 4 entries
      // v (vertex), e (edge) , c (color) and h (dht). 

      if (c == 'v') {
        // vertex.
        boost::iostreams::stream<boost::iostreams::array_source> 
                                    istrm(val.c_str(), val.length());   
        iarchive iarc(istrm);
        uint16_t owner;
        iarc >> owner;
        vertex_id_t vid = vertex_key_to_id(key);
        // if has data, set the data
        istrm.peek();
        if (!istrm.eof()) {
          VertexData vdata;
          iarc >> vdata;
          ASSERT_EQ(owner, idx);
          // this will overwrite all other "skipped" vertices
          atomout.add_vertex(vid, owner, vdata);
        }
        else {
          atomout.add_vertex_skip(vid, owner);
        }
      }
      else if (c == 'e') {
        // edge
        std::pair<vertex_id_t, vertex_id_t> edge = edge_key_to_id(key);
        if (val.length() > 0) {
          EdgeData edata;
          boost::iostreams::stream<boost::iostreams::array_source> 
                                    istrm(val.c_str(), val.length());   
          iarchive iarc(istrm);
          iarc >> edata;
          atomout.add_edge(edge.first, edge.second, edata);
          atomout.inc_numlocale();
        }
        else {
          atomout.add_edge_skip(edge.first, edge.second);
        }
      }
      else if (c == 'c') {
        // color
        vertex_id_t vid = vertex_key_to_id(key);       
        uint32_t color = *reinterpret_cast<const uint32_t*>(val.c_str());
        atomout.set_color(vid, color);
        while(max_color < color) {
          uint32_t old_max_color = max_color;
          uint32_t new_max_color = std::max(color, old_max_color);
          // no op if no change
          if (new_max_color == old_max_color || 
              atomic_compare_and_swap(max_color, old_max_color, new_max_color)) {
            break;
          }
        }
      }
      else if (c == 'h') {
        vertex_id_t vid = vertex_key_to_id(key);       
        uint16_t owner= *reinterpret_cast<const uint16_t*>(val.c_str());
        atomout.set_owner(vid, owner);
      }
    }
    delete cur;
  }
  atomout.synchronize();
  atom_properties ret;
  ret.adjacent_atoms = atomout.enumerate_adjacent_atoms();
  ret.num_local_vertices = atomout.num_local_vertices();
  ret.num_local_edges = atomout.num_local_edges();
  ret.max_color = max_color;
  ret.filename = output_disk_atom;
  //std::cout << idx << " " << ret.num_local_vertices << " " << ret.num_local_edges << "\n";

  for (size_t i = 0;i < disk_atom_files.size(); ++i) {
    delete atoms[i];
  }
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
