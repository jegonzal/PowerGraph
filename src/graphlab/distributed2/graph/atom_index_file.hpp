#ifndef GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE
#define GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE
#include <string>
#include <vector>
#include <graphlab/graph/graph.hpp>

namespace graphlab {


struct atom_file{
  std::string protocol;
  std::string file;
  vertex_id_t nverts;
  edge_id_t nedges;
  std::vector<vertex_id_t> adjatoms; 
};

struct atom_index_file {
  size_t nverts, nedges, natoms;
  std::vector<atom_file> atoms;
};

atom_index_file read_atom_index(std::string indexfile, size_t nparts);
std::vector<std::vector<size_t> >
  partition_atoms(const atom_index_file& atomindex);

} // end namespace graphlab

#endif
