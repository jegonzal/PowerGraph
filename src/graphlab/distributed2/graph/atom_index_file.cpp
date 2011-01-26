#include <fstream>
#include <graphlab/distributed2/atom_index_fie.hpp>

namespace graphlab {

atom_index_file read_atom_index(std::string indexfile) {
  atom_index_file idxfile;
  // open file
  std::ifstream fin(indexfile.c_str());
  assert(fin.good());
  size_t natoms;
  fin >> idxfile.nverts >> idxfile.nedges >> natoms;
  for (size_t i = 0;i < natoms; ++i) {
    
  }
}


std::vector<std::vector<size_t> >
  partition_atoms(const atom_index_file& atomindex);




}
