#include <fstream>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>

namespace graphlab {

  
atom_index_file read_atom_index(std::string indexfile) {
  
  atom_index_file idxfile;
  // open file
  std::ifstream fin(indexfile.c_str());
  assert(fin.good());
  size_t natoms;
/*
 [nverts] [nedges] [natoms]
 */
  fin >> idxfile.nverts >> idxfile.nedges >> natoms;
  
  
  for (size_t i = 0;i < natoms; ++i) {
    //   [atom_nverts] [atom_nedges]  [#adjacent atoms]  
    //           [space delimited list of adjacent atoms]    [atom_file_name] 

    atom_file atom;
    fin >> atom.nverts >> atom.nedges;
    size_t numadj;
    fin >> numadj;
    for (size_t j = 0; j < numadj; ++j) {
      size_t tmp;
      fin >> tmp;
      atom.adjatoms.push_back(tmp);
    }
    std::string afile;
    getline(fin, afile);
    afile = trim(afile);
    // read the protocol
    size_t protsep = afile.find_first_of("://");
    assert(protsep != std::string::npos);
    atom.protocol = afile.substr(0, protsep);
    atom.file = afile.substr(protsep + 3);
    idxfile.atoms.push_back(atom);
  }
  
  return idxfile;
}


size_t identity_function(const size_t &val) {
  return val;
}



std::vector<std::vector<size_t> >
partition_atoms(const atom_index_file& atomindex, size_t nparts) {
  // build the atom graph
  // vertex weight is #edges
  // edge weight is 1
  graph<size_t, size_t> atomgraph; 
  for (size_t i = 0;i < atomindex.atoms.size(); ++i) {
    atomgraph.add_vertex(atomindex.atoms[i].nedges);
  }
  for (size_t i = 0;i < atomindex.atoms.size(); ++i) {
    for (size_t j = 0; j < atomindex.atoms[i].adjatoms.size(); ++j) {
      atomgraph.add_edge(i, atomindex.atoms[i].adjatoms[j], 1);
    }
  }
  
  std::vector<uint32_t> retpart;
  atomgraph.metis_weighted_partition(nparts, retpart, 
                                     identity_function, identity_function);
  
  std::vector<std::vector<size_t> > ret;
  ret.resize(nparts);
  for (size_t i = 0;i < retpart.size(); ++i) {
    ret[retpart[i]].push_back(i);
  }
  return ret;

}



} // end namespace graphlab
