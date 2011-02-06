#include <fstream>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>
#include <graphlab/logger/assertions.hpp>

namespace graphlab {

  
void atom_index_file::read_from_file(std::string indexfile) {
  // open file
  std::ifstream fin(indexfile.c_str());
  assert(fin.good());
/*
 [nverts] [nedges] [natoms]
 */
  fin >> nverts >> nedges >> natoms;
  
  
  for (size_t i = 0;i < natoms; ++i) {
    //   [atom_nverts] [atom_nedges]  [#adjacent atoms]  
    //           [space delimited list of adjacent atoms]    [atom_file_name]
    // [space delimited list of adjacent atoms] could be of 2 forms
    // just a list:    1 3 5 10 2
    // with weights:   1:100 3:50: 5:19 10:105: 2:99   (format is atomid:weight)
    //
    // if no weight attached, default is 1

    atom_file_descriptor atom;
    fin >> atom.nverts >> atom.nedges;
    size_t numadj;
    fin >> numadj;
    for (size_t j = 0; j < numadj; ++j) {
      std::string tonextspace;
      getline(fin, tonextspace, ' ');
      tonextspace = trim(tonextspace);
      // check for a colon
      size_t colonpos = tonextspace.find(':');
      if (colonpos == std::string::npos) {
        // if no colon
        atom.adjatoms.push_back(atoi(tonextspace.c_str()));
        atom.optional_weight_to_adjatoms.push_back(1);
      }
      else {
        // split at colon
        tonextspace[colonpos] = 0;
        atom.adjatoms.push_back(atoi(tonextspace.c_str()));
        atom.optional_weight_to_adjatoms.push_back(atoi(&(tonextspace.c_str()[colonpos+1])));
      }
    }
    
    std::string afile;
    getline(fin, afile);
    afile = trim(afile);
    // read the protocol
    size_t protsep = afile.find_first_of("://");
    assert(protsep != std::string::npos);
    atom.protocol = afile.substr(0, protsep);
    atom.file = afile.substr(protsep + 3);
    atoms.push_back(atom);
  }
}


void atom_index_file::write_to_file(std::string outfilename) {
  // check if all weights are unit. and while we are at it
  // make sure that that the optional_weight_to_adjatoms is the right length
  // (either 0 or same length as adjatoms)
  bool allweightsunit = true;
  for (size_t i = 0;i < atoms.size(); ++i) {
    assert(atoms[i].optional_weight_to_adjatoms.size() == 0 ||
           atoms[i].optional_weight_to_adjatoms.size() == atoms[i].adjatoms.size());
    
    for (size_t j = 0;j < atoms[i].optional_weight_to_adjatoms.size(); ++j) {
      if (atoms[i].optional_weight_to_adjatoms[j] != 1) {
        allweightsunit = false;
        break;
      }
    }
    if (allweightsunit == false) break;
  }
  // if all weights are unity, clear the optional_weight_to_adjatoms table
  if (allweightsunit) {
    for (size_t i = 0; i < atoms.size(); ++i) atoms[i].optional_weight_to_adjatoms.clear();
  }

  // write file
  std::ofstream fout(outfilename.c_str());
  assert(fout.good());
  fout << nverts << " " << nedges << " " << natoms << "\n";
  for (size_t i = 0;i < natoms; ++i) {
    fout << atoms[i].nverts << "\t" << atoms[i].nedges << "\t"
         << atoms[i].adjatoms.size() << "\t";
    for (size_t j = 0;j < atoms[i].adjatoms.size(); ++j) {
      fout << atoms[i].adjatoms[j];
      if (atoms[i].optional_weight_to_adjatoms.size() != 0) {
        fout << ":" << atoms[i].optional_weight_to_adjatoms[j];
      }
      fout << " ";
    }
    fout << "\t" << atoms[i].protocol << "://" << atoms[i].file << "\n";
  }
}

size_t identity_function(const size_t &val) {
  return val;
}

size_t unity_function(const size_t &val) {
  return 1;
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
      size_t atomweight = 1;
      if (atomindex.atoms[i].optional_weight_to_adjatoms.size() != 0) {
        atomweight = atomindex.atoms[i].optional_weight_to_adjatoms[j];
      }
      atomgraph.add_edge(i, atomindex.atoms[i].adjatoms[j], atomweight);
    }
  }
  
  std::vector<uint32_t> retpart;
  atomgraph.metis_weighted_partition(nparts, retpart,
                                     unity_function, identity_function);
  //atomgraph.metis_partition(nparts, retpart);
  std::vector<std::vector<size_t> > ret;
  ret.resize(nparts);
  for (size_t i = 0;i < retpart.size(); ++i) {
    ret[retpart[i]].push_back(i);
  }
  return ret;

}



} // end namespace graphlab
