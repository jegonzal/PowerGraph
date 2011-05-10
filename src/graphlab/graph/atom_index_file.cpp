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

#include <fstream>
#include <algorithm>
#include <vector>
#include <set>


#include <boost/filesystem.hpp>

#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/fs_util.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/graph/atom_index_file.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

#include <graphlab/macros_def.hpp>


namespace graphlab {
  
  void atom_index_file::read_from_file(std::string indexfile) {
    // open file
    std::ifstream fin(indexfile.c_str());
    assert(fin.good());
    /*
      [nverts] [nedges] [natoms]
    */
    std::string label;
    size_t val;
    for (size_t i = 0;i < 4; ++i) {
      fin >> label >> val;
      if (trim(label) == "NumVerts:") {
        nverts = val;
      }
      else if (trim(label) == "NumEdges:") {
        nedges = val;
      }
      else if (trim(label) == "NumColors:") {
        ncolors = val;
      }
      else if (trim(label) == "NumAtoms:") {
        natoms = val;
      }
      else {
        logger(LOG_ERROR, "Unrecognized label in index file: %s", label.c_str());
      }
    }
  
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
          atom.optional_weight_to_adjatoms.
            push_back(atoi(&(tonextspace.c_str()[colonpos+1])));
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
      for (size_t i = 0; i < atoms.size(); ++i) 
        atoms[i].optional_weight_to_adjatoms.clear();
    }

    // write file
    std::ofstream fout(outfilename.c_str());
    assert(fout.good());
    fout << "NumVerts:\t" << nverts << "\n" 
         << "NumEdges:\t" << nedges << "\n" 
         << "NumColors:\t" << ncolors << "\n" 
         << "NumAtoms:\t" << natoms << "\n";
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

  /**
   * Takes in a subset of atoms and bisect its atom subgraph using metis
   */
  std::vector<std::vector<size_t> >
  bisect_atoms_metis(const atom_index_file& atomindex, 
                     std::vector<size_t> atomsubset) {
    // build a graph using only the atomsubset
    std::set<size_t> atomsubset_set;
    std::map<size_t, size_t> atomrevmap;
    std::copy(atomsubset.begin(), 
              atomsubset.end(), 
              std::inserter(atomsubset_set,atomsubset_set.end()));

  
    graph<size_t, size_t> atomgraph;
    // add vertices
    foreach(size_t i, atomsubset) {
      atomrevmap[i] = atomgraph.add_vertex(atomindex.atoms[i].nedges + 1);
    }
  
    foreach(size_t i, atomsubset) {
      for (size_t j = 0; j < atomindex.atoms[i].adjatoms.size(); ++j) {
        // only add edges which connect within the atomsubset
        if (atomsubset_set.find(atomindex.atoms[i].adjatoms[j]) != atomsubset_set.end()) {
          size_t nbr = atomrevmap[atomindex.atoms[i].adjatoms[j]];
          size_t atomweight = 1;
          if (atomindex.atoms[i].optional_weight_to_adjatoms.size() != 0) {
            atomweight = atomindex.atoms[i].optional_weight_to_adjatoms[j];
          }
          if (nbr == atomrevmap[i]) continue;
          atomgraph.add_edge(atomrevmap[i], nbr, atomweight);
          atomgraph.add_edge(nbr, atomrevmap[i], atomweight);
        }
      }
    }
    std::vector<uint32_t> retpart;
    atomgraph.metis_weighted_partition(2, retpart,
                                       identity_function, identity_function, true);

    std::vector<std::vector<size_t> > ret(2);
    for (size_t i = 0;i < retpart.size(); ++i) {
      ret[retpart[i]].push_back(atomsubset[i]);
    }
    return ret;
  }


  /**
   * Takes in a subset of atoms and bisect it down the middle
   */
  std::vector<std::vector<size_t> > 
  bisect_atoms_heuristic(const atom_index_file& atomindex, 
                         std::vector<size_t> atomsubset) {
    std::vector<std::vector<size_t> > ret(2);
    for (size_t i = 0;i < atomsubset.size(); ++i) {
      ret[int(i >= (atomsubset.size() / 2))].push_back(atomsubset[i]);
    }

    // we could do a simple heuristic refinement here.
    return ret;
  }

  std::vector<std::vector<size_t> > 
  partition_atoms_sliced(const atom_index_file& atomindex, size_t nparts) {
    std::vector<std::vector<size_t> > ret;
    ret.resize(nparts);
    for (size_t i = 0;i < atomindex.atoms.size(); ++i) {
      ret[i % nparts].push_back(i);
    }
    return ret;
  }

  std::vector<std::vector<size_t> > 
  partition_atoms(const atom_index_file& atomindex, size_t nparts) {
    // build the atom graph
    // vertex weight is #edges
    // edge weight is 1
    // I cannot ask for more parts tan atom
    ASSERT_LE(nparts, atomindex.atoms.size());
  
    graph<size_t, size_t> atomgraph; 
    for (size_t i = 0;i < atomindex.atoms.size(); ++i) {
      atomgraph.add_vertex(atomindex.atoms[i].nedges + 1);
    }
    for (size_t i = 0;i < atomindex.atoms.size(); ++i) {
      for (size_t j = 0; j < atomindex.atoms[i].adjatoms.size(); ++j) {
        size_t atomweight = 1;
        if (atomindex.atoms[i].optional_weight_to_adjatoms.size() != 0) {
          atomweight = atomindex.atoms[i].optional_weight_to_adjatoms[j];
        }
        if (atomindex.atoms[i].adjatoms[j] == i) continue;
        atomgraph.add_edge(i, atomindex.atoms[i].adjatoms[j], atomweight);
        atomgraph.add_edge(atomindex.atoms[i].adjatoms[j], i, atomweight);
      }
    }
    
    std::vector<uint32_t> retpart; 
    if (atomgraph.num_edges() >= atomgraph.num_vertices() * atomgraph.num_vertices() / 2) {
      logstream(LOG_WARNING) << "high atom edge density. Using random partition" << std::endl;
      atomgraph.random_partition(nparts, retpart);
    
    }
    else {
      //  std::cout << atomgraph;
      atomgraph.metis_weighted_partition(nparts, retpart,
                                         identity_function, identity_function, true);
      //atomgraph.metis_partition(nparts, retpart);
    }
    std::vector<std::vector<size_t> > ret;
    ret.resize(nparts);
    for (size_t i = 0;i < retpart.size(); ++i) {
      ret[retpart[i]].push_back(i);
    }
    // compute the current weight of each part
    std::vector<size_t> partweights(nparts, 0);
    size_t totalweight = 0;

    mutable_queue<size_t, size_t> partition_weights;
  
    for (size_t i = 0;i < ret.size(); ++i) {
      for (size_t j = 0; j < ret[i].size(); ++j) {
        partweights[i] += atomindex.atoms[ret[i][j]].nedges + 1;
        totalweight += atomindex.atoms[ret[i][j]].nedges + 1;
      }
      partition_weights.push(i, partweights[i]);
    }

    logstream(LOG_INFO) << "Balance factor (max / avg): " << float(partition_weights.top().second) /
      (float(totalweight) / nparts) << std::endl;


    // ok. the annoying part is that metis might actually give me a
    // smaller number of partitions.
    std::vector<size_t> missingparts;
    for (size_t i = 0;i < ret.size(); ++i) {
      if (ret[i].size() == 0)  missingparts.push_back(i);
    }


    if (missingparts.size() > 0) {
      logstream(LOG_WARNING) << "Metis generated only " << nparts - missingparts.size() <<
        " when " << nparts << " was requested. Attempting to repartition." << std::endl;
    }


    for (size_t i = 0;i < missingparts.size(); ++i) {
      // for each missing part, get the "heaviest"
      // partition and bisect it
      size_t repartidx = partition_weights.top().first;
      std::vector<std::vector<size_t> > bisect = bisect_atoms_metis(atomindex, ret[repartidx]);
      // if any the bisectinos are 0, use the heuristic bisection
      if (bisect[0].size() == 0 || bisect[1].size() == 0) {
        std::cout << "bisecting: " << repartidx << std::endl;
        bisect = bisect_atoms_heuristic(atomindex, ret[repartidx]);
      }
      // assert
      ASSERT_TRUE(bisect[0].size() != 0 && bisect[1].size() != 0);
      // update the partitions
      ret[repartidx] = bisect[0];
      ret[missingparts[i]] = bisect[1];
      // update the mutable queue
      partition_weights.update(repartidx, bisect[0].size());
      partition_weights.update(missingparts[0], bisect[1].size());
    }

    if (missingparts.size() > 0) {
      logstream(LOG_WARNING) << "Repartition successful. New balance factor (max / avg): " << float(partition_weights.top().second) /
        (float(totalweight) / nparts) << std::endl;

    }
    return ret;

  }



} // namespace graphlab








