#include <fstream>
#include <algorithm>
#include <vector>
#include <set>


#include <boost/filesystem.hpp>

#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/fs_util.hpp>
#include <graphlab/util/mutable_queue.hpp>
#include <graphlab/distributed2/graph/atom_file.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>


#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>


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
    //  std::cout << atomgraph;
    std::vector<uint32_t> retpart;
    atomgraph.metis_weighted_partition(nparts, retpart,
                                       identity_function, identity_function, true);
    //atomgraph.metis_partition(nparts, retpart);
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














  // /** This function computes the collection of local files to be
  //     loaded */
  // void find_local_atom_files(const std::string& pathname,
  //                            std::vector<std::string>& files) {
  //   namespace fs = boost::filesystem;
  //   fs::path path(pathname);
  //   assert(fs::exists(path));
  //   for(fs::directory_iterator iter( path ), end_iter; 
  //       iter != end_iter; ++iter) {
  //     if( ! fs::is_directory(iter->status()) ) {
  //       std::string filename(iter->path().filename());
  //       size_t period_loc = filename.rfind('.');
  //       if(period_loc != std::string::npos) {
  //         std::string ending(filename.substr(period_loc));
  //         if(ending == ".atom") {
  //           files.push_back(iter->path().filename());
  //         }
  //       }
  //     }
  //   }
  //   std::sort(files.begin(), files.end());
  // }




  struct atom_file_desc_extra {
    atom_file_descriptor desc;
    edge_id_t nlocalverts;
    edge_id_t ninedges;
    size_t max_color;
    atom_file_desc_extra() : 
      nlocalverts(0), ninedges(0) { }
    
    void save(oarchive &oarc) const{
      oarc << desc
           << nlocalverts
           << ninedges
           << max_color;
    }
  
    void load(iarchive &iarc) {
      iarc >> desc
           >> nlocalverts
           >> ninedges
           >> max_color;
    }
  };




  void build_atom_index_file(const std::string& path,
                             const std::string& aindex_fname) {
    // Create a dc comm layer
    logstream(LOG_INFO)
      << "Initializing distributed communication layer." << std::endl;
    dc_init_param param;         
    if( ! init_param_from_mpi(param) ) {
      logstream(LOG_FATAL) 
        << "Failed MPI laucher!" << std::endl;
      // if ( !init_param_from_env(param) ) {
      //   std::cout << "Failed to get environment variables" << std::endl;
      //   assert(false);
      // }
    }
    distributed_control dc(param);

    dc.full_barrier();
  
    // load the vector of filenames
    std::vector<std::string> local_fnames; 
    size_t num_atoms = 0;
    {
      logstream(LOG_INFO) 
        << dc.procid() << ": " 
        <<  "Computing local filenames"
        << std::endl;
      fs_util::list_files_with_suffix(path, ".atom", local_fnames);
      logstream(LOG_INFO) 
        << dc.procid() << ": " 
        <<  "Found files: "
        << std::endl;
      for(size_t i = 0; i < local_fnames.size(); ++i) {
        logstream(LOG_INFO) 
          << dc.procid() << ": " 
          <<  "(" << local_fnames[i] << ")"
          << std::endl;        
      }
      // compute the fnames that are used by this machine
      std::vector< std::vector< std::string > > partition_fnames;
      dc.gather_partition(local_fnames, partition_fnames);
      num_atoms = 0;
      for(size_t i = 0; i < partition_fnames.size(); ++i) 
        num_atoms += partition_fnames[i].size();
      logstream(LOG_INFO)
        << dc.procid() << ": " 
        <<  " number of atoms:    " << num_atoms
        << std::endl;
      // update the local fnames      
      local_fnames = partition_fnames[dc.procid()];
      logstream(LOG_INFO) 
        << dc.procid() << ": " 
        <<  "Assigned Files: "
        << std::endl;
      for(size_t i = 0; i < local_fnames.size(); ++i) {
        logstream(LOG_INFO) 
          << dc.procid() << ": " 
          <<  "(" << local_fnames[i] << ")"
          << std::endl;        
      };
    }
  
    dc.full_barrier();

    typedef std::map<procid_t, atom_file_desc_extra> atomid2desc_type;
    atomid2desc_type atomid2desc;
    {
      logstream(LOG_INFO) 
        << dc.procid() << ": " 
        <<  "Processing local atom files: "
        << std::endl;


      // hack: only need the type to read the atom header. The body is
      // not needed here
      typedef atom_file<bool, bool> fake_atom_file_type;
      foreach(const std::string& fname, local_fnames) {
        std::string absfname = path + "/" + fname;
        fake_atom_file_type afile;
        afile.input_filename("file", absfname);
        afile.load_structure(); // we do not need to load all
        procid_t atom_id = afile.atom_id();
        atom_file_desc_extra& afd_extra(atomid2desc[atom_id]);
        atom_file_descriptor& afd(afd_extra.desc);
        afd.file = absfname;
        afd.protocol = "file";
        afd.nverts = afile.globalvids().size();
        afd.nedges = afile.edge_src_dest().size();      
        // Compute adjacency counts for adjacent atoms as well as the
        // number of non ghost (local) vertices
        std::map<procid_t, size_t> adjatom_count;
        // Compute local vertices as well as update the adjacent atom
        // count
        afd_extra.nlocalverts = 0;
        foreach(const procid_t& other_atom_id, afile.atom())  {
          if(other_atom_id == atom_id) ++afd_extra.nlocalverts;
          adjatom_count[other_atom_id]++;
        }
        // Compute max color
        afd_extra.max_color = 0;
        foreach(const vertex_color_type& vcolor, afile.vcolor())  {
          afd_extra.max_color = std::max(afd_extra.max_color, size_t(vcolor));
        }
        // count in edges
        afd_extra.ninedges = 0;
        typedef std::pair<procid_t, procid_t> edge_pair_type;
        foreach(const edge_pair_type& edge, afile.edge_src_dest()) {
          const vertex_id_t target_vid(edge.second);
          if(afile.atom().at(target_vid) == atom_id) 
            ++afd_extra.ninedges;
        }
        // Fill in afd
        typedef std::pair<procid_t, size_t> count_pair;
        foreach(const count_pair& pair, adjatom_count) {
          afd.adjatoms.push_back(pair.first);
          afd.optional_weight_to_adjatoms.push_back(pair.second);
        }
      } // end of foreach
    } //end of loop over local atoms
  
    dc.full_barrier();

    // Gather all the machine counts
    std::vector<atomid2desc_type> all_descriptors(dc.numprocs());
    all_descriptors.at(dc.procid()) = atomid2desc;
    // Gather to root machine
    const size_t ROOT_NODE(0);
    const bool IS_ROOT(dc.procid() == ROOT_NODE);
    dc.gather(all_descriptors, ROOT_NODE);
    if(IS_ROOT) {
      atom_index_file aif;
      aif.atoms.resize(num_atoms);
      aif.natoms = num_atoms;
      aif.nverts = 0;
      aif.nedges = 0;
      aif.ncolors = 0;
      foreach(const atomid2desc_type& a2d, all_descriptors) {
        typedef atomid2desc_type::value_type pair_type;
        foreach(const pair_type& pair, a2d) {
          procid_t atomid = pair.first;
          assert(atomid < aif.atoms.size());
          // check that this atom has not already been added
          assert(aif.atoms[atomid].protocol.empty());
          // add the atom
          aif.atoms[atomid] = pair.second.desc;
          aif.nverts += pair.second.nlocalverts;
          aif.nedges += pair.second.ninedges;
          // use ncolors to track the max color;
          aif.ncolors = std::max(aif.ncolors, pair.second.max_color);
        }
      }
      // convert max color to actual number of colors
      aif.ncolors++;
      // Write results to file
      aif.write_to_file(aindex_fname);
    }
  

  } // and of build atom index file












  // void build_atom_index_file(const std::string& path) {
  //   // Create a dc comm layer
  //   std::cout << "Initializing distributed communication layer" << std::endl;
  
  //   // load the vector of filenames
  //   std::vector<std::string> local_fnames; 
  
  //   std::cout << "Computing local filenames"
  //             << std::endl;
  //   find_local_atom_files(path, local_fnames);
  //   size_t num_atoms = local_fnames.size();

  //   // Initialize the atom info file
  //   atom_index_file aif;
  //   aif.atoms.resize(num_atoms);
  //   aif.natoms = num_atoms;
  //   aif.nverts = 0;
  //   aif.nedges = 0;
  //   aif.ncolors = 0;

  //   // hack: only need the type to read the atom header. The body is
  //   // not needed here
  //   typedef atom_file<bool, bool> fake_atom_file_type;
  //   foreach(const std::string& fname, local_fnames) {
  //     std::string absfname = path + "/" + fname;
  //     fake_atom_file_type afile;
  //     afile.input_filename("file", absfname);
  //     afile.load_structure(); // we do not need to load all
  //     procid_t atom_id = afile.atom_id();
  //     atom_file_descriptor& afd(aif.atoms.at(atom_id));
  //     afd.file = absfname;
  //     afd.protocol = "file";
  //     afd.nverts = afile.globalvids().size();
  //     afd.nedges = afile.edge_src_dest().size();      
  //     // Compute adjacency counts for adjacent atoms as well as the
  //     // number of non ghost (local) vertices
  //     std::map<procid_t, size_t> adjatom_count;
  //     // Compute local vertices as well as update the adjacent atom
  //     // count
  //     foreach(const procid_t& other_atom_id, afile.atom())  {
  //       if(other_atom_id == atom_id) ++aif.nverts;
  //       adjatom_count[other_atom_id]++;
  //     }
  //     // Compute max color
  //     foreach(const vertex_color_type& vcolor, afile.vcolor())  {
  // //      std::cout << vcolor << '\t';
  //       aif.ncolors = std::max(aif.ncolors, size_t(vcolor));
  //     }
  //     // count in edges
  //     typedef std::pair<procid_t, procid_t> edge_pair_type;
  //     foreach(const edge_pair_type& edge, afile.edge_src_dest()) {
  //       vertex_id_t target_vid = edge.second;
  //       if(afile.atom().at(target_vid) == atom_id) 
  //         ++aif.nedges;
  //     }
  //     // Fill in afd
  //     typedef std::pair<procid_t, size_t> count_pair;
  //     foreach(const count_pair& pair, adjatom_count) {
  //       afd.adjatoms.push_back(pair.first);
  //       afd.optional_weight_to_adjatoms.push_back(pair.second);
  //     }
  //   } // end of foreach
   
  //   // convert max color to actual number of colors
  //   aif.ncolors++;
  //   // Write results to file
  //   aif.write_to_file(path + "/" + "atom_index.txt");
  // } // and of build atom index file



} // namespace graphlab


#include <graphlab/macros_undef.hpp>
