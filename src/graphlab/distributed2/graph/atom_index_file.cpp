#include <fstream>



#include <boost/filesystem.hpp>

#include <graphlab/logger/assertions.hpp>
#include <graphlab/util/stl_util.hpp>

#include <graphlab/distributed2/graph/atom_file.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>


#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>


#include <graphlab/macros_def.hpp>


using namespace graphlab;
  
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


std::vector<std::vector<size_t> > graphlab::
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


void find_local_atom_files(const std::string& pathname,
                           std::vector<std::string>& files) {
  namespace fs = boost::filesystem;
  fs::path path(pathname);
  assert(fs::exists(path));
  for(fs::directory_iterator iter( path ), end_iter; 
      iter != end_iter; ++iter) {
    if( ! fs::is_directory(iter->status()) ) {
      std::string filename(iter->path().filename());
      size_t period_loc = filename.rfind('.');
      if(period_loc != std::string::npos) {
        std::string ending(filename.substr(period_loc));
        if(ending == ".atom") {
          files.push_back(iter->path().filename());
        }
      }
    }
  }
  std::sort(files.begin(), files.end());

}


struct atom_file_desc_extra {
  atom_file_descriptor desc;
  edge_id_t nlocalverts;
  edge_id_t ninedges;
  atom_file_desc_extra() : 
    nlocalverts(0), ninedges(0) { }
    
  void save(oarchive &oarc) const{
    oarc << desc
         << nlocalverts
         << ninedges;
  }
  
  void load(iarchive &iarc) {
    iarc >> desc
         >> nlocalverts
         >> ninedges;
  }
};




void graphlab::
build_atom_index_file(const std::string& path) {
  // Create a dc comm layer
  std::cout << "Initializing distributed communication layer" << std::endl;
  dc_init_param param;      
  if ( !init_param_from_env(param) ) {
    std::cout << "Failed to get environment variables" << std::endl;
    assert(false);
  }
  distributed_control dc(param);
  dc.services().full_barrier();
  
  // load the vector of filenames
  std::vector<std::string> local_fnames; 
  size_t num_atoms = 0;
  {
    find_local_atom_files(path, local_fnames);
    // compute the fnames that are used by this machine
    std::vector< std::vector< std::string > > partition_fnames;
    dc.services().gather_partition(local_fnames, partition_fnames);
    for(size_t i = 0; i < partition_fnames.size(); ++i) 
      num_atoms += partition_fnames[i].size();
    // update the local fnames      
    local_fnames = partition_fnames[dc.procid()];
    std::cout << "Assigned Names: " << std::endl;
    foreach(std::string fname, local_fnames) 
      std::cout << "(" << fname << ")" << '\t';
    std::cout << std::endl;
  }
  
  dc.services().full_barrier();


  std::cout << "Computing atom files: " << std::endl;
  typedef std::map<procid_t, atom_file_desc_extra> atomid2desc_type;
  atomid2desc_type atomid2desc;
  {
    // hack: only need the type to read the atom payload which is not
    // needed here
    typedef atom_file<bool, bool> fake_atom_file_type;
    foreach(const std::string& fname, local_fnames) {
      std::string absfname = path + "/" + fname;
      fake_atom_file_type afile;
      afile.input_filename("file", absfname);
      afile.load_id_maps();
      procid_t atom_id = afile.atom_id();
      atom_file_desc_extra afd_extra(atomid2desc[atom_id]);
      atom_file_descriptor& afd(afd_extra.desc);
      afd.file = absfname;
      afd.protocol = "file";
      afd.nverts = afile.globalvids().size();
      afd.nedges = afile.edge_src_dest().size();      
      // Compute adjacency counts for adjacent atoms as well as the
      // number of non ghost (local) vertices
      std::map<procid_t, size_t> adjatom_count;
      afd_extra.nlocalverts = 0;
      foreach(const procid_t& other_atom_id, afile.atom())  {
        if(other_atom_id == atom_id) ++afd_extra.nlocalverts;
        adjatom_count[other_atom_id]++;
      }
      // count in edges
      afd_extra.ninedges = 0;
      typedef std::pair<procid_t, procid_t> edge_pair_type;
      foreach(const edge_pair_type& edge, afile.edge_src_dest()) {
        if(edge.second == atom_id) ++afd_extra.ninedges;
      }
      // Fill in afd
      typedef std::pair<procid_t, size_t> count_pair;
      foreach(const count_pair& pair, adjatom_count) {
        afd.adjatoms.push_back(pair.first);
        afd.optional_weight_to_adjatoms.push_back(pair.second);
      }
    } // end of foreach
  } //end of loop over local atoms
  
  dc.services().full_barrier();

  // Gather all the machine counts
  std::vector<atomid2desc_type> all_descriptors;
  all_descriptors[dc.procid()] = atomid2desc;
  dc.services().all_gather(all_descriptors);
  
  // join the maps to build the atom index file
  atom_index_file aif;
  aif.atoms.resize(num_atoms);
  aif.natoms = num_atoms;
  aif.nverts = 0;
  aif.nedges = 0;
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
    }
  }
  // Write results to file
  aif.write_to_file("atom_index.txt");
  

} // and of build atom index file






#include <graphlab/macros_undef.hpp>
