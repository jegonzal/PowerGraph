#ifndef GRAPHLAB_DISTRIBUTED_ATOM_FILE_HPP
#define GRAPHLAB_DISTRIBUTED_ATOM_FILE_HPP
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/util/fs_util.hpp>

//#include <graphlab/distributed2/graph/partitioning/adjacency_list.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {



  /**
   * The contents of an atom file.
   */ 
  template <typename VertexData, typename EdgeData>
  class atom_file {
  public:
    
    atom_file(procid_t atom_id = 0) : 
      atom_id_(atom_id), iarc(NULL), loadstage(0) {  }



    template<typename VFunction, typename EFunction>
    void load_from_txt(const std::string& vdata_fname,
                       const std::string& edata_fname,
                       VFunction vfun, EFunction efun) {
      clear();
      std::map<vertex_id_t, vertex_id_t> global2local;
      { // load all vertex data
        logstream(LOG_INFO) << "Reading " << vdata_fname 
                            << std::endl;
        std::ifstream vdata_fin(vdata_fname.c_str());
        if(!vdata_fin.good()) 
          logstream(LOG_FATAL) << "Invalid vdata filename: " 
                            << vdata_fname << std::endl;
         while(vdata_fin.good()) {           
          std::string line;
          std::getline(vdata_fin, line);
          if(vdata_fin.good()) {
            VertexData vdata_value;
            const vertex_id_t vid = vfun(line, vdata_value);
            vertex_id_t local_vid(globalvids_.size());
            globalvids_.push_back(vid);
            global2local[vid] = local_vid;
            vdata_.push_back(vdata_value);
          }
        }
      } // end of load all the vertex data
      { // Load all the edge data
        logstream(LOG_INFO) << "Reading " << edata_fname 
                            << std::endl;
        std::ifstream edata_fin(edata_fname.c_str());
        if(!edata_fin.good()) 
          logstream(LOG_FATAL) << "Invalid edata filename: " 
                            << edata_fname << std::endl;
        while(edata_fin.good()) {           
          std::string line;
          std::getline(edata_fin, line);
          if(edata_fin.good()) {
            EdgeData edata_value;
            const std::pair<vertex_id_t, vertex_id_t> src_dest_pair =
              efun(line, edata_value);
            const vertex_id_t gsrc( src_dest_pair.first );
            const vertex_id_t gdest( src_dest_pair.second );
            ASSERT_NE(gsrc, gdest);
            typedef std::map<vertex_id_t, vertex_id_t>::const_iterator 
              iterator;
            iterator dest_iter(global2local.find(gdest));
            if(dest_iter == global2local.end()) 
              logstream(LOG_FATAL) << "Invalid destination vid "
                                << gdest
                                << " in file " << edata_fname 
                                << std::endl;
            const vertex_id_t local_dest_vid(dest_iter->second);
            
            iterator src_iter(global2local.find(gsrc));
            vertex_id_t local_src_vid(-1);
            // if the source was not found then it is a ghost and we
            // will need to update the globalvid list
            if(src_iter == global2local.end()) {
              local_src_vid = globalvids_.size();
              globalvids_.push_back(gsrc);
              // update the map
              global2local[gsrc] = local_src_vid;              
            } else local_src_vid = src_iter->second;
            edata_.push_back(edata_value);
            ASSERT_NE(local_src_vid, local_dest_vid);
            edge_src_dest_.push_back(std::make_pair(local_src_vid, 
                                                    local_dest_vid));            
          } // end of if statement
        } // end of while loop
        edata_fin.close();
      } // end of load all the edge data
      ASSERT_EQ(global2local.size(), globalvids_.size());
      // Cleanup some local data structures
      vcolor_.resize(globalvids_.size());
      atom_.resize(globalvids_.size());      
    } // end of construct atom file from a text file



    /**
     * Associates this object with a particular atom file
     * and a particular protocol. The contents of the file (through
     * the accessors) are not yet loaded until the load functions are 
     * called. \see load_id_maps() \see load_structure() \see load_all()
     */
    void input_filename(std::string protocol, std::string filename) {
      if (protocol == "file") {
        fin.open(filename.c_str());
        iarc = new iarchive(fin);
        loadstage = 0;
      } else {
        logger(LOG_FATAL, "Unrecognized protocol: %s", protocol.c_str());
      }
    }
    /**
     * Only load the globalvids and the globaleids from the input.
     * input_filename() must be called first.
     * All remaining entries are not loaded.
     */
    void load_id_maps() {
      if (loadstage == 0) {
        (*iarc) >> atom_id_ >> globalvids_ >> globaleids_ 
                >> atom_ >> vcolor_;
        ++loadstage;
      }
    }
  
    /**
     * Loads all entries except the vdata and edata entries.
     * input_filename() must be called first.
     */
    void load_structure() {
      if (loadstage < 1) load_id_maps();
      if (loadstage == 1) {
        (*iarc) >> edge_src_dest_;
        ++loadstage;
      }
    }
  
    /**
     * Completely loads the file defined by the input_filename()
     * All datastructures are now accessible.
     */
    void load_all() {
      if (loadstage < 2) load_structure();
      if (loadstage == 2) {
        (*iarc) >> vdata_ >> edata_;
        ++loadstage;
      }

    }
  
    void write_to_file(std::string protocol, std::string outfilename) {
      assert(protocol == "file");
      std::ofstream fout;
      fout.open(outfilename.c_str());
      assert(fout.good());
      oarchive oarc(fout);
      oarc << atom_id_ << globalvids_ << globaleids_ << atom_ << vcolor_ 
           << edge_src_dest_
           << vdata_ << edata_;
    }
  
    /**
     * Clears all loaded data and closes the file.
     */
    void clear() {
      if (iarc != NULL) {
        delete iarc;
        fin.close();
        iarc = NULL;
      }
      globalvids_.clear();
      globaleids_.clear();
      atom_.clear();
      vcolor_.clear();
      edge_src_dest_.clear();
      vdata_.clear();
      edata_.clear();
    }
  
    ~atom_file() {
      clear();
    }
    
    inline const procid_t& atom_id() const { return atom_id_; }
    inline const std::string& protocol() const { return protocol_; }
    inline const std::string& filename() const { return filename_; }
    inline const std::vector<vertex_id_t>& globalvids() const { return globalvids_; }
    inline const std::vector<edge_id_t>& globaleids() const { return globaleids_; }
    inline const std::vector<procid_t>& atom() const { return atom_; }
    inline const std::vector<vertex_color_type>& vcolor() const { return vcolor_; }
    inline const std::vector< std::pair<vertex_id_t, vertex_id_t> >& edge_src_dest() const { return edge_src_dest_; }
    inline const std::vector<VertexData>& vdata() const { return vdata_; }
    inline const std::vector<EdgeData>& edata() const { return edata_; }

    inline procid_t& atom_id() { return atom_id_; }
    inline std::string& protocol() { return protocol_; }
    inline std::string& filename() { return filename_; }
    inline std::vector<vertex_id_t>& globalvids() { return globalvids_; }
    inline std::vector<edge_id_t>& globaleids() { return globaleids_; }
    inline std::vector<procid_t>& atom() { return atom_; }
    inline std::vector<vertex_color_type>& vcolor() { return vcolor_; }
    inline std::vector< std::pair<vertex_id_t, vertex_id_t> >& edge_src_dest() { return edge_src_dest_; }
    inline std::vector<VertexData>& vdata() { return vdata_; }
    inline std::vector<EdgeData>& edata() { return edata_; }


    void costly_checking() {
      ASSERT_EQ(globalvids().size(), atom().size());
      ASSERT_EQ(globalvids().size(), vcolor().size());
      std::set<vertex_id_t> allgvids;
      foreach(const vertex_id_t& gvid, globalvids()) {
        ASSERT_EQ(allgvids.count(gvid), 0);
        allgvids.insert(gvid);
      }
      typedef std::pair<vertex_id_t, vertex_id_t> pair_type;
      foreach(const pair_type pair, edge_src_dest()) {
        const vertex_id_t source_lvid = pair.first;
        const vertex_id_t target_lvid = pair.second;
        ASSERT_NE(source_lvid, target_lvid);
        ASSERT_LT(source_lvid, globalvids().size());
        ASSERT_LT(target_lvid, globalvids().size());
        ASSERT_TRUE(atom()[source_lvid] == atom_id() ||
                    atom()[target_lvid] == atom_id());        
      }

    }

  private:
    procid_t atom_id_;


    std::string protocol_;
    std::string filename_;
  
    // for file protocol
    std::ifstream fin;
    iarchive *iarc;
    size_t loadstage;
  

    std::vector<vertex_id_t> globalvids_;
    std::vector<edge_id_t> globaleids_;
    std::vector<procid_t> atom_;
    std::vector<vertex_color_type> vcolor_;
    std::vector< std::pair<vertex_id_t, vertex_id_t> > edge_src_dest_;
    std::vector<VertexData> vdata_;
    std::vector<EdgeData> edata_;
  };

}

#include <graphlab/distributed2/graph/atom_index_file.hpp>

namespace graphlab {
  /**
     Converts the partition partid as an atom.
     graph must be colored!
  */
  template <typename VertexData, typename EdgeData>
  void graph_partition_to_atom(const graph<VertexData, EdgeData> &g, 
                               const std::vector<uint32_t>& vertex2part,
                               uint32_t partid,
                               atom_file<VertexData, EdgeData> &atom,
                               bool noglobaleids = false) {
    atom.clear();
    // build the ID mappings
    // collect the set of vertices / edges that are within this atom
    // that would be all vertices with this partition ID as well as all neighbors
    boost::unordered_set<vertex_id_t> struct_vertices;
    boost::unordered_set<edge_id_t>  struct_edges;
    for (size_t v = 0;v < g.num_vertices(); ++v) {
      if (vertex2part[v] == partid) {
        // loop through all the edges
        struct_vertices.insert(v);
        foreach(edge_id_t eid, g.out_edge_ids(v)) {
          struct_edges.insert(eid);
          struct_vertices.insert(g.target(eid));
        }
        foreach(edge_id_t eid, g.in_edge_ids(v)) {
          struct_edges.insert(eid);
          struct_vertices.insert(g.source(eid));
        }
  
      }
    }

    boost::unordered_map<vertex_id_t, vertex_id_t> global2localvid;
    {
      // done. now construct the mappings
      foreach(vertex_id_t vid, struct_vertices) {
        global2localvid[vid] = atom.globalvids().size();
        atom.globalvids().push_back(vid);
        atom.atom().push_back(vertex2part[vid]);
        atom.vcolor().push_back(g.color(vid));
        //atom.vdata().push_back(g.vertex_data(vid));
        if (vertex2part[vid] == partid) atom.vdata().push_back(g.vertex_data(vid));
      }
    }
    {
      foreach(edge_id_t eid, struct_edges) {
        atom.globaleids().push_back(eid);
        atom.edge_src_dest().push_back(std::make_pair<vertex_id_t,
                                       vertex_id_t>(global2localvid[g.source(eid)],
                                                    global2localvid[g.target(eid)]));
        //atom.edata().push_back(g.edge_data(eid));
        
        if (vertex2part[g.target(eid)] == partid) atom.edata().push_back(g.edge_data(eid));
      }
    }
    if (noglobaleids) atom.globaleids().clear();
  }



  /**
     Converts a graph into an on disk atom representation.
     graph must be colored before calling this function.
     The index file is written to idxfilename while the atoms 
     are written to atombasename.0, atombasename.1, etc
  */
  template <typename VertexData, typename EdgeData>
  void graph_partition_to_atomindex(graph<VertexData, EdgeData> &graph,
                                    const std::vector<uint32_t>& vertex2part,
                                    std::string idxfilename,
                                    std::string atombasename,
                                    bool noglobaleids = false) {
    assert(graph.valid_coloring());
    // get the number of partitions
    uint32_t numparts = 0 ;
    for (size_t i = 0;i < vertex2part.size(); ++i) {
      numparts = std::max(numparts, vertex2part[i]);
    }
    ++numparts;
    atom_index_file idxfile;
    std::ofstream fout(idxfilename.c_str());
    idxfile.nverts = graph.num_vertices();
    idxfile.nedges = graph.num_edges();
    idxfile.natoms = numparts;
    idxfile.ncolors = 0;
    std::vector<std::vector<size_t> > edgeweights(numparts);
    for (size_t i = 0;i < edgeweights.size(); ++i) edgeweights[i].resize(numparts, 0);
    
    for (size_t i = 0;i < graph.num_vertices(); ++i) {
      idxfile.ncolors = std::max<size_t>(idxfile.ncolors, graph.color(i));
      foreach(edge_id_t e, graph.in_edge_ids(i)) {
        vertex_id_t srcatom = vertex2part[graph.source(e)];
        vertex_id_t destatom = vertex2part[graph.target(e)];
        edgeweights[srcatom][destatom]++;
        edgeweights[destatom][srcatom]++;
      }
    }
    
    // scale down the max edge weight to about 1000
    size_t maxedgeweight = 0;
    for (size_t i = 0; i < numparts; ++i) {
      for (size_t j = 0;j < numparts; ++j) {
        maxedgeweight = std::max(maxedgeweight, edgeweights[i][j]);
      }
    }
    
    for (size_t i = 0; i < numparts; ++i) {
      for (size_t j = 0;j < numparts; ++j) {
        edgeweights[i][j] = (size_t)(1000.0 * double(edgeweights[i][j]) / maxedgeweight) + 1;
      }
    }
    idxfile.ncolors++;
    for (size_t i = 0; i < numparts; ++i) {
      std::string atomfilename = atombasename + "." + tostr(i);
      atom_file<VertexData, EdgeData> atomfile;
      graph_partition_to_atom(graph, vertex2part, i, atomfile, noglobaleids);
      atomfile.write_to_file("file", atomfilename);
      
      atom_file_descriptor desc;
      desc.protocol = "file";
      desc.file = atomfilename;
      // get list of adjacent atoms
      std::set<size_t> adjatoms;
      for (size_t v = 0; v < atomfile.atom().size(); ++v) {
        if (atomfile.atom()[v] != i) adjatoms.insert(atomfile.atom()[v]);
      }
      foreach(size_t v, adjatoms) {
        desc.adjatoms.push_back(v);
        desc.optional_weight_to_adjatoms.push_back(edgeweights[i][v]);
      }
      desc.nverts = atomfile.globalvids().size();
      desc.nedges = atomfile.globaleids().size();
      idxfile.atoms.push_back(desc);
    }
    idxfile.write_to_file(idxfilename);
  }


  


  
  template<typename VertexData, typename EdgeData,
           typename VFunction, typename EFunction>
  void generate_atom_files(const std::string& path,
                           const std::string& atompath,
                           VFunction vfun, EFunction efun) {
    typedef atom_file<VertexData, EdgeData> atom_file_type;
    const std::string vdata_suffix(".vdata_txt");
    const std::string edata_suffix(".edata_txt");
    dc_init_param param;         
    if( ! init_param_from_mpi(param) ) {
      logstream(LOG_FATAL) 
        << "Failed MPI laucher!" << std::endl;
    }
    distributed_control dc(param);
    dc.full_barrier();      
    logstream(LOG_INFO) 
      << "Initializing distributed shuffler object." 
      << std::endl;    

    // Determine the local vdata files on this machine
    std::vector<std::string> local_fnames;
    fs_util::list_files_with_suffix(path, vdata_suffix,
                                    local_fnames);
    std::vector< std::vector<std::string> > partition_fnames;
    dc.gather_partition(local_fnames, partition_fnames);
    local_fnames = partition_fnames.at(dc.procid());   
    std::vector< procid_t > local_atomids(local_fnames.size(), 
                                          procid_t(-1));

    std::vector<procid_t> atomid2proc;
    { // comute local atom ids
      procid_t starting_atomid(0);
      for(size_t i = 0; i < dc.procid(); ++i)
        starting_atomid += partition_fnames.at(i).size();
      for(size_t i = 0; i < local_fnames.size(); ++i)
        local_atomids.at(i) = starting_atomid + i;      
      for(size_t i = 0; i < partition_fnames.size(); ++i) {
        for(size_t j = 0; j < partition_fnames[i].size(); ++j) {
          atomid2proc.push_back(i);
        }
      }
    }

    std::vector< std::string > local_atom_fnames(local_fnames.size());
    { // compute atom file names
      for(size_t i = 0; i < local_atom_fnames.size(); ++i) 
        local_atom_fnames.at(i) =  
          fs_util::change_suffix(local_fnames[i], ".atom");        
    }

    typedef std::pair<vertex_id_t, vertex_id_t> edge_t;
    typedef std::pair<procid_t, edge_t> atom_edge_t;
    std::vector< std::vector<atom_edge_t> > proc2out_edges(dc.numprocs());

    typedef std::map<procid_t, std::vector<vertex_id_t> > 
      atomid2vertexids_type;
    std::vector< atomid2vertexids_type > gather_vec(dc.numprocs());    
    atomid2vertexids_type& atomid2vertexids(gather_vec.at(dc.procid()));
    // Build out each atom file
    for(size_t i = 0; i < local_fnames.size(); ++i) {
      const std::string vdata_fname = path + "/" +
        fs_util::change_suffix(local_fnames[i], vdata_suffix);
      const std::string edata_fname = path + "/" +
        fs_util::change_suffix(local_fnames[i], edata_suffix);
      atom_file_type afile;
      afile.load_from_txt(vdata_fname, edata_fname, vfun, efun);
      afile.atom_id() = local_atomids[i];
      { // update the atomid2vertexids map
        std::vector<vertex_id_t>& 
          localverts(atomid2vertexids[afile.atom_id()]);
        localverts.resize(afile.globalvids().size());
        for(size_t j = 0; j < localverts.size(); ++j)
          localverts[j] = afile.globalvids().at(j);
        foreach(const edge_t& edge, afile.edge_src_dest()) {
          const procid_t atomid(afile.atom()[edge.first]);
          const vertex_id_t source(afile.globalvids()[edge.first]);
          const vertex_id_t target(afile.globalvids()[edge.second]);
          const edge_t global_edge(std::make_pair(source, target));
          if(atomid != afile.atom_id()) 
            proc2out_edges[atomid2proc[atomid]].
              push_back(std::make_pair(atomid, global_edge));
        }

        
      }
      afile.filename() = local_atom_fnames.at(i);
      afile.protocol() = "file";
      // construct the atom filename
      const std::string fname(path + "/" + afile.filename() +"_tmp");
      logstream(LOG_INFO) << "Saving atom file: " 
                          << fname
                          << std::endl;
      afile.write_to_file("file", fname );
    }

    std::vector< std::vector<edge_t> > atom2out_edges(atomid2proc.size());
    {
      std::vector< std::vector<atom_edge_t> > result;
      mpi_tools::all2all(proc2out_edges, result);
      for(size_t i = 0; i < result.size(); ++i) {
        for(size_t j = 0; j < result[i].size(); ++j) {
          atom2out_edges[result[i][j].first].push_back(result[i][j].second);
        }
      }
    }


    const size_t ROOT_NODE(0);
    dc.gather(gather_vec, ROOT_NODE);
    std::vector<procid_t> vid2atomid;
    if(dc.procid() == 0) {
      logstream(LOG_INFO) << "Computing vertex to atom map."  << std::endl;
      // compute max vertex id
      vertex_id_t max_vid(0);
      typedef atomid2vertexids_type::value_type pair_type;
      for(size_t i = 0; i < gather_vec.size(); ++i) {         
        foreach(const pair_type& pair, gather_vec[i]) {
          for(size_t j = 0; j < pair.second.size(); ++j) {
            max_vid = std::max(max_vid, pair.second[j]);
          }
        }
      }
      const vertex_id_t nverts(max_vid + 1);
      logstream(LOG_INFO) << "Num vertices " << nverts << std::endl;
      // Construct inverse map
      vid2atomid.resize(nverts, procid_t(-1));
      // Invert the map
      for(size_t i = 0; i < gather_vec.size(); ++i) {
        foreach(const pair_type& pair, gather_vec[i]) {
          const procid_t atomid(pair.first);
          for(size_t j = 0; j < pair.second.size(); ++j) {
            const vertex_id_t vid(pair.second[j]);
            ASSERT_LT(vid, nverts);
            ASSERT_EQ(vid2atomid[vid], vertex_id_t(-1));
            vid2atomid[vid] = atomid;
          }
        }
      }
      // check all the entries
      for(size_t i = 0; i < vid2atomid.size(); ++i)
        ASSERT_NE(vid2atomid[i], procid_t(-1));
      // broadcast the map to all machines
      dc.broadcast(vid2atomid, true);
    } else {
      dc.broadcast(vid2atomid, false);
    }

    // Loop through local atom files
    logstream(LOG_INFO) << "Updating local atoms."  << std::endl;
    



    // Build out each atom file
    for(size_t i = 0; i < local_atom_fnames.size(); ++i) {
      const std::string tmp_fname(path + "/" + local_atom_fnames[i] + "_tmp");
      logstream(LOG_INFO) << "Reloading " << tmp_fname << std::endl;
      atom_file_type afile;
      afile.input_filename("file", tmp_fname);
      afile.load_all();
      std::map<vertex_id_t, vertex_id_t> global2local;
      // Update the atom location information
      for(size_t j = 0; j < afile.globalvids().size(); ++j) {
        const vertex_id_t gvid(afile.globalvids()[j]);
        global2local[gvid] = j;
        ASSERT_LT(gvid, vid2atomid.size());
        const procid_t atomid(vid2atomid[gvid]);
        ASSERT_NE(atomid, procid_t(-1));
        afile.atom()[j] = atomid;
      }
      // Adde the extra edges
      for(size_t j = 0; j < atom2out_edges[i].size(); ++j) {
        const edge_t edge(atom2out_edges[i][j]);
        ASSERT_TRUE(global2local.find(edge.first) != global2local.end());
        const vertex_id_t source_lvid(global2local[edge.first]);        
        if(global2local.find(edge.second) == global2local.end()) {
          const vertex_id_t target_lvid = afile.globalvids().size();
          global2local[edge.second] = target_lvid;
          afile.globalvids().push_back(edge.second);          
        }
        const vertex_id_t target_lvid(global2local[edge.second]);
        afile.edge_src_dest().push_back
          (std::make_pair(source_lvid, target_lvid));

      }

      const std::string fname( atompath + "/" + local_atom_fnames[i]);
      // Resave the atom file
      logstream(LOG_INFO) << "Final save of atom file: " 
                          << fname
                          << std::endl;
      afile.write_to_file("file", fname );
    }

    

  } // end of generate atom file










} // end namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif
