#ifndef GRAPHLAB_DISTRIBUTED_ATOM_FILE_HPP
#define GRAPHLAB_DISTRIBUTED_ATOM_FILE_HPP
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <boost/unordered_map.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_env.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/distributed2/graph/graph_fragment.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * The contents of an atom file.
   */ 
  template <typename VertexData, typename EdgeData>
  class atom_file {
  public:
    atom_file(): iarc(NULL), loadstage(0) {  }

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
        (*iarc) >> globalvids_ >> globaleids_ >> atom_;
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
        (*iarc) >> vcolor_ >> edge_src_dest_;
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
      oarc << globalvids_ << globaleids_ << atom_ << vcolor_ << edge_src_dest_
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
  
  
    inline const std::string& protocol() const { return protocol_; }
    inline const std::string& filename() const { return filename_; }
    inline const std::vector<vertex_id_t>& globalvids() const { return globalvids_; }
    inline const std::vector<edge_id_t>& globaleids() const { return globaleids_; }
    inline const std::vector<procid_t>& atom() const { return atom_; }
    inline const std::vector<vertex_color_type>& vcolor() const { return vcolor_; }
    inline const std::vector< std::pair<vertex_id_t, vertex_id_t> >& edge_src_dest() const { return edge_src_dest_; }
    inline const std::vector<VertexData>& vdata() const { return vdata_; }
    inline const std::vector<EdgeData>& edata() const { return edata_; }

    inline std::string& protocol() { return protocol_; }
    inline std::string& filename() { return filename_; }
    inline std::vector<vertex_id_t>& globalvids() { return globalvids_; }
    inline std::vector<edge_id_t>& globaleids() { return globaleids_; }
    inline std::vector<procid_t>& atom() { return atom_; }
    inline std::vector<vertex_color_type>& vcolor() { return vcolor_; }
    inline std::vector< std::pair<vertex_id_t, vertex_id_t> >& edge_src_dest() { return edge_src_dest_; }
    inline std::vector<VertexData>& vdata() { return vdata_; }
    inline std::vector<EdgeData>& edata() { return edata_; }


  private:
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
    dense_bitset goodvertices(g.num_vertices()), goodedges(g.num_edges());

    for (size_t v = 0;v < g.num_vertices(); ++v) {
      if (vertex2part[v] == partid) {
        // add myself and all neighbors
        goodvertices.set_bit_unsync(v);
        // loop through all the edges
        foreach(edge_id_t eid, g.out_edge_ids(v)) {
          goodedges.set_bit_unsync(eid);
          goodvertices.set_bit_unsync(g.target(eid));
        }
        foreach(edge_id_t eid, g.in_edge_ids(v)) {
          goodedges.set_bit_unsync(eid);
          goodvertices.set_bit_unsync(g.source(eid));
        }
  
      }
    }

    boost::unordered_map<vertex_id_t, vertex_id_t> global2localvid;
    {
      // done. now construct the mappings
      uint32_t vid;
      if (goodvertices.first_bit(vid)) {
        do {
          global2localvid[vid] = atom.globalvids().size();
          atom.globalvids().push_back(vid);
          atom.atom().push_back(vertex2part[vid]);
          atom.vcolor().push_back(g.color(vid));
          atom.vdata().push_back(g.vertex_data(vid));
        } while(goodvertices.next_bit(vid));
      }
    }
    {
      uint32_t eid;
      if (goodedges.first_bit(eid)) {
        do {
          atom.globaleids().push_back(eid);
          atom.edge_src_dest().push_back(std::make_pair<vertex_id_t,
                                         vertex_id_t>(global2localvid[g.source(eid)],
                                                      global2localvid[g.target(eid)]));
          atom.edata().push_back(g.edge_data(eid));
        } while(goodedges.next_bit(eid));
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
  void graph_partition_to_atomindex(const graph<VertexData, EdgeData> &graph,
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
    std::ofstream fout(idxfilename.c_str());
    fout << graph.num_vertices() << "\t" << graph.num_edges() << "\t" << numparts << "\n";
  
    for (size_t i = 0; i < numparts; ++i) {
      std::string atomfilename = atombasename + "." + tostr(i);
      atom_file<VertexData, EdgeData> atomfile;
      graph_partition_to_atom(graph, vertex2part, i, atomfile, noglobaleids);
      atomfile.write_to_file("file", atomfilename);
    
      // get list of adjacent atoms
      std::set<size_t> adjatoms;
      for (size_t v = 0; v < atomfile.atom().size(); ++v) {
        if (atomfile.atom()[v] != i) adjatoms.insert(atomfile.atom()[v]);
      }
    
      fout << atomfile.globalvids().size() << "\t" << atomfile.globaleids().size() 
           << "\t" << adjatoms.size()  << "\t";
      foreach(size_t v, adjatoms) {
        fout << v << " ";
      }
      fout << "\t"<<"file://" << atomfilename << "\n";
    }
  }





  namespace atom_file_impl {




    template <typename VertexData, typename EdgeData>
    class atom_shuffler {
    public:
      typedef VertexData vertex_data_type;
      typedef EdgeData edge_data_type;
      typedef atom_shuffler<VertexData, EdgeData> atom_shuffler_type;
      typedef atom_file<VertexData, EdgeData> atom_file_type;

      struct atom_info {
        graphlab::mutex mutex;
        atom_file_type atom_file;
        std::map<vertex_id_t, vertex_id_t> global2local;
        std::string vdatafn;
        std::string edatafn;
        std::ofstream vdatastream;
        std::ofstream edatastream;        
        vertex_id_t get_local_vid(vertex_id_t gvid) {
          std::map<vertex_id_t, vertex_id_t>::iterator iter 
            = global2local.find(gvid);
          assert(iter != global2local.end());
          return iter->second;
        }
      };

      typedef std::vector< atom_info* > atom_map_type;

    private:

      procid_t owning_machine(procid_t atomid) {
        return atomid % rmi.numprocs();
      }
      
      bool is_local(procid_t atomid) {
        return owning_machine(atomid) == rmi.procid();
      }
      
      atom_info& get_atom_info(procid_t atomid, bool create = false) {
        // Check that this is the owning machine and the file has been allocated
        assert(atomid < atomid2info.size());
        if(create) {
          assert(atomid2info[atomid] == NULL);
          atomid2info[atomid] = new atom_info();
        }
        assert(atomid2info[atomid] != NULL);
        return *atomid2info[atomid];
      }


      dc_dist_object< atom_shuffler_type > rmi;
      size_t num_atoms;
      std::string path;
      atom_map_type atomid2info;

    public:
      atom_shuffler(distributed_control& dc,
                    const size_t& num_atoms,
                    const std::string& path) : 
        rmi(dc, this), num_atoms(num_atoms), path(path),
        atomid2info(num_atoms, NULL) { 
        // Initialize the atom files for this machine
        for(procid_t atomid = 0; atomid < num_atoms; ++atomid) {         
          // if this is the owning machine setup the map and files
          if(is_local(atomid)) {
            // get the atom info (and create it)
            atom_info& ainfo(get_atom_info(atomid, true));
            ainfo.mutex.lock();
            { // open vdata temporary storage file
              std::stringstream strm;
              strm << path << "/"
                   << "tmp_vdata_"
                   << std::setw(3) << std::setfill('0')
                   << atomid
                   << ".bin";
              ainfo.vdatafn = strm.str();
              ainfo.vdatastream.open(ainfo.vdatafn.c_str(), 
                                     std::ios::binary | 
                                     std::ios::out |
                                     std::ios::trunc );
              assert(ainfo.vdatastream.good());
            }
            { // open edata temporary storage file
              std::stringstream strm;
              strm << path << "/"
                   << "tmp_edata_"
                   << std::setw(3) << std::setfill('0')
                   << atomid
                   << ".bin";
              ainfo.edatafn = strm.str();
              ainfo.edatastream.open(ainfo.edatafn.c_str(), 
                                       std::ios::binary | 
                                       std::ios::out |
                                       std::ios::trunc );
              assert(ainfo.edatastream.good());
            }
            ainfo.mutex.unlock();
          } //end of if owning machine
        } //end of for loop

        dc.services().barrier();
      } //end of constructor
      


      ~atom_shuffler() {
        // Clear the atom info map
        for(size_t i = 0; i < atomid2info.size(); ++i) {
          if(atomid2info[i] != NULL) {
            delete atomid2info[i];
            atomid2info[i] = NULL;
          }
        }
      } // end of destructor
      

      void add_vertex_local(const procid_t& to_atomid,
                            const vertex_id_t& gvid,
                            const procid_t& atomid,
                            const vertex_color_type& vcolor,
                            const vertex_data_type& vdata) {
        assert(is_local(to_atomid));
        // Get the atom info
        atom_info& ainfo(get_atom_info(to_atomid));
        ainfo.mutex.lock();
        // test to see if the vertex has been added already
        if(ainfo.global2local.find(gvid) == ainfo.global2local.end()) { // first add
          vertex_id_t localvid = ainfo.atom_file.globalvids().size();
          ainfo.atom_file.globalvids().push_back(gvid);
          ainfo.atom_file.vcolor().push_back(vcolor);
          ainfo.atom_file.atom().push_back(atomid);
          ainfo.global2local[gvid] = localvid;
          assert(ainfo.vdatastream.good());
          oarchive oarc(ainfo.vdatastream);
          oarc << vdata;          
        }
        ainfo.mutex.unlock();
      }

      void add_vertex(const procid_t& to_atomid,
                      const vertex_id_t& gvid,
                      const procid_t& atomid,
                      const vertex_color_type& vcolor,
                      const vertex_data_type& vdata) {

        if(is_local(to_atomid)) {
          // local add
          add_vertex_local(to_atomid, gvid, atomid, vcolor, vdata);
        } else {
          // remote add
          rmi.remote_request(owning_machine(to_atomid),
                             &atom_shuffler_type::add_vertex_local,
                             to_atomid, gvid, atomid, vcolor, vdata);
        }
      } // end of add_vertex



      void add_edge_local(const procid_t& to_atomid,
                          const vertex_id_t& source_gvid,
                          const vertex_id_t& target_gvid,
                          const edge_data_type& edata) {
        assert(is_local(to_atomid));
        // Get the atom info
        atom_info& ainfo(get_atom_info(to_atomid));
        ainfo.mutex.lock();
        vertex_id_t source_lvid(ainfo.get_local_vid(source_gvid));
        vertex_id_t target_lvid(ainfo.get_local_vid(target_gvid));        
        ainfo.atom_file.edge_src_dest().push_back(std::make_pair(source_lvid, 
                                                                 target_lvid));
        assert(ainfo.edatastream.good());
        oarchive oarc(ainfo.edatastream);
        oarc << edata;
        ainfo.mutex.unlock();
      } // end of add edge local


      void add_edge(const procid_t& to_atomid,
                    const vertex_id_t& source_gvid,
                    const vertex_id_t& target_gvid,
                    const edge_data_type& edata) {
        if(is_local(to_atomid)) {
          add_edge_local(to_atomid, source_gvid, target_gvid, edata);
        } else {
          rmi.remote_request(owning_machine(to_atomid),
                             &atom_shuffler_type::add_edge_local,
                             to_atomid, source_gvid, target_gvid, edata);
        }
      } // end of add edge
    }; // end of atom shuffler
  }; // end of namespace atom_file_impl;






  template <typename VertexData, typename EdgeData>
  void graph_partition_to_atoms(const std::string& path) {
    typedef VertexData vertex_data_type;
    typedef EdgeData edge_data_type;
    std::cout << "Building atom shuffler object. " << std::endl;
    typedef atom_file_impl::atom_shuffler<VertexData, EdgeData>
      atom_shuffler_type;
    // startup distributed control
    dc_init_param param;
    if ( !init_param_from_env(param) ) {
      std::cout << "Failed to get environment variables" << std::endl;
      assert(false);
    }
    distributed_control dc(param);
    dc.services().barrier();
    
    

    // load the maps
    vertex_color_type max_color = 0;
    std::vector<vertex_color_type> vertex2color;
    {
      std::cout << "Loading coloring file" << std::endl;
      std::string absfname = path + "/coloring.txt";
      std::ifstream fin(absfname.c_str());
      while(fin.good()) {
        vertex_color_type vcolor;
        fin >> vcolor;
        max_color = std::max(max_color, vcolor);
        if(fin.good()) vertex2color.push_back(vcolor);
      }
      fin.close();
    }
    //size_t num_colors( max_color + 1 );

    procid_t max_atomid = 0;
    std::vector<procid_t> vertex2atomid;
    {
      std::cout << "Loading paritioning file" << std::endl;
      std::string absfname = path + "/partitioning.txt";
      std::ifstream fin(absfname.c_str());
      while(fin.good()) {
        procid_t atomid;
        fin >> atomid;
        max_atomid = std::max(max_atomid, atomid);
        if(fin.good()) vertex2atomid.push_back(atomid);
      }
      fin.close();
    }
    size_t num_atoms( max_atomid + 1 );

    // get the filenames accessible by this process
    std::vector<std::string> fnames;
    {
      graph_fragment::list_structure_files(path, fnames);
      // compute the fnames that are used by this machine
      std::vector< std::vector< std::string > > partition_fnames;
      dc.services().gather_partition(fnames, partition_fnames);

      // update the local fnames      
      fnames = partition_fnames[dc.procid()];
      std::cout << "Assigned Names: " << std::endl;
      foreach(std::string fname, fnames) 
        std::cout << "(" << fname << ")" << '\t';
      std::cout << std::endl;
    }

    // Create the atom shuffler
    atom_shuffler_type atom_shuffler(dc, num_atoms, path);
    std::cout << "Waiting at barrier " << std::endl;
    dc.services().comm_barrier();
    dc.services().barrier();


    std::cout << "Adding all vertices: " << std::endl;
    foreach(const std::string fname, fnames) {
      graph_fragment::file_description desc;
      graph_fragment::structure_description structure;
      { // Get the local structure description
        std::string absfname = path + "/" + fname;
        std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
        assert(fin.good());
        graphlab::iarchive iarc(fin);
        iarc >> desc >> structure;
        fin.close();
      }
      { // read all the vdata
        std::string absfname = path + "/" + desc.vertex_data_filename;
        std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
        assert(fin.good());
        graphlab::iarchive iarc(fin);        
        // Loop through all the structures reading in the vertex data
        for(vertex_id_t i = 0, vid = desc.begin_vertex; 
            vid < desc.end_vertex; ++vid, ++i) {
          assert(vid < vertex2atomid.size());
          assert(vid < vertex2color.size());
          vertex_data_type vdata;
          assert(fin.good());
          iarc >> vdata;
          assert(fin.good());
          // send the vertex data to the owning atoms
          atom_shuffler.add_vertex(vertex2atomid[vid],
                                   vid,
                                   vertex2atomid[vid],
                                   vertex2color[vid],
                                   vdata);
          // Loop over the vertex neighbors and send to neighbor atoms
          foreach(vertex_id_t neighbor_vid, structure.neighbor_ids[i]) {
            assert(neighbor_vid < vertex2atomid.size());
            // if the neighbor is stored in a different atom file then
            // send the vertex data to the neighbor for ghosting purposes
            if(vertex2atomid[neighbor_vid] != vertex2atomid[vid]) {
              // send the vertex data to the owning atoms
              atom_shuffler.add_vertex(vertex2atomid[neighbor_vid],
                                       vid,
                                       vertex2atomid[vid],
                                       vertex2color[vid],
                                       vdata);
            } // end of if neighbor is in different atom
          }
        }
        fin.close();
      } // end of read all vdata
    } // end of adding all vertices


    std::cout << "Waiting at barrier " << std::endl;
    dc.services().comm_barrier();
    dc.services().barrier();
    
    std::cout << "Adding all edges" << std::endl;
    
    foreach(const std::string fname, fnames) {
      graph_fragment::file_description desc;
      graph_fragment::structure_description structure;
      
      { // Get the local structure description
        std::string absfname = path + "/" + fname;
        std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
        assert(fin.good());
        graphlab::iarchive iarc(fin);
        iarc >> desc >> structure;
        fin.close();
      }


      { // read all the edge data
        std::string absfname = path + "/" + desc.edge_data_filename;
        std::ifstream fin(absfname.c_str(), std::ios::binary | std::ios::in);
        assert(fin.good());
        graphlab::iarchive iarc(fin);        
        // Loop through all the structures reading in the vertex data
        for(vertex_id_t i = 0, target = desc.begin_vertex; 
            target < desc.end_vertex; ++target, ++i) {
          assert(target < vertex2atomid.size());
          assert(i < structure.in_neighbor_ids.size());
          foreach(vertex_id_t source, structure.in_neighbor_ids[i]) {
            assert(source < vertex2atomid.size());
            edge_data_type edata;
            assert(fin.good());
            iarc >> edata;
            assert(fin.good());
            atom_shuffler.add_edge(vertex2atomid[source],
                                   source, target, edata);            
            if(vertex2atomid[source] != vertex2atomid[target])
              atom_shuffler.add_edge(vertex2atomid[target],
                                     source, target, edata);                       
          } // end of loop over out neighbors
        } // end of loop over vertices in file
        fin.close();
      } // end of read all vdata
    } // end of for loop

    std::cout << "Waiting at barrier " << std::endl;
    dc.services().comm_barrier();
    dc.services().barrier();
    

  } // end of graph partition to atom index



} // end namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif
