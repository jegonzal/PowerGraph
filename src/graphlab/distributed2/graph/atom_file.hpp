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
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/dense_bitset.hpp>
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
    private:

      dc_dist_object< atom_shuffler_type > rmi;
      std::map<procid_t, atom_file_type> atom2atomfile;
      std::string path;
      std::map<procid_t, std::string> atom2vdatafn;
      std::map<procid_t, std::string> atom2edatafn;
      std::map<procid_t, std::ofstream> atom2vdatastream;
      std::map<procid_t, std::ofstream> atom2edatastream;
    public:
      atom_shuffler(distributed_control& dc,
                    const std::string& path) : 
        rmi(dc, this), path(path) { 
        // Initialize the atom files for this machine
        for(procid_t i = 0; i < rmi.numprocs(); ++i) {
          procid_t owningmachine = i % rmi.numprocs();
          // if this is the owning machine setup the map and files
          if(owningmachine == rmi.procid()) {
            atom2atomfile[i] = atom_file_type();
            { // open vdata temporary storage file
              std::stringstream strm;
              strm << path << "/"
                   << "tmp_vdata_"
                   << std::setw(3) << std::setfill('0')
                   << rmi.procid()
                   << ".bin";
              std::string fullfn = strm.str();
              atom2vdatafn[i] = fullfn;
              atom2vdatastream[i].open(fullfn.c_str(), 
                                       std::ios::binary | 
                                       std::ios::out |
                                       std::ios::trunc );
              assert(atom2vdatastream[i].good());
            }
            { // open edata temporary storage file
              std::stringstream strm;
              strm << path << "/"
                   << "tmp_edata_"
                   << std::setw(3) << std::setfill('0')
                   << rmi.procid()
                   << ".bin";
              std::string fullfn = strm.str();
              atom2edatafn[i] = fullfn;
              atom2edatastream[i].open(fullfn.c_str(), 
                                       std::ios::binary | 
                                       std::ios::out |
                                       std::ios::trunc );
              assert(atom2edatastream[i].good());
            }               
          } //end of if owning machine
        } //end of for loop
      } //end of constructor
      
      

      void add_vertex_local(const vertex_id_t& vid,
                            const procid_t& atomid,
                            const vertex_color_type& vcolor,
                            const vertex_data_type& vdata) {
        assert(atomid % rmi.numprocs() == rmi.procid());
        assert(atom2atomfile.find(atomid) != atom2atomfile.end());
        atom2atomfile[atomid].globalvids().push_back(vid);
        atom2atomfile[atomid].vcolor().push_back(vcolor);
        atom2atomfile[atomid].atom().puash_back(atomid);
        assert(atom2vdatastream.find(atomid) != atom2vdatastream.end());
        assert(atom2vdatastream[atomid].good());        
        oarchive oarc(atom2vdatastream[atomid]);
        oarc << vdata;        
      }

      void add_vertex(const vertex_id_t& vid,
                      const procid_t& atomid,
                      const vertex_color_type& vcolor,
                      const vertex_data_type& vdata) {
        procid_t owning_machine = atomid % rmi.numprocs();
        if(owning_machine == rmi.procid()) {
          // local add
          add_vertex_local(vid, atomid, vcolor, vdata);
        } else {
          // remote add
          rmi.remote_request(owning_machine, 
                             &atom_shuffler_type::add_vertex_local,
                             vid, atomid, vcolor, vdata);
        }
      } // end of add_vertex


      void add_vertex_boundary_local(const vertex_id_t& vid,
                                     const procid_t& atomid,
                                     const procid_t& from_atomid,
                                     const vertex_color_type& vcolor,
                                     const vertex_data_type& vdata) {
        assert(atomid % rmi.numprocs() == rmi.procid());
        assert(atom2atomfile.find(atomid) != atom2atomfile.end());
        atom2atomfile[atomid].globalvids().push_back(vid);
        atom2atomfile[atomid].vcolor().push_back(vcolor);
        atom2atomfile[atomid].atom().push_back(from_atomid);
        assert(atom2vdatastream.find(atomid) != atom2vdatastream.end());
        assert(atom2vdatastream[atomid].good());        
        oarchive oarc(atom2vdatastream[atomid]);
        oarc << vdata;        
      } // add vertex boundary local


      void add_vertex_boundary(const vertex_id_t& vid,
                               const procid_t& atomid,
                               const procid_t& from_atomid,
                               const vertex_color_type& vcolor,
                               const vertex_data_type& vdata) {
        procid_t owning_machine = atomid % rmi.numprocs();
        if(owning_machine == rmi.procid()) {
          // local add
          add_vertex_boundary_local(vid, atomid, from_atomid, vcolor, vdata);
        } else {
          // remote add
          rmi.remote_request(owning_machine, 
                             &atom_shuffler_type::add_vertex_boundary_local,
                             vid, atomid, from_atomid, vcolor, vdata);
        }
      } // end of add_vertex_boundary

    }; // end of atom shuffler

  }; // end of namespace atom_file_impl;


  // assume 1 per machine
  template <typename VertexData, typename EdgeData>
  void graph_partition_to_atomindex(const std::string& path) {
    typedef VertexData vertex_data_type;
    typedef EdgeData edge_data_type;
    typedef atom_file_impl::atom_shuffler<VertexData, EdgeData>
      atom_shuffler_type;
    // startup distributed control
    dc_init_param param;
    if (init_param_from_env(param) == false) {
      assert(false);
    }
    distributed_control dc(param);
    dc.services().barrier();
    
    // load the maps
    std::vector<vertex_color_type> vertex2color;
    {
      std::string absfname = path + "/coloring.txt";
      std::ifstream fin(absfname);
      while(fin.good()) {
        vertex_color_type vcolor;
        fin >> vcolor;
        if(fin.good()) vertex2color.push_back(vcolor);
      }
      fin.close();
    }
    std::vector<procid_t> vertex2atomid;
    {
      std::string absfname = path + "/partitioning.txt";
      std::ifstream fin(absfname);
      while(fin.good()) {
        procid_t atomid;
        fin >> atomid;
        if(fin.good()) vertex2atomid.push_back(atomid);
      }
      fin.close();
    }

    std::vector<std::string> fnames;
    graph_fragment::list_structure_files(path, fnames);


    // Create the atom shuffler
    atom_shuffler_type atom_shuffler(dc, path);

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
        graphlab::iarchive iarc(fin);        
        // Loop through all the structures reading in the vertex data
        for(vertex_id_t vid = desc.begin_vertex; vid < desc.end_vertex; ++vid) {
          assert(vid < vertex2atomid.size());
          assert(vid < vertex2color.size());
          vertex_data_type vdata;
          iarc >> vdata;
          assert(fin.good());
          atom_shuffler.add_vertex(vid,
                                   vertex2atomid[vid],
                                   vertex2color[vid],
                                   vdata);
        }
        fin.close();
      } // end of read all vdata




    }

  } // end of graph partition to atom index



} // end namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif
