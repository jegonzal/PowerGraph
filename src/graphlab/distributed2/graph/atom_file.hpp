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
#include <graphlab/distributed2/graph/partitioning/raw_fragment.hpp>




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
    idxfile.ncolors = graph.compute_coloring();
    
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
      }
      desc.nverts = atomfile.globalvids().size();
      desc.nedges = atomfile.globaleids().size();
      idxfile.atoms.push_back(desc);
    }
    idxfile.write_to_file(idxfilename);
  }









} // end namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif
