#ifndef GRAPHLAB_DISTRIBUTED_ATOM_FILE_HPP
#define GRAPHLAB_DISTRIBUTED_ATOM_FILE_HPP
#include <vector>
#include <string>
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

namespace graphlab {

/**
 * The contents of an atom file.
 */ 
template <typename VertexData, typename EdgeData>
class atom_file {
 public:
  atom_file():protocol_(protocol), filename_(filename), iarc(NULL), loadstage(0){
  }

  /**
   * Associates this object with a particular atom file
   * and a particular protocol. The contents of the file (through
   * the accessors) are not yet loaded until the load functions are 
   * called. \see load_id_maps() \see load_structure() \see load_all()
   */
  void set_filename(std::string protocol, std::string filename) {
    if (protocol == "file") {
      fin.open(filename.c_str(), "b");
      iarc = new iarchive(fin);
      loadstage = 0;
    }
    else {
      logger(LOG_FATAL, "Unrecognized protocol: %s", protocol.c_str());
    }
  }
  /**
   * Only load the globalvids and the globaleids
   * All remaining entries are not loaded.
   */
  void load_id_maps() {
    if (loadstage == 0) {
      (*iarc) >> globalvids_ >> globaleids_;
      ++loadstage;
    }
  }
  
  /**
   * Loads all entries except the vdata and edata entries.
   */
  void load_structure() {
    if (loadstage < 1) load_id_maps();
    if (loadstage == 1) {
      (*iarc) >> atom_ >> vcolor_ >> in_edges_ >> out_edges_ >> edge_src_dest_;
      ++loadstage;
    }
  }
  
  /**
   * Completely loads the file. All datastructures are now accessible.
   */
  void load_all() {
    if (loadstage < 2) load_structure();
    if (loadstage == 2) {
      (*iarc) >> vdata_ >> edata_;
      ++loadstage;
    }

  }
  
  /**
   * Clears all loaded data and closes the file.
   */
  void clear() {
    delete iarc;
    fin.close();
    globalvids_.clear();
    globaleids_.clear();
    atom_.clear();
    vcolor_.clear();
    in_edges_.clear();
    out_edges_.clear();
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
  inline const std::vector< std::vector<edge_id_t> >& in_edges() const { return in_edges_; }
  inline const std::vector< std::vector<edge_id_t> >& out_edges() const { return out_edges_; }
  inline const std::vector< std::pair<vertex_id_t, vertex_id_t> >& edge_src_dest() const { return edge_src_dest_; }
  inline const std::vector<VertexData>& vdata() const { return vdata_; }
  inline const std::vector<EdgeData>& edata() const { return edata_; }


private:
  std::string protocol_;
  std::string filename_;
  
  // for file protocol
  ifstream fin;
  iarchive *iarc;
  size_t loadstage;
  
  std::vector<vertex_id_t> globalvids_;
  std::vector<edge_id_t> globaleids_;
  std::vector<procid_t> atom_;
  std::vector<vertex_color_type> vcolor_;
  std::vector< std::vector<edge_id_t> >  in_edges_;
  std::vector< std::vector<edge_id_t> >  out_edges_;
  std::vector< std::pair<vertex_id_t, vertex_id_t> > edge_src_dest_;
  std::vector<VertexData> vdata_;
  std::vector<EdgeData> edata_;
};


} // end namespace graphlab
#endif
