#ifndef GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE
#define GRAPHLAB_DISTRIBUTED_ATOM_INDEX_FILE

#include <string>
#include <vector>
#include <graphlab/graph/graph.hpp>

namespace graphlab {


  struct atom_file_descriptor{
    std::string protocol;
    std::string file;
    vertex_id_t nverts;
    edge_id_t nedges;
    std::vector<vertex_id_t> adjatoms;
    std::vector<vertex_id_t> optional_weight_to_adjatoms;
    void save(oarchive &oarc) const{
      oarc << protocol
           << file
           << nverts
           << nedges
           << adjatoms
           << optional_weight_to_adjatoms;
    }
    
    void load(iarchive &iarc) {
      iarc >> protocol
           >> file
           >> nverts
           >> nedges
           >> adjatoms
           >> optional_weight_to_adjatoms;
    }
  };

  struct atom_index_file {

    size_t nverts, nedges, natoms;
    std::vector<atom_file_descriptor> atoms;

    void read_from_file(std::string indexfile);
    void write_to_file(std::string outfilename);
  };



  std::vector<std::vector<size_t> >
  partition_atoms(const atom_index_file& atomindex, size_t nparts);


  /**
   * This parallel function constructs an atom index by reading all
   * the atom files stored at path.  This relies on the distributed
   * comm layer.
   */ 
  void build_atom_index_file(const std::string& path);


} // end namespace graphlab

#endif
