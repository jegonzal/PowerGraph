#ifndef GRAPHLAB_RAW_FRAGMENT
#define GRAPHLAB_RAW_FRAGMENT

/**
 * This file contains the structures used to represent the raw
 * distributed graph information on disk ("prior to running graphlab")
 *
 */


#include <vector>
#include <graphlab/graph/graph.hpp>
#include <graphlab/serialization/serialization_includes.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab {


  struct raw_fragment {

    static const std::string structure_suffix;
    static const std::string vdata_suffix;
    static const std::string edata_suffix;


    std::string structure_filename;
    std::string edge_data_filename;
    std::string vertex_data_filename;

    vertex_id_t id;
    vertex_id_t nfragments;
    vertex_id_t nverts;
    vertex_id_t nedges;
    vertex_id_t fragment_size;
    vertex_id_t fragment_remainder;
    vertex_id_t begin_vertex, end_vertex;
    vertex_id_t num_local_verts;      
    edge_id_t   num_local_edges;


   
    std::vector< std::vector<vertex_id_t> >      in_neighbor_ids;
   
    raw_fragment() { } 
    
    raw_fragment(const std::string& base,
                 vertex_id_t id, 
                 vertex_id_t nfragments,
                 vertex_id_t nverts,
                 vertex_id_t nedges);

    bool is_local(vertex_id_t vid) const;
    vertex_id_t local_id(vertex_id_t gvid) const;
    vertex_id_t owning_fragment(vertex_id_t vid) const;
    void vids_to_fragmentids(const std::vector<vertex_id_t>& vids,
                             std::vector<vertex_id_t>& fragmentids) const;



    bool add_edge(const vertex_id_t source, const vertex_id_t target);   
    
    void finalize_adjacency_lists();

    static void list_structure_files(const std::string& pathname, 
                                      std::vector<std::string>& files);

    void save(graphlab::oarchive& oarc) const;
    void load(graphlab::iarchive& iarc);
    void partial_load(graphlab::iarchive& iarc);
    
  };


  /** Display a graph fragment description */
  std::ostream& operator<<(std::ostream& out, const raw_fragment& desc);







}; // end namespace
#include <graphlab/macros_undef.hpp>
#endif





