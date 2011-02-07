#ifndef GRAPHLAB_GRAPH_FRAGMENT
#define GRAPHLAB_GRAPH_FRAGMENT

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


  namespace graph_fragment {

    extern const std::string structure_suffix;
    extern const std::string vdata_suffix;
    extern const std::string edata_suffix;

    /**
     * Description of the fragment file. This is the first record in the
     * file.
     */
    struct file_description {
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
      edge_id_t num_local_edges;
      
      file_description();
      file_description(const std::string& base,
                       vertex_id_t id, 
                       vertex_id_t nfragments,
                       vertex_id_t nverts,
                       vertex_id_t nedges);

      bool is_local(vertex_id_t vid) const;
      vertex_id_t local_id(vertex_id_t gvid) const;
      vertex_id_t owning_fragment(vertex_id_t vid) const;
      void vids_to_fragmentids(const std::vector<vertex_id_t>& vids,
                               std::vector<vertex_id_t>& fragmentids) const;

      void save(graphlab::oarchive& oarc) const;
      void load(graphlab::iarchive& iarc); 
    }; // end of fragment_description



    struct structure_description {
      file_description desc;
      std::vector< std::vector<vertex_id_t> >      neighbor_ids;
      std::vector< std::vector<vertex_id_t> >      in_neighbor_ids;
      std::vector< std::vector<edge_id_t> >        in_edge_ids;


      structure_description();
      structure_description(const file_description& desc);
      
      /**
       * add an edge to the neighbor ids and edge ids.  if the target
       * is local then return true. 
       *
       */
      bool add_edge(const vertex_id_t source, const vertex_id_t target);   
      void finalize();


      void save(graphlab::oarchive& oarc) const;
      void load(graphlab::iarchive& iarc); 
    }; // end of structure_description
      

  }; // End of graph fragment




}; // end namespace




/** Display a graph fragment description */
std::ostream& 
operator<<(std::ostream& out, 
           const graphlab::graph_fragment::file_description& desc);








#include <graphlab/macros_undef.hpp>
#endif





