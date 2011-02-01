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

    /**
     * Description of the fragment file. This is the first record in the
     * file.
     */
    struct description {
      std::string name;
      vertex_id_t id;
      vertex_id_t nfragments;
      vertex_id_t nverts;
      vertex_id_t nedges;
  
      vertex_id_t fragment_size;
      vertex_id_t fragment_remainder;
      vertex_id_t begin_vertex, end_vertex;

 
      description();
      description(const std::string& base,
                  vertex_id_t id, 
                  vertex_id_t nfragments,
                  vertex_id_t nverts,
                  vertex_id_t nedges);

      bool is_local(vertex_id_t vid) const;
      vertex_id_t owning_fragment(vertex_id_t vid) const;
      vertex_id_t num_local_verts() const;
      void save(graphlab::oarchive& oarc) const;
      void load(graphlab::iarchive& iarc); 
      void vids_to_fragmentids(const std::vector<vertex_id_t>& vids,
                               std::vector<vertex_id_t>& fragmentids) const;

    }; // end of fragment_description




    template<typename VertexData, typename EdgeData>
    struct vertex_record {
      typedef VertexData vertex_data_type;
      typedef EdgeData edge_data_type;

      vertex_id_t vid;
      vertex_data_type vdata;
      
      std::vector<vertex_id_t> out_neighbors;
      std::vector<vertex_id_t> out_neighbors_fragmentid;
      
      std::vector<vertex_id_t> in_neighbors;
      std::vector<vertex_id_t> in_neighbors_fragmentid;
      std::vector<edge_data_type> in_edge_data;
      
      
      void save(graphlab::oarchive& oarc) const {
        oarc << vid
             << vdata
             << out_neighbors
             << out_neighbors_fragmentid
             << in_neighbors
             << in_neighbors_fragmentid
             << in_edge_data;
      }      
      void load(graphlab::iarchive& iarc) {
        iarc >> vid
             >> vdata
             >> out_neighbors
             >> out_neighbors_fragmentid
             >> in_neighbors
             >> in_neighbors_fragmentid
             >> in_edge_data;
      }     
    }; // end of vertex record

  }; // End of graph fragment




}; // end namespace




/** Display a graph fragment description */
std::ostream& 
operator<<(std::ostream& out, 
           const graphlab::graph_fragment::description& desc);








#include <graphlab/macros_undef.hpp>
#endif





