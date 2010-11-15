#ifndef PGIBBS_JUNCTION_TREE_HPP
#define PGIBBS_JUNCTION_TREE_HPP


/**
 *
 * Represents a junction tree
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>





#include <iostream>
#include <iomanip>

#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cassert>


#include <boost/unordered_map.hpp>



#include <graphlab.hpp>


#include "factorized_model.hpp"
#include "mrf.hpp"





struct jtree_vertex_data {
  vertex_id_t parent;
  domain_t variables;
  bool calibrated;
  bool sampled;
  std::vector<factor_id_t> factor_ids;
  factor_t factor;
  assignment_t asg;
  
  jtree_vertex_data() : 
    parent(NULL_VID),
    calibrated(false), 
    sampled(false) { }
}; // End of vertex data


struct jtree_edge_data {
  domain_t variables;
  factor_t message;
  bool calibrated;
  bool received;
  jtree_edge_data() : 
    calibrated(false), 
    received(false) { }
}; // End of edge data


//! define the graph type:
typedef graphlab::graph< jtree_vertex_data, jtree_edge_data> jtree_graph_type;
typedef graphlab::types<jtree_graph_type> jtree_gl;



//// Junction tree construction code
//// =====================================================================>

//! The fast set used in junction tree construction
typedef graphlab::fast_set<2*MAX_DIM, vertex_id_t> vertex_set;




struct jtree_list {
  struct elim_clique {
    //! The parent of this elim clique in the jtree_list
    vertex_id_t parent;
    //! The vertex eliminated when this clique was created 
    vertex_id_t elim_vertex;
    //! The vertices created in this clique EXCLUDING elim_vertex
    vertex_set vertices; 
    elim_clique() : parent(NULL_VID) { }
  };
  typedef std::vector<elim_clique> clique_list_type;
  typedef boost::unordered_map<vertex_id_t, vertex_id_t> elim_time_type;
  //! The collection of cliques
  clique_list_type cliques;
  //! the time variable i was eliminated
  elim_time_type   elim_time;

  bool is_eliminated(const vertex_id_t vid) const {
    return elim_time.find(vid) != elim_time.end();
  }

  void clear() {
    cliques.clear();
    elim_time.clear();
  }

};






/**
 * Extend the jtree_list data structure by eliminating the vertex.  If
 * the jtree list can be extended then it is extended and this
 * function returns true.
 *
 **/
bool extend_jtree_list(const vertex_id_t elim_vertex,
                       const mrf_graph_type& mrf,
                       const size_t max_tree_width,
                       const size_t max_factor_size,
                       jtree_list& jt_list);





/**
 * Convert a jtree_list into a jtree_graph 
 */
void jtree_list_to_jtree_graph(const jtree_list& jt_list,
                               const mrf_graph_type& mrf,
                               const size_t num_factors,
                               jtree_graph_type& jt_graph);





enum JTREE_SDT_KEYS {
  JTREE_FACTOR_MAP_KEY,
  JTREE_MRF_GRAPH_KEY
};


//////////////////////////////////////////////////////////////////////
// routines for managing the factors and the MRF in a junction tree
// graph
void set_factor_map(const factorized_model::factor_map_t& factors,
                    jtree_gl::ishared_data_manager& shared_data_manager);
void set_mrf_graph(mrf_graph_type& mrf_graph,
                   jtree_gl::ishared_data_manager& shared_data_manager);


inline const factorized_model::factor_map_t& 
get_factor_map(const jtree_gl::ishared_data& shared_data) {
  typedef factorized_model::factor_map_t factor_map_t;
  const factor_map_t* fmap_ptr =  
    shared_data.get_constant(JTREE_FACTOR_MAP_KEY).as<const factor_map_t*>();
  assert(fmap_ptr != NULL);
  return *fmap_ptr;
}

inline const factor_t& get_factor(const jtree_gl::ishared_data& shared_data,
                                  const factor_id_t factor_id) {
  return get_factor_map(shared_data)[factor_id];
}

inline size_t get_num_factors(const jtree_gl::ishared_data& shared_data) {
  return get_factor_map(shared_data).size();
}

inline mrf_graph_type& get_mrf_graph(const jtree_gl::ishared_data& shared_data) {
  mrf_graph_type* mrf_graph_ptr =
    shared_data.get_constant(JTREE_MRF_GRAPH_KEY).as<mrf_graph_type*>();
  assert(mrf_graph_ptr != NULL);
  return *mrf_graph_ptr;
}

inline mrf_vertex_data& get_mrf_vdata(const jtree_gl::ishared_data& shared_data,
                                      const vertex_id_t vid) {
  return get_mrf_graph(shared_data).vertex_data(vid);
}







#endif
