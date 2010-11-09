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
  size_t changes;
  
  jtree_vertex_data() : 
    parent(NULL_VID),
    calibrated(false), 
    sampled(false),
    changes(0)  { }
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






#endif
