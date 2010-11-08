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






#include <graphlab/macros_def.hpp>




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
                       jtree_list& jt_list) {
  // sanity check: The vertex to eliminate should not have already
  // been eliminated
  assert(jt_list.elim_time.find(elim_vertex) == jt_list.elim_time.end());


  /// =====================================================================
  // Construct the elimination clique for the new vertex
  // 1) Fill out clique
  // 2) Track the corresponding factor size and treewidth
  // 3) Find the parent of this clique
  jtree_list::elim_clique clique;
  clique.elim_vertex = elim_vertex;
  // the factor must at least have the eliminated vertex
  size_t factor_size = 
    std::max(mrf.vertex_data(elim_vertex).variable.arity,
             uint32_t(1));
  foreach(const edge_id_t ineid, mrf.in_edge_ids(elim_vertex)) {
    const vertex_id_t vid = mrf.source(ineid);
    const bool is_in_jtree = 
      jt_list.elim_time.find(vid) != jt_list.elim_time.end();
    // if the neighbor is in the set of vertices being eliminated
    if(in_in_jtree) {      
      clique.vertices += vid;
      factor_size *= 
        std::max(mrf.vertex_data(vid).variable.arity, uint32_t(1) );
    }
    // if the clique ever gets too large then teminate
    // the + 1 is because we need to include the elim vertex
    if(clique.vertices.size() > max_tree_width) return false;
    if(factor_size > max_factor_size) return false;
  }

  // Determine the parent of this clique -------------------------
  vertex_id_t parent_id = 0;
  foreach(vertex_id_t vid, clique.vertices)
    parent_id = std::max(parent_id, jt_list.elim_time[vid]);
  clique.parent = parent_id;

  /// =====================================================================
  // Simulate injecting vertices in parent cliques back to when RIP is
  // satisfied
  vertex_set rip_verts = clique.vertices;
  for(vertex_id_t parent_vid = clique.parent; 
      !rip_verts.empty() && parent_vid < jt_list.cliques.size(); ) {
    const elim_clique& parent_clique = jt_list.cliques[parent_vid];    
    
    // Remove the parent vertex
    rip_verts -= parent_clique.elim_vertex;

    // Construct the new vertices that would normally be stored at
    // this vertes
    const vertex_set tmp_verts = rip_verts + parent_clique.vertices;

    // Check that the expanded clique is still within tree width
    if(tmp_verts.size()  > max_tree_width) return false;

    // Compute the factor size
    size_t factor_size = 
      std::max(mrf.vertex_data(parent_clique.elim_vertex).variable.arity,
               uint32_t(1));
    foreach(vertex_id_t vid, tmp_verts) {
      factor_size *= 
        std::max(mrf.vertex_data(vid).variable.arity, uint32_t(1));
    }
    if(factor_size > max_factor_size) return false;

    // Find the new parent
    vertex_id_t new_parent_vid = 0;
    foreach(vertex_id_t vid, tmp_verts) {
      new_parent_vid = 
        std::max(new_parent_vid, jt_list.elim_time[vid]);
    }
    // if the parent changes then we may need to update RIP with
    // tmp_verts otherwise we use rip_verts
    if(new_parent_vid != parent_clique.parent) 
      rip_verts = tmp_verts;
    else
      rip_verts -= parent_clique.vertices;
    
    // move up the tree
    parent_vid = new_parent_vid;
  }


  /// =====================================================================
  // Assert that if we reached this point RIP can be satisfied safely
  // so proceed to update local data structures
  const size_t elim_time = jt_list.cliques.size();
  jt_list.cliques.push_back(clique);
  jt_list.elim_time[clique.elim_vertex] = elim_time;

  /// =====================================================================
  // Satisfy RIP
  rip_verts = clique.vertices;
  for(vertex_id_t parent_vid = clique.parent; 
      !rip_verts.empty() && parent_vid < jt_list.cliques.size(); ) {
    // get the parent clique
    elim_clique& parent_clique = jt_list.cliques[parent_vid];       

    // otherwise update that the rip_verts
    rip_verts -= parent_clique.elim_vertex;

    // Update the clique
    parent_clique.vertices += rip_verts;

    // Determine the new parent (except first vertex)
    vertex_id_t new_parent_vid = 0;
    foreach(vertex_id_t vid, parent_clique.vertices) {
      new_parent_vid = 
        std::max(new_parent_vid, jt_list.elim_time[vid]);
    }

    //! if the parent changes we must update the rip_verts and the parent value
    if(new_parent_vid != parent_clique.parent) { 
      rip_verts = parent_clique.vertices;
      parent_clique.parent = new_parent_vid;
    } else {
      // If the parent is unchanged then we can remove all the
      // variables stored locally from the rip_verts since they all
      // already satisfy RIP.
      rip_verts -= parent_clique.vertices;
    }
    // Move up tree
    parent_vid = new_parent_vid;
  }

  // Ensure that the parent of the first clique is the null VID
  jt_list.cliques.front().parent = NULL_VID;

  // Add successfully
  return true;
} // end of extend clique list














void jtree_list_to_jtree_graph(const jtree_list& jt_list,
                               const mrf_graph_type& mrf,
                               const size_t num_factors,
                               jtree_graph_type& jt_graph) {
  
  //! Todo: Turn this into stack allocated boolean vector
  std::vector<bool> assigned_factors(num_factors, false);
  

  {  // Construct the junction tree
    // Ensure that the parent of the root is identifiable
    assert(jt_list.cliques.front().parent == NULL_VID); 
    foreach(const elim_clique& clique, jt_list.cliques) {      
      const mrf_vertex_data& elim_vertex_vdata = 
        mrf.vertex_data(clique.elim_vertex);
      // Create the vertex data
      jtree_vertex_data vdata;
      // set the vertex parent
      vdata.parent = clique.parent;
      // add the eliminated vertex
      vdata.variables = elim_vertex_vdata.variable;
      // add all the other variables in the clique
      foreach(vertex_id_t vid, clique.vertices) 
        vdata.variables += mrf.vertex_data(vid).variable;      
      // Add the vertex to the junction tree
      vertex_id_t child_id = jt_graph.add_vertex(vdata);
      // get the cliques parent
      vertex_id_t parent_id = clique.parent;
      // Add the edge to parent if not root
      if(parent_id != NULL_VID) {
        // Get the parent vertex data
        const jtree_vertex_data& parent_vdata =
          jt_graph.vertex_data(parent_id);
        jtree_edge_data edata;
        edata.variables = 
          vdata.variables.intersect(parent_vdata.variables);
        // Add the actual edges
        jt_graph.add_edge(child_id, parent_id, edata);
        jt_graph.add_edge(parent_id, child_id, edata);
      }
    } // end of for each
  } // End of construct cliques


  { // Assign factors 
    // Very important that these be assigned in reverse order
    size_t jt_vid = jt_graph.num_vertices() - 1;
    rev_foreach(const elim_clique& clique, jt_list.cliques) {
      assert(jt_vid < jt_graph.num_vertices());
      jtree_vertex_data& jt_vdata = jt_graph.vertex_data(jt_vid--);
      const mrf_vertex_data& mrf_vdata = 
        mrf.vertex_data(clique.elim_vertex);
      foreach(factor_id_t fid, mrf_vdata.factor_ids) {
        if(!assigned_factors[fid]) {
          jt_vdata.factor_ids.push_back(fid);
          assigned_factors[fid] = true;
        }
      }
    }
  }
} // end of build junction tree
















#include <graphlab/macros_undef.hpp>
#endif
