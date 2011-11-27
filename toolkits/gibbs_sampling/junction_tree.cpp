/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */



#include "junction_tree.hpp"


#include <graphlab/util/stl_util.hpp>

#include <graphlab/macros_def.hpp>

/**
 * Extend the jtree_list data structure by eliminating the vertex.  If
 * the jtree list can be extended then it is extended and this
 * function returns true.
 *
 **/
bool jtree_list::
extend(const mrf_graph_type::vertex_id_type elim_vertex,
       const mrf_graph_type& mrf,
       const size_t max_tree_width,
       const size_t max_factor_size) {
  typedef mrf_graph_type::edge_id_type edge_id_type;
  typedef mrf_graph_type::vertex_id_type vertex_id_type;
  // sanity check: The vertex to eliminate should not have already
  // been eliminated
  ASSERT_FALSE( contains(elim_vertex) );
  /// =====================================================================
  // Construct the elimination clique for the new vertex
  // 1) Fill out clique
  // 2) Track the corresponding factor size and treewidth
  // 3) Find the parent of this clique
  jtree_list::elim_clique clique;
  clique.elim_vertex = elim_vertex;
  // the factor must at least have the eliminated vertex
  size_t factor_size = 
    std::max(mrf.vertex_data(elim_vertex).variable.size(),
             uint32_t(1));
  foreach(const edge_id_type ineid, mrf.in_edge_ids(elim_vertex)) {
    const vertex_id_type vid = mrf.source(ineid);
    const bool is_in_jtree = contains(vid);
    // if the neighbor is in the set of vertices being eliminated
    if(is_in_jtree) {      
      clique.vertices += vid;
      factor_size *= 
        std::max(mrf.vertex_data(vid).variable.size(), uint32_t(1) );
    }
    // if the clique ever gets too large then teminate
    // the + 1 is because we need to include the elim vertex
    if((clique.vertices.size() > max_tree_width) || 
       (max_factor_size > 0 && factor_size > max_factor_size)) 
      return false;
  }

  // Determine the parent of this clique -------------------------
  vertex_id_type parent_id = 0;
  foreach(vertex_id_type vid, clique.vertices)
    parent_id = std::max(parent_id, elim_time_lookup(vid));
  clique.parent = parent_id;

  /// =====================================================================
  // Simulate injecting vertices in parent cliques back to when RIP is
  // satisfied
  vertex_set rip_verts = clique.vertices;
  for(vertex_id_t parent_vid = clique.parent; 
      !rip_verts.empty() && parent_vid < cliques.size(); ) {
    const jtree_list::elim_clique& parent_clique = cliques[parent_vid];    
    
    // Remove the parent vertex
    rip_verts -= vertex_set(parent_clique.elim_vertex);

    // Construct the new vertices that would normally be stored at
    // this vertes
    const vertex_set tmp_verts = rip_verts + parent_clique.vertices;

    // Check that the expanded clique is still within tree width
    if(tmp_verts.size()  > max_tree_width) return false;

    // If we care about the maximum factor size Compute the factor
    // size and fail if the factor is too large
    if(max_factor_size > 0) {
      size_t factor_size = 
        std::max(mrf.vertex_data(parent_clique.elim_vertex).variable.size(),
                 uint32_t(1));
      foreach(vertex_id_t vid, tmp_verts) {
        factor_size *= 
          std::max(mrf.vertex_data(vid).variable.size(), uint32_t(1));
      }
      if(factor_size > max_factor_size) return false;
    }

    // Find the new parent
    vertex_id_t new_parent_vid = 0;
    foreach(vertex_id_t vid, tmp_verts) {
      new_parent_vid = 
        std::max(new_parent_vid, elim_time_lookup(vid));
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
  const size_t new_elim_time = cliques.size();
  cliques.push_back(clique);
  elim_time[clique.elim_vertex] = new_elim_time;

  /// =====================================================================
  // Satisfy RIP
  rip_verts = clique.vertices;
  for(vertex_id_t parent_vid = clique.parent; 
      !rip_verts.empty() && parent_vid < cliques.size(); ) {
    // get the parent clique
    jtree_list::elim_clique& parent_clique = cliques[parent_vid];       

    // otherwise update that the rip_verts
    rip_verts -= vertex_set(parent_clique.elim_vertex);


    // Construct the new vertices that would normally be stored at
    // this vertes
    const vertex_set tmp_verts = 
      rip_verts + parent_clique.vertices;

    // Determine the new parent (except first vertex)
    vertex_id_t new_parent_vid = 0;
    foreach(vertex_id_t vid, tmp_verts) {
      new_parent_vid = 
        std::max(new_parent_vid, elim_time_lookup(vid));
    }

    //! if the parent changes we must update the rip_verts and the parent value
    if(new_parent_vid != parent_clique.parent) { 
      rip_verts = tmp_verts;
      parent_clique.parent = new_parent_vid;
    } else {
      // If the parent is unchanged then we can remove all the
      // variables stored locally from the rip_verts since they all
      // already satisfy RIP.
      rip_verts -= parent_clique.vertices;
    }
    // update the local vertices
    parent_clique.vertices = tmp_verts;

    // Move up tree
    parent_vid = new_parent_vid;
  }

  // Ensure that the parent of the first clique is the null VID
  cliques.front().parent = NULL_VID;
  // Add successfully
  return true;
} // end of extend clique list




/**
 * Convert a jtree_list into a jtree_graph 
 */
void jtree_list::
load_graph(const mrf_graph_type& mrf,
           const size_t num_factors,
           jtree_graph_type& jt_graph) const {
  //! Todo: Turn this into stack allocated boolean vector
  std::vector<bool> assigned_factors(num_factors, false);
  

  {  // Construct the junction tree
    // Ensure that the parent of the root is identifiable
    ASSERT_EQ(cliques.front().parent, NULL_VID); 
    foreach(const jtree_list::elim_clique& clique, cliques) {      
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
    rev_foreach(const jtree_list::elim_clique& clique, cliques) {
      ASSERT_LT(jt_vid, jt_graph.num_vertices());
      jtree_vertex_data& jt_vdata = jt_graph.vertex_data(jt_vid--);
      const mrf_vertex_data& mrf_vdata = mrf.vertex_data(clique.elim_vertex);
      foreach(factor_id_t fid, mrf_vdata.factor_ids) {
        if(!assigned_factors[fid]) {
          jt_vdata.factor_ids.push_back(fid);
          assigned_factors[fid] = true;
        }
      }
    }
  }
} // end of build junction tree



struct jtgraph_node {
  vertex_id_t vid;
  std::set<vertex_id_t> vars;
  std::set<vertex_id_t> neighbors;
};

typedef std::vector<jtgraph_node> jtgraph_type;


std::ostream& 
operator<<(std::ostream& out, const std::set<vertex_id_t>& set) {
  using namespace graphlab;
  return out << set;
  // out << "{";
  // size_t i = 0; 
  // foreach(const vertex_id_t vid, set) {
  //   out << vid;
  //   if(i + 1 < set.size()) out << ", ";
  //   i++;
  // }
  // out << "}";
  // return out;
}

void test_rip(const jtgraph_type& graph) {
  namespace gl = graphlab;
  std::map<vertex_id_t, 
    std::map<vertex_id_t, std::set< vertex_id_t> > > edgedata;
  std::vector< std::set<vertex_id_t> > reachable(graph.size());
  foreach(const jtgraph_node& node, graph) {
    reachable[node.vid] = node.vars;
    // initiailize the out edges;
    foreach(vertex_id_t nvid, node.neighbors) {
      ASSERT_NE(nvid, node.vid);
      edgedata[node.vid][nvid] = node.vars;
    }
  }
  rev_foreach(const jtgraph_node& node, graph) {
    // Receive in reachable
    foreach(vertex_id_t nvid, node.neighbors) {
      reachable[node.vid].insert(edgedata[nvid][node.vid].begin(),
                                 edgedata[nvid][node.vid].end());
    }
    // write out reachable
    foreach(vertex_id_t nvid, node.neighbors) {
      std::set<vertex_id_t> tmpset =         
        gl::set_difference(reachable[node.vid], edgedata[nvid][node.vid]);
      edgedata[node.vid][nvid].insert(tmpset.begin(), tmpset.end());

    }
  }
  foreach(const jtgraph_node& node, graph) {
    // Receive in reachable
    foreach(vertex_id_t nvid, node.neighbors) {
      reachable[node.vid].insert(edgedata[nvid][node.vid].begin(),
                                 edgedata[nvid][node.vid].end());
    }
    // write out reachable
    foreach(vertex_id_t nvid, node.neighbors) {
      std::set<vertex_id_t> tmpset =         
        gl::set_difference(reachable[node.vid], edgedata[nvid][node.vid]);
      edgedata[node.vid][nvid].insert(tmpset.begin(), tmpset.end());
    }
  }

  // Check the running intersection property
  foreach(const jtgraph_node& node, graph) {
    // std::cout << node.vid << ": " << node.vars << std::endl     
    //           << "\t" << node.neighbors << std::endl;
    std::set<vertex_id_t> local_sep_set;
    foreach(const vertex_id_t n1, node.neighbors) {
      //      std::cout << "\t" << n1 << "--" << edgedata[n1][node.vid] << std::endl;
      foreach(const vertex_id_t n2, node.neighbors) {
        if(n1 != n2) {
          local_sep_set =
            gl::set_union(local_sep_set,
                          gl::set_intersect(edgedata[n1][node.vid],
                                            edgedata[n2][node.vid]));
        }
      }
    }
    //    std::cout << "\t" << local_sep_set << std::endl;
    ASSERT_TRUE(gl::is_subset(local_sep_set, node.vars));
  }
  //  getchar();

}


/**
 * Scan the junction tree list to ensure that all invariants hold.
 */
void jtree_list::
validate() const {
  jtgraph_type jtgraph(cliques.size());
  // validate the junction tree list data structure
  for(size_t i = 0; i < cliques.size(); ++i) {
    const elim_clique& clique = cliques[i];
    ASSERT_EQ(i, elim_time_lookup(clique.elim_vertex));
    ASSERT_FALSE(clique.vertices.contains(clique.elim_vertex));
    const bool is_root = (clique.parent == NULL_VID);
    if(is_root) {
      ASSERT_EQ(i, 0);
    } else {
      ASSERT_GT(i, 0);
      // ensure that the parent is elimated later
      ASSERT_LT(clique.parent, jtgraph.size());
      const elim_clique& parent_clique = cliques[clique.parent];
      ASSERT_TRUE(clique.vertices.contains(parent_clique.elim_vertex));
      ASSERT_LE(clique.vertices, 
                parent_clique.vertices + parent_clique.elim_vertex);
    }
    // populate the jtgraph node
    jtgraph_node& node(jtgraph[i]);
    node.vid = i;
    node.vars.insert(clique.elim_vertex);
    node.vars.insert(clique.vertices.begin(), clique.vertices.end());
    if(!is_root) {
      node.neighbors.insert(clique.parent);   
      jtgraph[clique.parent].neighbors.insert(i);
    }
  }
  // test running intersection property.  
  test_rip(jtgraph);
}







#include <graphlab/macros_undef.hpp>
