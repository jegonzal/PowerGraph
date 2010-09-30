#ifndef PGIBBS_SEQUENTIAL_JT_GIBBS_HPP
#define PGIBBS_SEQUENTIAL_JT_GIBBS_HPP


#include <iostream>
#include <map>
#include <set>
#include <vector>


#include "data_structures.hpp"


#include <graphlab/macros_def.hpp>
  
typedef std::set<vertex_id_t> vertex_set;
typedef std::set<variable_t> variable_set;


size_t build_junction_tree(const mrf::graph_type& mrf, 
                           const vertex_set& in_tree,
                           junction_tree::graph_type& jt) {
  size_t tree_width = 0;
  jt.clear();

  // build the neighbor sets for the variables considered in the blocks
  std::map<vertex_id_t, vertex_set> neighbors;
  std::set<size_t> all_factors;
  // build active neighbors (induced graph of in_tree)
  foreach(vertex_id_t v, in_tree) {
    foreach(edge_id_t eid, mrf.out_edge_ids(v)) {
      vertex_id_t neighbor_v = mrf.target(eid);      
      if(in_tree.count(neighbor_v))  neighbors[v].insert(neighbor_v);
    }
    // Add self edge for simplicity
    neighbors[v].insert(v);
    const mrf::vertex_data& vdata = mrf.vertex_data(v);
    all_factors.insert(vdata.factor_ids.begin(),
                       vdata.factor_ids.end());
  }


  // Construct an elimination ordering:
  graphlab::mutable_queue<vertex_id_t, size_t> elim_order;
  typedef std::pair<vertex_id_t, vertex_set> neighborhood_type;
  foreach(const neighborhood_type& hood, neighbors) {
    elim_order.push(hood.first, in_tree.size() - hood.second.size());
  }

  // Track the child cliques 
  std::map<vertex_id_t, vertex_set> child_cliques;
    
  // Run the elimination;
  while(!elim_order.empty()) {
    const std::pair<vertex_id_t, size_t> top = elim_order.pop();
    const vertex_id_t elim_vertex = top.first;
    const vertex_set& clique_verts = neighbors[elim_vertex];

    // if all the vars in the domain exceed the max factor graph
    // repesentation then we fail.
    if(clique_verts.size() >= MAX_DIM) {
      std::cout << "Clique too large!! ";
      foreach(vertex_id_t vid, clique_verts) std::cout << vid << " ";
      std::cout << std::endl;
      return -1;
    }
    tree_width = std::max(tree_width, clique_verts.size());


    // Start building up the clique data structure
    vertex_id_t clique_id = jt.add_vertex(junction_tree::vertex_data());
    junction_tree::vertex_data& clique = jt.vertex_data(clique_id);

    // Print the state
    std::cout << clique_id << ": [[";

    
    // Add all the variables to the domain
    variable_set vars;
    foreach(vertex_id_t vid, clique_verts) { 
      vars.insert(mrf.vertex_data(vid).variable);
      std::cout << vid << " ";
    }
    std::cout << "]] : ";
    clique.variables = domain_t(vars);
    clique.factor_ids = mrf.vertex_data(elim_vertex).factor_ids;
    

    // construct the union of all possible child cliques
    vertex_set union_child_cliques = child_cliques[elim_vertex];


    // If this clique contains all the remaining variables then we can
    // do all the remaining elimination
    bool finished_elimination = 
      clique_verts.size() > elim_order.size();
    if(finished_elimination ) {
      elim_order.clear();
      // add all the factors
      foreach(vertex_id_t vid, clique_verts) { 
        const mrf::vertex_data& vdata = mrf.vertex_data(vid);
        clique.factor_ids.insert(vdata.factor_ids.begin(), 
                                 vdata.factor_ids.end());
        union_child_cliques.insert(child_cliques[vid].begin(), 
                                   child_cliques[vid].end());
      }
    }

    // Clean up factors
    clique.factor_ids = graphlab::set_intersect(clique.factor_ids,
                                                all_factors);
    
    foreach(size_t fid, clique.factor_ids) all_factors.erase(fid);

 
    // Connect the clique to all of its children cliques
    foreach(vertex_id_t child_id,  union_child_cliques) {
      std::cout << child_id << " ";              
      const junction_tree::vertex_data& child_clique = jt.vertex_data(child_id);
      junction_tree::edge_data edata;
      edata.variables = child_clique.variables.intersect(clique.variables);
      jt.add_edge(child_id, clique_id, edata);
      jt.add_edge(clique_id, child_id, edata);
    }

    std::cout << "( ";
    foreach(size_t fid, clique.factor_ids) std::cout << fid << " ";
    std::cout << ")" << std::endl;

    // If we still have variables to eliminate
    if(!elim_order.empty()) {

      // Disconnect variable and connect neighbors and mark their
      // children cliques
      foreach(const vertex_id_t n_vertex, clique_verts) {
        if(n_vertex != elim_vertex) {
          // connect neighbors
          neighbors[n_vertex].insert(clique_verts.begin(), clique_verts.end());
          // disconnect the variable
          neighbors[n_vertex].erase(elim_vertex);
          // Update the fill count
          elim_order.update(n_vertex, in_tree.size() - neighbors[n_vertex].size());
          // Update the clique neighbors
          child_cliques[n_vertex].insert(clique_id);
          child_cliques[n_vertex] = 
            graphlab::set_difference(child_cliques[n_vertex], 
                                     child_cliques[elim_vertex]);
        }
      }
    }
  }
  // We must have assigned all the factors
  assert(all_factors.empty());

  return tree_width;

} // end of build junction tree



void sample_once(const factorized_model& factor_graph,
                 mrf::graph_type& mrf,
                 size_t root) {


  

  


}











#include <graphlab/macros_undef.hpp>
#endif
