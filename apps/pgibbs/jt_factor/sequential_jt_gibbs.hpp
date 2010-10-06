#ifndef PGIBBS_SEQUENTIAL_JT_GIBBS_HPP
#define PGIBBS_SEQUENTIAL_JT_GIBBS_HPP


#include <iostream>
#include <map>
#include <set>
#include <vector>


#include "data_structures.hpp"
#include "update_function.hpp"



#include <graphlab/macros_def.hpp>
  

size_t compute_tree_width(const mrf::graph_type& mrf, 
                          const vertex_set& in_tree) {
  size_t tree_width = 0;

  // build the neighbor sets for the variables considered in the blocks
  std::map<vertex_id_t, vertex_set> neighbors;
  // build active neighbors (induced graph of in_tree)
  foreach(vertex_id_t v, in_tree) {
    foreach(edge_id_t eid, mrf.out_edge_ids(v)) {
      vertex_id_t neighbor_v = mrf.target(eid);      
      if(in_tree.count(neighbor_v))  neighbors[v].insert(neighbor_v);
    }
    // Add self edge for simplicity
    neighbors[v].insert(v);    
  }


  // Construct an elimination ordering:
  graphlab::mutable_queue<vertex_id_t, size_t> elim_order;
  typedef std::pair<vertex_id_t, vertex_set> neighborhood_type;
  foreach(const neighborhood_type& hood, neighbors) {
    elim_order.push(hood.first, in_tree.size() - hood.second.size());
  }
  
  // Run the elimination;
  while(!elim_order.empty()) {
    const std::pair<vertex_id_t, size_t> top = elim_order.pop();
    const vertex_id_t elim_vertex = top.first;
    const vertex_set& clique_verts = neighbors[elim_vertex];

    tree_width = std::max(tree_width, clique_verts.size());
    // if all the vars in the domain exceed the max factor graph
    // repesentation then we fail.
    if(tree_width > MAX_DIM) {
      //      std::cout << "Clique too large!! ";
      // foreach(vertex_id_t vid, clique_verts) std::cout << vid << " ";
      // std::cout << std::endl;
      return -1;
    }

    // If this clique contains all the remaining variables then we can
    // do all the remaining elimination
    bool finished_elimination = clique_verts.size() > elim_order.size();
    if(finished_elimination ) {  elim_order.clear();  }


    // std::cout << "( ";
    // foreach(size_t fid, clique.factor_ids) std::cout << fid << " ";
    // std::cout << ")";
    //    std::cout << std::endl;

    // If we still have variables to eliminate
    if(!elim_order.empty()) {
      foreach(const vertex_id_t n_vertex, clique_verts) {
        if(n_vertex != elim_vertex) {
          vertex_set& neighbor_set = neighbors[n_vertex];
          // connect neighbors
          neighbor_set.insert(clique_verts.begin(), clique_verts.end());
          // disconnect the variable
          neighbor_set.erase(elim_vertex);
          size_t clique_size = neighbor_set.size();
          // if(clique_size > MAX_DIM) return -1;
          // Update the fill count
          elim_order.update(n_vertex, in_tree.size() - clique_size);
        }
      }
    } // if the elim order is not empty
  }
  return tree_width;

} // end of build junction tree



size_t build_junction_tree(const mrf::graph_type& mrf, 
                           const vertex_set& in_tree,
                           junction_tree::graph_type& jt) {
  size_t tree_width = 0;
  jt.clear();

  // build the neighbor sets for the variables considered in the blocks
  std::map<vertex_id_t, vertex_set> neighbors;
  // build active neighbors (induced graph of in_tree)
  foreach(vertex_id_t v, in_tree) {
    foreach(edge_id_t eid, mrf.out_edge_ids(v)) {
      vertex_id_t neighbor_v = mrf.target(eid);      
      if(in_tree.count(neighbor_v))  neighbors[v].insert(neighbor_v);
    }
    // Add self edge for simplicity
    neighbors[v].insert(v);    
  }


  // Construct an elimination ordering:
  graphlab::mutable_queue<vertex_id_t, size_t> elim_order;
  typedef std::pair<vertex_id_t, vertex_set> neighborhood_type;
  foreach(const neighborhood_type& hood, neighbors) {
    elim_order.push(hood.first, in_tree.size() - hood.second.size());
  }
    // Track the child cliques 
  std::map<vertex_id_t, vertex_set> child_cliques;
  // track the used factors
  std::set<size_t> used_factors;
  // Run the elimination;
  while(!elim_order.empty()) {
    const std::pair<vertex_id_t, size_t> top = elim_order.pop();
    const vertex_id_t elim_vertex = top.first;
    const vertex_set& clique_verts = neighbors[elim_vertex];

    tree_width = std::max(tree_width, clique_verts.size());
    // if all the vars in the domain exceed the max factor graph
    // repesentation then we fail.
    if(tree_width > MAX_DIM) {
      //      std::cout << "Clique too large!! ";
      // foreach(vertex_id_t vid, clique_verts) std::cout << vid << " ";
      // std::cout << std::endl;
      return -1;
    }


    // Start building up the clique data structure
    vertex_id_t clique_id = jt.add_vertex(junction_tree::vertex_data());
    junction_tree::vertex_data& clique = jt.vertex_data(clique_id);

    // Print the state
    //    std::cout << clique_id << ": [[";

    
    // Add all the variables to the domain
    foreach(vertex_id_t vid, clique_verts) { 
      clique.variables += mrf.vertex_data(vid).variable;     
      //      std::cout << vid << " ";
    }
    //    std::cout << "]] : ";

    // Add all the unused factors associated with the eliminated variable
    foreach(size_t fid, mrf.vertex_data(elim_vertex).factor_ids) {
      if(used_factors.count(fid) == 0) {
        clique.factor_ids.insert(fid);
      }
    }
    used_factors.insert(clique.factor_ids.begin(), clique.factor_ids.end());
       

    // Define the union of child cliques to the all the child cliques
    // of the eliminated variable.
    vertex_set& union_child_cliques = child_cliques[elim_vertex];


    // If this clique contains all the remaining variables then we can
    // do all the remaining elimination
    bool finished_elimination = clique_verts.size() > elim_order.size();
    if(finished_elimination ) {
      elim_order.clear();
      // add all the remaining factors
      foreach(vertex_id_t vid, clique_verts) { 
        const mrf::vertex_data& vdata = mrf.vertex_data(vid);
        // add all the unsued factor ids
        foreach(size_t fid, vdata.factor_ids) {
          if(used_factors.count(fid) == 0) { 
            clique.factor_ids.insert(fid); 
          }          
        }
        union_child_cliques.insert(child_cliques[vid].begin(), 
                                   child_cliques[vid].end());
      }
    }

 
    // Connect the clique to all of its children cliques
    foreach(vertex_id_t child_id,  union_child_cliques) {
      //      std::cout << child_id << " ";              
      const junction_tree::vertex_data& child_clique = jt.vertex_data(child_id);
      junction_tree::edge_data edata;
      edata.variables = 
        child_clique.variables.intersect(clique.variables);
      // edata.variables = 
      //   graphlab::set_intersect(child_clique.variables,
      //                           clique.variables);
      jt.add_edge(child_id, clique_id, edata);
      jt.add_edge(clique_id, child_id, edata);
    }

    // std::cout << "( ";
    // foreach(size_t fid, clique.factor_ids) std::cout << fid << " ";
    // std::cout << ")";
    //    std::cout << std::endl;

    // If we still have variables to eliminate
    if(!elim_order.empty()) {
      // Disconnect variable and connect neighbors and mark their
      // children cliques
      vertex_set& elim_clique_set = child_cliques[elim_vertex];
      foreach(const vertex_id_t n_vertex, clique_verts) {
        if(n_vertex != elim_vertex) {
          vertex_set& neighbor_set = neighbors[n_vertex];
          // connect neighbors
          neighbor_set.insert(clique_verts.begin(), clique_verts.end());
          // disconnect the variable
          neighbor_set.erase(elim_vertex);
          // if the neighbor clique size becomes too large then we are done
          size_t clique_size = neighbor_set.size();
          // if(clique_size > MAX_DIM) return -1;
          // Update the fill count
          elim_order.update(n_vertex, in_tree.size() - clique_size);
          // Update the clique neighbors
          vertex_set& child_clique_set = child_cliques[n_vertex];
          child_clique_set.insert(clique_id);
          foreach(vertex_id_t vid, elim_clique_set) {
            child_clique_set.erase(vid);
          }
        }
      }
    } // if the elim order is not empty
  }
  return tree_width;

} // end of build junction tree




size_t build_junction_tree(const mrf::graph_type& mrf,
                           vertex_id_t root,
                           junction_tree::graph_type& jt) {

  vertex_set block;
  std::queue<vertex_id_t> bfs_queue;
  vertex_set visited;
  jt.clear();

  size_t tree_width = 0;
  //  size_t looplimit = 1000;
  // add the root
  bfs_queue.push(root);
  visited.insert(root);
  while(!bfs_queue.empty()) {
    //    looplimit--;
    //    if (looplimit == 0) break;
    vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop();
    block.insert(next_vertex);
    // build a junction tree
    tree_width = compute_tree_width(mrf, block);
    if(tree_width <= MAX_DIM) {
      // add the neighbors to the search queue
      foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
        vertex_id_t neighbor_vid = mrf.target(eid);
        if(visited.count(neighbor_vid) == 0) {
          bfs_queue.push(neighbor_vid);
          visited.insert(neighbor_vid);
        }
      }
    } else {
      // remove variable since tree is too large
      block.erase(next_vertex);
    }
  }
  std::cout << "Varcount: " << block.size() << std::endl;
  // Rebuild the junction tree 
  tree_width = build_junction_tree(mrf, block, jt);
  return tree_width;
}








void sample_once(const factorized_model& factor_graph,
                 mrf::graph_type& mrf,
                 vertex_id_t root) {
 
  junction_tree::gl::core jt_core;
  std::cout << "Building Tree" << std::endl;
  size_t tree_width = build_junction_tree(mrf, root, jt_core.graph());
  std::cout << "Root:  " << root << " ----------------" << std::endl;
  std::cout << "Tree width: " << tree_width << std::endl;

  // Setup the core
  jt_core.set_scheduler_type("fifo");
  jt_core.set_scope_type("edge");
  jt_core.set_ncpus(8);
  jt_core.set_engine_type("async");
 

  // Setup the shared data
  typedef factorized_model::factor_map_t factor_map_t;
  const factor_map_t* ptr = & factor_graph.factors();
  jt_core.shared_data().set_constant(junction_tree::FACTOR_KEY, 
                                     ptr);
  jt_core.shared_data().set_constant(junction_tree::MRF_KEY, 
                                     &mrf);

  std::cout << "Running Inference" << std::endl;
  // Calibrate the tree
  jt_core.add_task_to_all(junction_tree::calibrate_update, 1.0);
  jt_core.start();


  // for(vertex_id_t i = 0; i < jt_core.graph().num_vertices(); ++i) {
  //   const junction_tree::vertex_data& vdata = jt_core.graph().vertex_data(i);
  //   std::cout << i << ": " << vdata.sampled << " :--> ";
  //   foreach(edge_id_t eid, jt_core.graph().out_edge_ids(i)) {
  //     std::cout << jt_core.graph().target(eid) << " ";
  //   }
  //   std::cout << std::endl;
  // }




  // Ensure entire tree is sampled
  for(vertex_id_t i = 0; i < jt_core.graph().num_vertices(); ++i) {
    assert(jt_core.graph().vertex_data(i).sampled);
  }



  //   // Schedule sampling starting at the lass vetex
  //   jt_core.add_task(jt_core.graph().num_vertices()-1, 
  //                    junction_tree::sample_update, 1.0);
 

  
  

}








#include <graphlab/macros_undef.hpp>
#endif
