#ifndef PGIBBS_SEQUENTIAL_JT_GIBBS_HPP
#define PGIBBS_SEQUENTIAL_JT_GIBBS_HPP


#include <iostream>
#include <map>
#include <set>
#include <vector>


#include "data_structures.hpp"
#include "update_function.hpp"
#include "image.hpp"


#include <graphlab/macros_def.hpp>
  


// typedef graphlab::fast_set<2*MAX_DIM, vertex_id_t> vertex_set;
typedef std::set<vertex_id_t> vertex_set;


struct clique_type {
  vertex_set factors;
  vertex_set vars;
  vertex_set children;
};



template<typename A, typename B>
const B& safe_find(const std::map<A, B>& const_map, const A& key) {
  typedef typename std::map<A, B>::const_iterator const_iterator;
  const_iterator iter = const_map.find(key);
  assert(iter != const_map.end());
  return iter->second;
}

typedef std::map<vertex_id_t, vertex_set> vset_map;






size_t compute_tree_width(const vset_map& var2factors_const,
                          const vset_map& factor2vars_const,
                          std::vector<vertex_id_t>* elim_order = NULL,
                          std::vector<vertex_set>* clique_sets = NULL) {

  // Reset elim order and clique set (if they are available)
  if(elim_order != NULL) elim_order->clear();
  if(clique_sets != NULL) clique_sets->clear();

  // Make local copies of the maps
  vset_map var2factors(var2factors_const);
  vset_map factor2vars(factor2vars_const);

  // track the treewidth
  size_t tree_width = 0;

  // Construct an elimination ordering:
  graphlab::mutable_queue<vertex_id_t, int> elim_priority_queue;
  typedef vset_map::value_type vset_map_pair;
  foreach(const vset_map_pair& pair, var2factors) {
    vertex_id_t var = pair.first;
    const vertex_set& factors = pair.second;
    size_t fill_edges = 0;
    foreach(vertex_id_t fid, factors) 
      fill_edges += ( factor2vars[fid].size() - 1);
    elim_priority_queue.push(var, -fill_edges);
    //std::cout << "Fill: " << var << ":  " << fill_edges << std::endl;
  }

  // keep track of the next unique factor id value (used to created
  // temporary factors along the way).
  size_t next_new_factor_id = factor2vars.rbegin()->first + 1;

  // vertex_set used_factors;
  // std::map<vertex_id_t, clique_type> cliques;

  // Run the elimination;
  while(!elim_priority_queue.empty()) {
    const std::pair<vertex_id_t, float> top = elim_priority_queue.pop();
    const vertex_id_t elim_vertex = top.first;
    const vertex_set& factorset = safe_find(var2factors, elim_vertex);
    
    //    std::cout << "Eliminating: " << elim_vertex << " in ";
    //     foreach(vertex_id_t fid, factorset)
    //       std::cout << fid << " ";
    //     std::cout << std::endl;

    // Track the affected vertices
    vertex_set affected_vertices;

    // Erase the variable from all of its factors
    foreach(vertex_id_t fid, factorset) {
      vertex_set& vset = factor2vars[fid];
      vset.erase(elim_vertex);
      affected_vertices.insert(vset.begin(), vset.end()); 
      // First make sure that the treewidth hasn't gotten too large
      tree_width = std::max(tree_width, affected_vertices.size() + 1);
      if(tree_width > MAX_DIM) {
        return tree_width;
      }
    }

    // if necessary store the elimination ordering and the clique set
    if(elim_order != NULL) elim_order->push_back(elim_vertex);
    if(clique_sets != NULL) clique_sets->push_back(affected_vertices);



    // Merge any factors
    if(factorset.size() > 1) {
      //      std::cout << "Merging: ---------------------------------" << std::endl;
      
      // Build the new factor
      size_t new_factor_id = next_new_factor_id++;
      factor2vars[new_factor_id] = affected_vertices;

      // Merge all the nonconstant factors and create a new factor
      foreach(vertex_id_t fid, factorset) {
        const vertex_set& other_factor_vset = factor2vars[fid];
        // Disconnect the old factor and reconnect the new factor
        foreach(vertex_id_t vid, other_factor_vset) {
          var2factors[vid].erase(fid);
          var2factors[vid].insert(new_factor_id);
        }
      }

      //       // Display the new factor
      //       std::cout << std::endl;
      //       std::cout << "Factor: (" << new_factor_id << ") " << elim_vertex << " ";
      //       foreach(vertex_id_t vid, affected_vertices) 
      //         std::cout << vid << " ";      
      //       std::cout << "[ ";
      //       foreach(vertex_id_t fid, factorset) 
      //         std::cout << fid << " ";
      //       std::cout << "]" << std::endl;
      //       std::cout << std::endl;

    } // end of merge

    // Update the fill order
    foreach(vertex_id_t v, affected_vertices) {
      size_t fill_edges = 0;
      const vertex_set& factors = var2factors[v];
      foreach(vertex_id_t fid, factors)
        fill_edges += (factor2vars[fid].size() - 1);
      elim_priority_queue.update(v, -fill_edges);
    }

  }
  return tree_width;
} // end of compute treewidth






size_t build_junction_tree(const mrf::graph_type& mrf,
                           vertex_id_t root,
                           junction_tree::graph_type& jt) {
  jt.clear();

  vset_map var2factors;
  vset_map factor2vars;


  std::queue<vertex_id_t> bfs_queue;
  std::set<vertex_id_t> visited;

  std::vector<vertex_id_t> elim_order;

  size_t tree_width = 0;
  //  size_t looplimit = 1000;
  // add the root
  bfs_queue.push(root);
  visited.insert(root);

  while(!bfs_queue.empty()) {
    //    looplimit--;
    //    if (looplimit == 0) break;
    const vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop();

    // Update data structures 
    const mrf::vertex_data& vdata = mrf.vertex_data(next_vertex);
    var2factors[next_vertex] = vdata.factor_ids;
    foreach(vertex_id_t fid, vdata.factor_ids) {
      factor2vars[fid].insert(next_vertex); 
    }

    // build a junction tree
    tree_width = compute_tree_width(var2factors, factor2vars, &elim_order);
    std::cout << "Tree_width: " << tree_width;
    if(tree_width <= MAX_DIM) {
      std::cout << std::endl;
      // add the neighbors to the search queue
      foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
        vertex_id_t neighbor_vid = mrf.target(eid);
        if(visited.count(neighbor_vid) == 0) {
          bfs_queue.push(neighbor_vid);
          visited.insert(neighbor_vid);
        }
      }
    } else {
      std::cout << "-----------FAILED-----------------" << std::endl;
      // remove the variable if we decide not to use it
      const mrf::vertex_data& vdata = mrf.vertex_data(next_vertex);
      var2factors.erase(next_vertex);
      foreach(vertex_id_t fid, vdata.factor_ids) {
        factor2vars[fid].erase(next_vertex); 
      } 

    }
  }

  tree_width = compute_tree_width(var2factors, factor2vars, &elim_order);

  image img(200, 200);
//   typedef vset_map::value_type pair_type;
//   foreach(const pair_type& pair, var2factors) {
  size_t index = 1;
  foreach(vertex_id_t vid, elim_order) {
    img.pixel(vid) = index++;
  }
  img.save("tree.pgm");


  std::cout << "Varcount: " << var2factors.size() << std::endl;
  //   // Rebuild the junction tree 
  //   tree_width = build_junction_tree(mrf, block, jt);
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

  std::cout << "Done!!!" << std::endl;
  return;

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
  

} // sample once











// size_t compute_tree_width(const mrf::graph_type& mrf, 
//                           const vertex_set& in_tree) {
//   size_t tree_width = 0;

//   // build the neighbor sets for the variables considered in the blocks
//   std::map<vertex_id_t, vertex_set> neighbors;
//   // build active neighbors (induced graph of in_tree)
//   foreach(vertex_id_t v, in_tree) {
//     foreach(edge_id_t eid, mrf.out_edge_ids(v)) {
//       vertex_id_t neighbor_v = mrf.target(eid);      
//       if(in_tree.count(neighbor_v))  neighbors[v].insert(neighbor_v);
//     }
//     // Add self edge for simplicity
//     neighbors[v].insert(v);    
//   }


//   // Construct an elimination ordering:
//   graphlab::mutable_queue<vertex_id_t, size_t> elim_order;
//   typedef std::pair<vertex_id_t, vertex_set> neighborhood_type;
//   foreach(const neighborhood_type& hood, neighbors) {
//     elim_order.push(hood.first, in_tree.size() - hood.second.size());
//   }
  
//   // Run the elimination;
//   while(!elim_order.empty()) {
//     const std::pair<vertex_id_t, size_t> top = elim_order.pop();
//     const vertex_id_t elim_vertex = top.first;
//     const vertex_set& clique_verts = neighbors[elim_vertex];

//     tree_width = std::max(tree_width, clique_verts.size());
//     // if all the vars in the domain exceed the max factor graph
//     // repesentation then we fail.
//     if(tree_width > MAX_DIM) {
//       //      std::cout << "Clique too large!! ";
//       // foreach(vertex_id_t vid, clique_verts) std::cout << vid << " ";
//       // std::cout << std::endl;
//       return -1;
//     }

//     // If this clique contains all the remaining variables then we can
//     // do all the remaining elimination
//     bool finished_elimination = clique_verts.size() > elim_order.size();
//     if(finished_elimination ) {  elim_order.clear();  }


//     // std::cout << "( ";
//     // foreach(size_t fid, clique.factor_ids) std::cout << fid << " ";
//     // std::cout << ")";
//     //    std::cout << std::endl;

//     // If we still have variables to eliminate
//     if(!elim_order.empty()) {
//       foreach(const vertex_id_t n_vertex, clique_verts) {
//         if(n_vertex != elim_vertex) {
//           vertex_set& neighbor_set = neighbors[n_vertex];
//           // connect neighbors
//           neighbor_set.insert(clique_verts.begin(), clique_verts.end());
//           // disconnect the variable
//           neighbor_set.erase(elim_vertex);
//           size_t clique_size = neighbor_set.size();
//           // if(clique_size > MAX_DIM) return -1;
//           // Update the fill count
//           elim_order.update(n_vertex, in_tree.size() - clique_size);
//         }
//       }
//     } // if the elim order is not empty
//   }
//   return tree_width;
// } // end of compute_tree_width



















// size_t build_junction_tree(const mrf::graph_type& mrf, 
//                            const vertex_set& in_tree,
//                            junction_tree::graph_type& jt) {
//   size_t tree_width = 0;
//   jt.clear();

//   // build the neighbor sets for the variables considered in the blocks
//   std::map<vertex_id_t, vertex_set> neighbors;
//   // build active neighbors (induced graph of in_tree)
//   foreach(vertex_id_t v, in_tree) {
//     foreach(edge_id_t eid, mrf.out_edge_ids(v)) {
//       vertex_id_t neighbor_v = mrf.target(eid);      
//       if(in_tree.count(neighbor_v))  neighbors[v].insert(neighbor_v);
//     }
//     // Add self edge for simplicity
//     neighbors[v].insert(v);    
//   }


//   // Construct an elimination ordering:
//   graphlab::mutable_queue<vertex_id_t, size_t> elim_order;
//   typedef std::pair<vertex_id_t, vertex_set> neighborhood_type;
//   foreach(const neighborhood_type& hood, neighbors) {
//     elim_order.push(hood.first, in_tree.size() - hood.second.size());
//   }
//     // Track the child cliques 
//   std::map<vertex_id_t, vertex_set> child_cliques;
//   // track the used factors
//   std::set<size_t> used_factors;
//   // Run the elimination;
//   while(!elim_order.empty()) {
//     const std::pair<vertex_id_t, size_t> top = elim_order.pop();
//     const vertex_id_t elim_vertex = top.first;
//     const vertex_set& clique_verts = neighbors[elim_vertex];

//     tree_width = std::max(tree_width, clique_verts.size());
//     // if all the vars in the domain exceed the max factor graph
//     // repesentation then we fail.
//     if(tree_width > MAX_DIM) {
//       //      std::cout << "Clique too large!! ";
//       // foreach(vertex_id_t vid, clique_verts) std::cout << vid << " ";
//       // std::cout << std::endl;
//       return -1;
//     }


//     // Start building up the clique data structure
//     vertex_id_t clique_id = jt.add_vertex(junction_tree::vertex_data());
//     junction_tree::vertex_data& clique = jt.vertex_data(clique_id);

//     // Print the state
//     //    std::cout << clique_id << ": [[";

    
//     // Add all the variables to the domain
//     foreach(vertex_id_t vid, clique_verts) { 
//       clique.variables += mrf.vertex_data(vid).variable;     
//       //      std::cout << vid << " ";
//     }
//     //    std::cout << "]] : ";

//     // Add all the unused factors associated with the eliminated variable
//     foreach(size_t fid, mrf.vertex_data(elim_vertex).factor_ids) {
//       if(used_factors.count(fid) == 0) {
//         clique.factor_ids.insert(fid);
//       }
//     }
//     used_factors.insert(clique.factor_ids.begin(), clique.factor_ids.end());
       

//     // Define the union of child cliques to the all the child cliques
//     // of the eliminated variable.
//     vertex_set& union_child_cliques = child_cliques[elim_vertex];


//     // If this clique contains all the remaining variables then we can
//     // do all the remaining elimination
//     bool finished_elimination = clique_verts.size() > elim_order.size();
//     if(finished_elimination ) {
//       elim_order.clear();
//       // add all the remaining factors
//       foreach(vertex_id_t vid, clique_verts) { 
//         const mrf::vertex_data& vdata = mrf.vertex_data(vid);
//         // add all the unsued factor ids
//         foreach(size_t fid, vdata.factor_ids) {
//           if(used_factors.count(fid) == 0) { 
//             clique.factor_ids.insert(fid); 
//           }          
//         }
//         union_child_cliques.insert(child_cliques[vid].begin(), 
//                                    child_cliques[vid].end());
//       }
//     }

 
//     // Connect the clique to all of its children cliques
//     foreach(vertex_id_t child_id,  union_child_cliques) {
//       //      std::cout << child_id << " ";              
//       const junction_tree::vertex_data& child_clique = jt.vertex_data(child_id);
//       junction_tree::edge_data edata;
//       edata.variables = 
//         child_clique.variables.intersect(clique.variables);
//       // edata.variables = 
//       //   graphlab::set_intersect(child_clique.variables,
//       //                           clique.variables);
//       jt.add_edge(child_id, clique_id, edata);
//       jt.add_edge(clique_id, child_id, edata);
//     }

//     // std::cout << "( ";
//     // foreach(size_t fid, clique.factor_ids) std::cout << fid << " ";
//     // std::cout << ")";
//     //    std::cout << std::endl;

//     // If we still have variables to eliminate
//     if(!elim_order.empty()) {
//       // Disconnect variable and connect neighbors and mark their
//       // children cliques
//       vertex_set& elim_clique_set = child_cliques[elim_vertex];
//       foreach(const vertex_id_t n_vertex, clique_verts) {
//         if(n_vertex != elim_vertex) {
//           vertex_set& neighbor_set = neighbors[n_vertex];
//           // connect neighbors
//           neighbor_set.insert(clique_verts.begin(), clique_verts.end());
//           // disconnect the variable
//           neighbor_set.erase(elim_vertex);
//           // if the neighbor clique size becomes too large then we are done
//           size_t clique_size = neighbor_set.size();
//           // if(clique_size > MAX_DIM) return -1;
//           // Update the fill count
//           elim_order.update(n_vertex, in_tree.size() - clique_size);
//           // Update the clique neighbors
//           vertex_set& child_clique_set = child_cliques[n_vertex];
//           child_clique_set.insert(clique_id);
//           foreach(vertex_id_t vid, elim_clique_set) {
//             child_clique_set.erase(vid);
//           }
//         }
//       }
//     } // if the elim order is not empty
//   }
//   return tree_width;

// } // end of build junction tree




#include <graphlab/macros_undef.hpp>
#endif
