#ifndef PGIBBS_SEQUENTIAL_JT_GIBBS_HPP
#define PGIBBS_SEQUENTIAL_JT_GIBBS_HPP


#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <boost/unordered_map.hpp>

#include "data_structures.hpp"
#include "update_function.hpp"
#include "image.hpp"


#include <graphlab/macros_def.hpp>
  


typedef graphlab::fast_set<16, vertex_id_t> vertex_set;
// typedef std::set<vertex_id_t> vertex_set;

typedef std::map<vertex_id_t, vertex_set> vset_map;
//typedef boost::unordered_map<vertex_id_t, vertex_set> vset_map;


struct elim_clique {
  vertex_id_t parent;
  vertex_id_t elim_vertex;
  vertex_set vertices; 
  vertex_set factor_ids;
  elim_clique() : parent(-1) { }
};


typedef std::vector<elim_clique> clique_vector;


template<typename A, typename B>
const B& safe_find(const std::map<A, B>& const_map, const A& key) {
  typedef typename std::map<A, B>::const_iterator const_iterator;
  const_iterator iter = const_map.find(key);
  assert(iter != const_map.end());
  return iter->second;
}




size_t build_minfill_elim_order(const vset_map& var2factors_const,
                                const vset_map& factor2vars_const,
                                const vertex_id_t max_factor_id,
                                std::vector<vertex_id_t>* elim_order = NULL,
                                clique_vector* cliques = NULL) {

  // Reset elim order and clique set (if they are available)
  if(elim_order != NULL) elim_order->clear();
  if(cliques != NULL) cliques->clear();

  // Make local copies of the maps
  vset_map var2factors(var2factors_const);
  vset_map factor2vars(factor2vars_const);

  // track the treewidth
  size_t tree_width = 0;

  // keep track of the next unique factor id value (used to created
  // temporary factors along the way).
  size_t next_new_factor_id = max_factor_id + 1;
  
  // Construct an elimination ordering:
  graphlab::mutable_queue<vertex_id_t, int> elim_priority_queue;
  typedef vset_map::value_type vset_map_pair;
  foreach(const vset_map_pair& pair, var2factors) {
    vertex_id_t vid = pair.first;
    const vertex_set& factors = pair.second;
    vertex_set fill_clique;
    foreach(vertex_id_t fid, factors) {
      const vertex_set& verts = factor2vars[fid];
      if(!verts.empty()) {
        fill_clique.insert(verts.begin(), verts.end());
      }
      assert(fid != next_new_factor_id);
    }
    int clique_edges = fill_clique.size() * (fill_clique.size() - 1) / 2;
      //  existing_edges;
    assert(clique_edges >= 0);
    elim_priority_queue.push(vid, -clique_edges);
  }



  
  // vertex_set used_factors;
  // std::map<vertex_id_t, clique_type> cliques;
  
  // Run the elimination;
  while(!elim_priority_queue.empty()) {
    const std::pair<vertex_id_t, int> top = elim_priority_queue.pop();
    const vertex_id_t elim_vertex = top.first;
    const vertex_set& factorset = var2factors[elim_vertex];
    
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
      assert(fid != next_new_factor_id);
    }

    // if necessary store the elimination ordering and the clique set
    if(elim_order != NULL) elim_order->push_back(elim_vertex);
    if(cliques != NULL) {
      elim_clique elim_clique;
      elim_clique.elim_vertex = elim_vertex;
      elim_clique.vertices = affected_vertices;
      cliques->push_back(elim_clique);
    }



    // Merge any factors
    if(factorset.size() > 1) {    
      // Build the new factor
      size_t new_factor_id = next_new_factor_id++;
      factor2vars[new_factor_id] = affected_vertices;
      // Remove all the factors from the affected vertices and add the
      // new factor
      foreach(vertex_id_t vid, affected_vertices) {
        var2factors[vid] -= factorset;
        var2factors[vid] += new_factor_id;
      }
    } // end of merge

    // Update the fill order
    foreach(vertex_id_t vid, affected_vertices) {
      const vertex_set& factors = var2factors[vid];
      vertex_set fill_clique;
      foreach(vertex_id_t fid, factors) {
        const vertex_set& verts = factor2vars[fid];
        if(!verts.empty()) {
          fill_clique.insert(verts.begin(), verts.end());
        }
      }
      int clique_edges = fill_clique.size() * (fill_clique.size() - 1) / 2;
        //   existing_edges;
      assert(clique_edges >= 0);
      elim_priority_queue.update(vid, -clique_edges);
    }

  }
  return tree_width;
} // end of build_minfill_elim_order









size_t eval_elim_order(const vset_map& var2factors_const,
                       const vset_map& factor2vars_const,
                       const vertex_id_t max_factor_id,
                       const std::vector<vertex_id_t>& elim_order,
                       clique_vector* cliques = NULL) {

  // Reset elim order and clique set (if they are available)
  if(cliques != NULL) cliques->clear();

  // Make local copies of the maps
  vset_map var2factors(var2factors_const);
  vset_map factor2vars(factor2vars_const);

  // track the treewidth
  size_t tree_width = 0;

  // keep track of the next unique factor id value (used to created
  // temporary factors along the way).
  size_t next_new_factor_id = max_factor_id + 1;

  // Run the elimination;
  rev_foreach(vertex_id_t elim_vertex, elim_order) {
    const vertex_set& factorset = var2factors[elim_vertex];
    
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
      assert(fid != next_new_factor_id);
    }

    // if necessary store the elimination ordering and the clique set
    if(cliques != NULL) {
      elim_clique elim_clique;
      elim_clique.elim_vertex = elim_vertex;
      elim_clique.vertices = affected_vertices;
      cliques->push_back(elim_clique);
    }


    // Merge any factors
    if(factorset.size() > 1) {
      // Build the new factor
      vertex_id_t new_factor_id = next_new_factor_id++;
      factor2vars[new_factor_id] = affected_vertices;
      // Remove all the factors from the affected vertices and 
      // add the new factor
      foreach(vertex_id_t vid, affected_vertices) {
        var2factors[vid] -= factorset;
        var2factors[vid] += new_factor_id;
      }
    } // end of merge
  }
  return tree_width;
} // end of eval elim order













template<typename T>
void junction_tree_from_cliques(const mrf::graph_type& mrf, 
                                const T& begin_iter,
                                const T& end_iter,
                                junction_tree::graph_type& jt)  {

  // Convert the iterators to a range and size
  const std::pair<T,T> cliques_range = 
    std::make_pair(begin_iter, end_iter);


  std::map<vertex_id_t, vertex_id_t> elim_time_map;
  std::set<vertex_id_t> assigned_factors;

  { // Compute the elimination time for each vertex that is eliminated
    size_t elim_time = 0;
    foreach(const elim_clique& clique, cliques_range) {
      elim_time_map[clique.elim_vertex] = elim_time++;
    }
  }


  { // Assign factors  
    rev_foreach(elim_clique& clique, cliques_range) {
      const mrf::vertex_data& vdata = mrf.vertex_data(clique.elim_vertex);
      foreach(vertex_id_t fid, vdata.factor_ids) {
        if(assigned_factors.count(fid) == 0) {
          clique.factor_ids += fid;
          assigned_factors.insert(fid);
        }
      }
    }
  }

 
  // Compute the parent of each clique
  foreach(elim_clique& clique, cliques_range) {
    vertex_id_t parent = 0;
    foreach(vertex_id_t vid, clique.vertices) {
      parent =  std::max(parent, elim_time_map[vid]);
    }
    clique.parent = parent;
  }





  // { // compute Subsumption
  //   begin_iter->parent = -1;
  //   rev_foreach(elim_clique& clique, cliques_range) {      
  //     for(vertex_id_t parent_id = clique.parent; 
  //         parent_id < (end_iter - begin_iter); 
  //         parent_id = clique.parent ) {         
  //       elim_clique& parent = *(begin_iter + parent_id);          
  //       if((parent.vertices + parent.elim_vertex) <= 
  //          (clique.vertices + clique.elim_vertex)) {
  //         clique.parent = parent.parent;
  //         parent.parent = parent_id;
  //         std::swap(parent, clique);
  //       } else break;
  //     }
  //   }
  // }




  {  // Construct the junction tree
    // Ensure that the parent of the root is identifiably undefined
    begin_iter->parent = -1; 
    foreach(elim_clique& clique, cliques_range) {      
      // Create the vertex data
      junction_tree::vertex_data vdata;
      // add the eliminated vertex
      vdata.variables = mrf.vertex_data(clique.elim_vertex).variable;
      foreach(vertex_id_t vid, clique.vertices) 
        vdata.variables += mrf.vertex_data(vid).variable;      
      foreach(vertex_id_t fid, clique.factor_ids) 
        vdata.factor_ids.insert(fid);
      // Add the vertex
      vertex_id_t child_id = jt.add_vertex(vdata);

      vertex_id_t parent_id = clique.parent;
      // Add the edge to parent if not root
      if(parent_id < (end_iter - begin_iter)) {
        // Get the parent vertex data
        const junction_tree::vertex_data& parent_vdata =
          jt.vertex_data(parent_id);
        junction_tree::edge_data edata;
        edata.variables = 
          vdata.variables.intersect(parent_vdata.variables);
        // Add the actual edges
        jt.add_edge(child_id, parent_id, edata);
        jt.add_edge(parent_id, child_id, edata);
      }
    } // end of for each
  } // End of construct cliques

  //   { // Print out the clique list
  //     size_t i = 0;
  //     foreach(const elim_clique& clique, cliques_range) {
  //       std::cout << i << " --> " << clique.parent << "   \t" 
  //                 << clique.elim_vertex << " : "
  //                 << (clique.vertices + clique.elim_vertex) 
  //                 << "    Factors"  << clique.factor_ids
  //                 <<std::endl;
  //       i++;
  //     }
  //   }
} // end of build junction tree









/**
 *  Extend the clique tree with the next vertex
 *
 **/
bool extend_clique_list(const mrf::graph_type& mrf,
                        vertex_id_t elim_vertex,
                        std::map<vertex_id_t, vertex_id_t>& elim_time_map,
                        clique_vector& cliques) {
  // sanity check: The vertex to eliminate should not have already
  // been eliminated
  assert(elim_time_map.find(elim_vertex) == elim_time_map.end());

  // Construct the elimination clique for the new vertex
  elim_clique clique;
  clique.elim_vertex = elim_vertex;
  foreach(edge_id_t ineid, mrf.in_edge_ids(elim_vertex)) {
    vertex_id_t vid = mrf.source(ineid);
    // if the neighbor is in the set of vertices being eliminated
    if(elim_time_map.find(vid) != elim_time_map.end()) {
      clique.vertices += mrf.source(ineid);
    }
    // if the clique ever gets too large than teminate
    if(clique.vertices.size() + 1 > MAX_DIM) return false;
  }

  // Determine the parent of this clique
  vertex_id_t parent_id = 0;
  foreach(vertex_id_t vid, clique.vertices)
    parent_id = std::max(parent_id, elim_time_map[vid]);
  clique.parent = parent_id;


  // Simulate injecting vertices in parent cliques back to when RIP is
  // satisfied
  vertex_set rip_verts = clique.vertices;
  for(vertex_id_t parent_vid = clique.parent; 
      !rip_verts.empty() && parent_vid < cliques.size(); ) {
    const elim_clique& parent_clique = cliques[parent_vid];    

    // otherwise update that the rip_verts
    rip_verts += parent_clique.vertices;
    rip_verts -= parent_clique.elim_vertex;

    // Check that the expanded clique is still within tree width
    vertex_set tmp_vset = rip_verts + parent_clique.elim_vertex;
    if(tmp_vset.size() > MAX_DIM) return false;

    // Find the new parent
    vertex_id_t new_parent_vid = 0;
    foreach(vertex_id_t vid, rip_verts) {
      new_parent_vid = 
        std::max(new_parent_vid, elim_time_map[vid]);
    }
    parent_vid = new_parent_vid;
  }

  // Assert that if we reached this point RIP can be satisfied safely
  // so proceed to update local data structures
  size_t elim_time = cliques.size();
  cliques.push_back(clique);
  elim_time_map[clique.elim_vertex] = elim_time;

  // Satisfy RIP
  rip_verts = clique.vertices;
  for(vertex_id_t parent_vid = clique.parent; 
      !rip_verts.empty() && parent_vid < cliques.size(); ) {
    // get the parent clique
    elim_clique& parent_clique = cliques[parent_vid];       

    // otherwise update that the rip_verts
    rip_verts += parent_clique.vertices;
    rip_verts -= parent_clique.elim_vertex;

    // Update the clique
    parent_clique.vertices = rip_verts;

    // Determine the new parent (except first vertex)
    vertex_id_t new_parent_vid = 0;
    foreach(vertex_id_t vid, parent_clique.vertices) {
      new_parent_vid = 
        std::max(new_parent_vid, elim_time_map[vid]);
    }
    parent_vid = new_parent_vid;

    // Update the parent for this clique
    parent_clique.parent = new_parent_vid;
  }
  // Add successfully
  return true;
}











///// Old code ----------------------------------------------
/////////////////////////////////////////////////



















size_t incremental_build_junction_tree(const mrf::graph_type& mrf,
                                       vertex_id_t root,
                                       junction_tree::graph_type& jt) {
  jt.clear();



  // Local data structures
  std::map<vertex_id_t, vertex_id_t> elim_time_map;
  clique_vector cliques;

  std::vector<vertex_id_t> elim_order;

  std::queue<vertex_id_t> bfs_queue;
  std::set<vertex_id_t> visited;

  // add the root
  bfs_queue.push(root);
  visited.insert(root);



  while(!bfs_queue.empty()) {
    // Take the top element
    const vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop(); 

    // test the 
    bool safe_extension = 
      extend_clique_list(mrf, next_vertex,
                         elim_time_map,
                         cliques);
    if(safe_extension) {   
      // Save the elimited vertex
      elim_order.push_back(next_vertex);
      // add the neighbors to the search queue
      foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
        vertex_id_t neighbor_vid = mrf.target(eid);
        if(visited.count(neighbor_vid) == 0) {
          bfs_queue.push(neighbor_vid);
          visited.insert(neighbor_vid);
        }
      }
    }

//     std::cout << "========================================="
//               << std::endl
//               << next_vertex << ":  " << safe_extension
//               << std::endl;    
//     { // Print out the clique list
//       size_t i = 0;
//       foreach(const elim_clique& clique, cliques) {
//         std::cout << i << " --> " << clique.parent << "   \t" 
//                   << clique.elim_vertex << " : "
//                   << (clique.vertices + clique.elim_vertex) 
//                   << "    Factors"  << clique.factor_ids
//                   << std::endl;
//         i++;
//       }

//       if(elim_order.size() == 2000) break;
  } // end of while loop
  



  std::cout << "Varcount: " << elim_order.size() << std::endl;  
  
  junction_tree_from_cliques(mrf, 
                             cliques.begin(), cliques.end(), 
                             jt);

 
  // image img(200, 200);
  // size_t index = elim_order.size();
  // foreach(vertex_id_t vid, elim_order) {
  //   img.pixel(vid) = index++;
  // }
  // img.save("tree.pgm");

  std::cout << "Finished build junction tree" << std::endl;
  return 1;
}






















size_t bfs_build_junction_tree(const mrf::graph_type& mrf,
                               vertex_id_t root,
                               junction_tree::graph_type& jt) {
  jt.clear();
  vset_map var2factors;
  vset_map factor2vars;


  std::queue<vertex_id_t> bfs_queue;
  std::set<vertex_id_t> visited;

  std::vector<vertex_id_t> elim_order;
  size_t tree_width = 0;

  // add the root
  bfs_queue.push(root);
  visited.insert(root);
  vertex_id_t max_factor_id = 0;
  while(!bfs_queue.empty()) {
    // Take the top element
    const vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop(); 

    // Update data structures 
    const mrf::vertex_data& vdata = mrf.vertex_data(next_vertex);
    var2factors[next_vertex] = vdata.factor_ids;
    vertex_id_t tmp_max_factor_id = max_factor_id;
    foreach(vertex_id_t fid, vdata.factor_ids) {
      factor2vars[fid].insert(next_vertex); 
      tmp_max_factor_id = std::max(tmp_max_factor_id, fid);
    }

    // build a junction tree using min fill
    elim_order.push_back(next_vertex);
    tree_width = eval_elim_order(var2factors, factor2vars,
                                 tmp_max_factor_id, elim_order);
    
    if(tree_width <= MAX_DIM) {   
      // add the neighbors to the search queue
      foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
        vertex_id_t neighbor_vid = mrf.target(eid);
        if(visited.count(neighbor_vid) == 0) {
          bfs_queue.push(neighbor_vid);
          visited.insert(neighbor_vid);
        }
      }
      // keep the current max factor id
      max_factor_id = tmp_max_factor_id;
    } else {
      // remove the variable if we decide not to use it
      const mrf::vertex_data& vdata = mrf.vertex_data(next_vertex);
      var2factors.erase(next_vertex);
      foreach(vertex_id_t fid, vdata.factor_ids) {
        factor2vars[fid].erase(next_vertex); 
      } 
      elim_order.pop_back();
    }
    if(var2factors.size() == 2000) break;
  }
  
  clique_vector cliques;
  tree_width = eval_elim_order(var2factors, factor2vars, 
                               max_factor_id + 1, 
                               elim_order,
                               &cliques);


  std::cout << "Varcount: " << var2factors.size() << std::endl;
  junction_tree_from_cliques(mrf, 
                             cliques.rbegin(), cliques.rend(),
                             jt);
  
//   std::cout << "Elim Tree Width: " << tree_width << std::endl;
//   image img(200, 200);
//   size_t index = elim_order.size();
//   foreach(vertex_id_t vid, elim_order) {
//     img.pixel(vid) = index++;
//   }
//   img.save("tree.pgm");

  return tree_width;
}







size_t min_fill_build_junction_tree(const mrf::graph_type& mrf,
                                    vertex_id_t root,
                                    junction_tree::graph_type& jt) {
  jt.clear();

  vset_map var2factors;
  vset_map factor2vars;


  std::queue<vertex_id_t> bfs_queue;
  std::set<vertex_id_t> visited;

  std::vector<vertex_id_t> elim_order;
  size_t tree_width = 0;

  // add the root
  bfs_queue.push(root);
  visited.insert(root);
  vertex_id_t max_factor_id = 0;
  while(!bfs_queue.empty()) {
    // Take the top element
    const vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop(); 

    // Update data structures 
    const mrf::vertex_data& vdata = mrf.vertex_data(next_vertex);
    var2factors[next_vertex] = vdata.factor_ids;
    vertex_id_t tmp_max_factor_id = max_factor_id;
    foreach(vertex_id_t fid, vdata.factor_ids) {
      factor2vars[fid].insert(next_vertex); 
      tmp_max_factor_id = std::max(tmp_max_factor_id, fid);
    }

    // build a junction tree using min fill
    tree_width = build_minfill_elim_order(var2factors, factor2vars, 
                                          tmp_max_factor_id);

    // if the treewidth is below the max 
    if(tree_width <= MAX_DIM) {
      // add the neighbors to the search queue
      foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
        vertex_id_t neighbor_vid = mrf.target(eid);
        if(visited.count(neighbor_vid) == 0) {
          bfs_queue.push(neighbor_vid);
          visited.insert(neighbor_vid);
        }
      }
      // keep the current max factor id
      max_factor_id = tmp_max_factor_id;
    } else {
      // remove the variable if we decide not to use it
      const mrf::vertex_data& vdata = mrf.vertex_data(next_vertex);
      var2factors.erase(next_vertex);
      foreach(vertex_id_t fid, vdata.factor_ids) {
        factor2vars[fid].erase(next_vertex); 
      } 
    }

    if(var2factors.size() == 2000) break;

  }


  
  clique_vector cliques;
  tree_width = 
    build_minfill_elim_order(var2factors, factor2vars, 
                             max_factor_id, &elim_order,
                             &cliques);


  std::cout << "Varcount: " << var2factors.size() << std::endl;
  junction_tree_from_cliques(mrf, 
                             cliques.rbegin(), cliques.rend(),
                             jt);

//   std::cout << "Min Fill Tree Width: " << tree_width << std::endl;


//   image img(200, 200);
//   size_t index = elim_order.size();
//   foreach(vertex_id_t vid, elim_order) {
//     img.pixel(vid) = index++;
//   }
//   img.save("tree.pgm");

  

  return tree_width;
}








  
void sample_once(const factorized_model& factor_graph,
                 mrf::graph_type& mrf,
                 vertex_id_t root) {

  junction_tree::gl::core jt_core;
  std::cout << "Building Tree" << std::endl;
  size_t tree_width = incremental_build_junction_tree(mrf, root, jt_core.graph());
  std::cout << "Root:  " << root << " ----------------" << std::endl;
  std::cout << "Tree width: " << tree_width << std::endl;

  std::cout << "Done!!!" << std::endl;


  // Setup the core
  jt_core.set_scheduler_type("fifo");
  jt_core.set_scope_type("edge");
  jt_core.set_ncpus(2);
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
