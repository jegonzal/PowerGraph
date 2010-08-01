#ifndef SEQUENTIAL_TREE_GIBBS
#define SEQUENTIAL_TREE_GIBBS


#include <cassert>
#include <vector>
#include <queue>
#include <stack>
#include <string>
#include <set> 
#include <algorithm>
#include <limits>
#include <cmath>


#include "data_structures.hpp"


// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>

// ===============================================================
// =============================================================
// Tree growing  =============================================
// =======================================================


/** ******************************************************
 * Grow comb tree
 */
void grow_comb_tree(gl::graph& graph,
                    bool first_comb,
                    size_t rows, size_t cols,
                    std::vector<vertex_id_t>& tree_order) {
  tree_order.clear();
  if(first_comb) {
    // left column and first row
    // Rxxxxxxxx
    // x
    // xxxxxxxxx
    // x
    
    // add root
    tree_order.push_back(image::vertid(rows, cols, 0, 0));
    // Add left column
    for(size_t r = 1; r < rows; ++r) {
      vertex_id_t vid = image::vertid(rows, cols, r, 0);
      tree_order.push_back(vid);
      vertex_id_t parent_vid = image::vertid(rows, cols, r-1, 0);
      graph.vertex_data(vid).parent = parent_vid;      
    }
    // Add all the rows
    for(size_t r = 0; r < rows; r += 2) {
      for(size_t c = 1; c < cols-1; ++c) {
        vertex_id_t vid = image::vertid(rows, cols, r, c);
        tree_order.push_back(vid);
        vertex_id_t parent_vid = image::vertid(rows, cols, r, c-1);
        graph.vertex_data(vid).parent = parent_vid;      
      }
    }    
  } else {
    // right column and second row
    //         R
    // xxxxxxxxx
    //         x
    // xxxxxxxxx

    // add root
    tree_order.push_back(image::vertid(rows, cols, 0, cols-1));
    // Add left column
    for(size_t r = 1; r < rows; ++r) {
      vertex_id_t vid = image::vertid(rows, cols, r, cols-1);
      tree_order.push_back(vid);
      vertex_id_t parent_vid = image::vertid(rows, cols, r-1, cols-1);
      graph.vertex_data(vid).parent = parent_vid;      
    }
    // Add all the rows
    for(size_t r = 1; r < rows; r += 2) {
      for(size_t c = cols-2; c > 0 ; --c) {
        vertex_id_t vid = image::vertid(rows, cols, r, c);
        tree_order.push_back(vid);
        vertex_id_t parent_vid = image::vertid(rows, cols, r, c+1);
        graph.vertex_data(vid).parent = parent_vid;      
      }
    }    
  }
} // end of comb_tree




/** ******************************************************
 * Grow breadth first search tree
 */
void grow_bfs_tree(gl::graph& graph,
                   vertex_id_t root,
                   size_t max_tree_size,
                   std::vector<vertex_id_t>& tree_order) {
  std::queue<vertex_id_t> bfs_queue;
  std::set<vertex_id_t> in_tree, visited;

  bfs_queue.push(root);
  visited.insert(root);
  
  while(!bfs_queue.empty()) {
    vertex_id_t v = bfs_queue.front(); 
    bfs_queue.pop();

    // Get the potential parents in the tree
    size_t parents = 0;
    vertex_id_t last_parent = 0;
    foreach(edge_id_t eid, graph.out_edge_ids(v) ) {
      vertex_id_t neighbor = graph.target(eid);
      if(in_tree.count(neighbor) > 0) {
        last_parent = neighbor;
        parents++;
      }
    }
    
    // Safely added if at one most possible parrent
    if(parents == 1 || v == root) {
      tree_order.push_back(v);
      in_tree.insert(v);
      
      // Set parent
      if(v != root) {
        graph.vertex_data(v).parent = last_parent;
      }

      // Stop growing if necessary
      if(tree_order.size() == max_tree_size) return;
      
      // Add any neighbors as potential candidates
      graphlab::edge_list out_edges = graph.out_edge_ids(v);
      std::vector<vertex_id_t> neighbors(out_edges.size());
      
      // Add the neighbors in a random order
      for(size_t i = 0; i < out_edges.size(); ++i) 
        neighbors[i] = graph.target(out_edges[i]);
      
      // Shuffle the neighbors to generate different trees each time
      std::random_shuffle(neighbors.begin(), neighbors.end());
      
      // Add the viable neighbors
      foreach(vertex_id_t neighbor, neighbors) {
        //std::cout << neighbor;
        if( visited.count(neighbor) == 0 ) {
          visited.insert(neighbor);
          bfs_queue.push(neighbor);
        }
      }
    } 
  }
} // End of grow_bfs_tree



/** ******************************************************
 * Grow weighted bfs tree
 */
void grow_weighted_bfs_tree(gl::graph& graph,
                            vertex_id_t root,
                            size_t max_tree_size,
                            std::vector<vertex_id_t>& tree_order) {
  typedef std::pair<float, vertex_id_t> elem_type;
  std::priority_queue<elem_type> bfs_queue;
  std::set<vertex_id_t> in_tree, visited;

  double root_priority = -float(graph.vertex_data(root).updates);
  bfs_queue.push( std::make_pair(root_priority, root) );
  visited.insert(root);
  
  while(!bfs_queue.empty()) {
    vertex_id_t v = bfs_queue.top().second; 
    bfs_queue.pop();

    // Get the potential parents in the tree
    size_t parents = 0;
    vertex_id_t last_parent = 0;
    foreach(edge_id_t eid, graph.out_edge_ids(v) ) {
      vertex_id_t neighbor = graph.target(eid);
      if(in_tree.count(neighbor) > 0) {
        last_parent = neighbor;
        parents++;
      }
    }
    
    // Safely added if at one most possible parrent
    if(parents == 1 || v == root) {
      tree_order.push_back(v);
      in_tree.insert(v);
      
      // Set parent
      if(v != root) {
        graph.vertex_data(v).parent = last_parent;
      }

      // Stop growing if necessary
      if(tree_order.size() == max_tree_size) return;

      
      // Add the viable neighbors
      graphlab::edge_list out_edges = graph.out_edge_ids(v);
      foreach(edge_id_t eid, out_edges) {
        vertex_id_t neighbor = graph.target(eid);
        // If the neighbor has not been evaluated yet
        if( visited.count(neighbor) == 0 ) {
          double priority = -float(graph.vertex_data(neighbor).updates);
          visited.insert(neighbor);
          bfs_queue.push( std::make_pair(priority, neighbor) );
        }
      }
    } // end of if safe to add 
  }
} // End of grow_weighted_bfs_tree





/** ******************************************************
 * Grow depth first tree
 */
void grow_dfs_tree(gl::graph& graph,
                   vertex_id_t root,
                   size_t max_tree_size,
                   std::vector<vertex_id_t>& tree_order) {
  std::stack<vertex_id_t> dfs_queue;
  std::set<vertex_id_t> in_tree, visited;

  dfs_queue.push(root);
  visited.insert(root);
  
  while(!dfs_queue.empty()) {
    vertex_id_t v = dfs_queue.top(); 
    dfs_queue.pop();

    // Get the potential parents in the tree
    size_t parents = 0;
    vertex_id_t last_parent = 0;
    foreach(edge_id_t eid, graph.out_edge_ids(v) ) {
      vertex_id_t neighbor = graph.target(eid);
      if(in_tree.count(neighbor) > 0) {
        last_parent = neighbor;
        parents++;
      }
    }
    
    // Safely added if at one most possible parrent
    if(parents == 1 || v == root) {
      tree_order.push_back(v);
      in_tree.insert(v);
      
      // Set parent
      if(v != root) {
        graph.vertex_data(v).parent = last_parent;
      }

      // Stop growing if necessary
      if(tree_order.size() == max_tree_size) return;
      
      // Add any neighbors as potential candidates
      graphlab::edge_list out_edges = graph.out_edge_ids(v);
      std::vector<vertex_id_t> neighbors(out_edges.size());
      
      // Add the neighbors in a random order
      for(size_t i = 0; i < out_edges.size(); ++i) 
        neighbors[i] = graph.target(out_edges[i]);
      
      // Shuffle the neighbors to generate different trees each time
      std::random_shuffle(neighbors.begin(), neighbors.end());
      
      // Add the viable neighbors
      foreach(vertex_id_t neighbor, neighbors) {
        //std::cout << neighbor;
        if( visited.count(neighbor) == 0 ) {
          visited.insert(neighbor);
          dfs_queue.push(neighbor);
        }
      }
    } 
  }
} // end of graph_dfs_tree




// ===============================================================
// =============================================================
// Sampling =================================================
// =======================================================


void tree_sample(gl::graph& graph,
                 const binary_factor& edge_factor,
                 std::vector<vertex_id_t>& tree_order) {
  graphlab::unary_factor marginal, cavity;

  // Upward sweep to the root (excluding root)
  for(size_t i = tree_order.size() - 1; i > 0; --i) {
    vertex_id_t vid = tree_order[i];
    const vertex_data& vdata = graph.vertex_data(vid);
    // Assert that the vertex has a parent
    assert(vdata.parent < graph.num_vertices());
    // construct the marginal (excluding the parent)
    marginal = vdata.potential;
    // Collect all in messages from children and condition on
    // assignments to non-children neighbors EXCEPT the parent.
    foreach(edge_id_t eid, graph.in_edge_ids(vid)) {
      vertex_id_t  neighbor = graph.source(eid);
      // if the neighbor is not the parent
      if(neighbor != vdata.parent) {
        const edge_data&   edata  = graph.edge_data(eid);
        const vertex_data& ndata  = graph.vertex_data(neighbor);
        // Determine if the neighbor is a child of this vertex
        if(vid == ndata.parent) {
          // Receive the message 
          marginal.times( edata.message );
        } else {
          // standard gibbs conditional
          marginal.condition(edge_factor, ndata.asg);
        }
      }
    }
    marginal.normalize();    
    // Send the message upward to this vertex parents
    edge_data& edata = graph.edge_data(vid, vdata.parent);
    edata.message.convolve(edge_factor, marginal);
    edata.message.normalize();
  } // end of upward sweep


  // Downward sweep starting at root (Sampling is done here as well)
  for(size_t i = 0; i < tree_order.size(); ++i) {
    vertex_id_t vid = tree_order[i];
    vertex_data& vdata = graph.vertex_data(vid);    

    // Get the in and out edges for the vertex
    graphlab::edge_list in_edges = graph.in_edge_ids(vid);
    graphlab::edge_list out_edges = graph.out_edge_ids(vid);    
    
    // Collect ALL in messages (and condition on all assignments to
    // neighbors that are not in the tree)
    marginal = vdata.potential;   
    foreach(edge_id_t eid, in_edges) {
      vertex_id_t  neighborid            = graph.source(eid);
      const vertex_data& ndata    = graph.vertex_data(neighborid);
      const edge_data&   edata    = graph.edge_data(eid);
      // Determine if the neighbor is either a
      //  1) child:    vid == ndata.parent        \ __ BP UPdate
      //  2) parent:   vdata.parent = neighborid  /
      //  3) External:  ----> Gibbs conditional
      if(vid == ndata.parent || neighborid == vdata.parent) {
        // BP conditional
        marginal.times( edata.message );
      } else {
        // standard gibbs conditional
        marginal.condition(edge_factor, ndata.asg);
      }      
    }
    marginal.normalize();

    // Send messages to all children   
    foreach(edge_id_t eid, out_edges) {
      vertex_id_t  neighborid  = graph.target(eid);
      const vertex_data& ndata = graph.vertex_data(neighborid);
      // If this edge is to a child in the tree
      if(vid == ndata.parent) {
        edge_data& edata       = graph.edge_data(eid);      
        // Construct the cavity
        cavity = marginal;
        // Get the reverse edge data
        edge_id_t rev_eid = graph.rev_edge_id(eid);
        const edge_data& rev_edata = graph.edge_data(rev_eid);
        // Divide out the upward message
        cavity.divide(rev_edata.message);
        // Convolve the cavity with the edge factor
        edata.message.convolve(edge_factor, cavity);
        edata.message.normalize();
      }
    }

    // Save the marginal as the new Rao-Blackwellized belief estimate
    vdata.belief.plus(marginal);
    // vdata.belief = marginal;

    // If this is not the root vertex then when drawing a sample we
    // need to condition on the parents assignment so we must divide
    // out the bp message and replace it with the conditional
    if( i != 0 ) { // if not the root vertex
      assert(vdata.parent < graph.num_vertices());
      // Get the bp message from the parent
      const edge_data& edata = graph.edge_data(vdata.parent, vid);
      // Divide out the message
      marginal.divide(edata.message);
      // Multiply in the conditional
      const vertex_data& pdata = graph.vertex_data(vdata.parent);
      marginal.condition(edge_factor, pdata.asg);
      marginal.normalize();
    }
    // Make this vertex the root
    vdata.parent = -1;

    
    // Now we have the correct marginal so sample from it:
    vdata.asg = marginal.sample();
    vdata.updates++;


  } // end of downward sweep
} // end of tree_bp






void tree_sample(gl::graph& graph,
                 const binary_factor& edge_factor,
                 size_t max_tree_size,
                 size_t max_samples) {
  size_t nsamples = 0;
  std::vector<vertex_id_t> tree_order;
  while(nsamples < max_samples) {
    // select a root
    vertex_id_t root =
      graphlab::random::rand_int(graph.num_vertices() - 1);
    // Build a tree order
    tree_order.clear();
    size_t next_tree_size = std::min(max_tree_size,
                                     max_samples - nsamples);
    // grow_bfs_tree(graph, root, next_tree_size, tree_order);
    //    grow_dfs_tree(graph, root, next_tree_size, tree_order);
    grow_weighted_bfs_tree(graph, root, next_tree_size, tree_order);
    nsamples += tree_order.size();    
    tree_sample(graph, edge_factor, tree_order);
  }
  // we should have exactly max_samples when done
  assert(nsamples == max_samples);
}




void comb_sample(gl::graph& graph,
                 size_t rows, size_t cols,
                 const binary_factor& edge_factor) {
  std::vector<vertex_id_t> tree_order;
  
  grow_comb_tree(graph, true, rows, cols, tree_order);
  tree_sample(graph, edge_factor, tree_order);

  tree_order.clear();
  
  grow_comb_tree(graph, false, rows, cols, tree_order);
  tree_sample(graph, edge_factor, tree_order);

}










#include <graphlab/macros_undef.hpp>
#endif




