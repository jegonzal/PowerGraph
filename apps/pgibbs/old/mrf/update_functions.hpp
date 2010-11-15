#ifndef UPDATE_FUNCTIONS_HPP
#define UPDATE_FUNCTIONS_HPP

#include "data_structures.hpp"



// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>


/** Important!! edge_factor_id must be last */
enum constants {MAX_HEIGHT_ID, 
                NSAMPLES_ID,
                MAX_NSAMPLES_ID,
                USE_WEIGHTS_ID,
                PRUNING_ID,
                NROOTS_ID,
                EDGE_FACTOR_ID};



graph_type* global_graph;



const double grow_root_residual = 4.0;
const double grow_tree_residual = 1.0; /* between + [0,1) */
const double up_tree_residual = 2.0;
const double down_tree_residual = 3.0;



void grow_tree_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data);

void up_tree_update(gl::iscope& scope, 
                    gl::icallback& scheduler,
                    gl::ishared_data* shared_data);

void down_tree_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data);






//! Compute the edge weight
double compute_edge_weight(gl::iscope& scope,
                           gl::ishared_data* shared_data,
                           vertex_id_t parent,
                           edge_id_t parent_eid) {  
  assert(shared_data != NULL);
  const vertex_data& vdata = scope.vertex_data();
  const vertex_data& parent_vdata = scope.neighbor_vertex_data(parent);
  const edge_data& edata = scope.edge_data(parent_eid);
  const binary_factor& edge_factor =
    shared_data->get_constant(EDGE_FACTOR_ID + 
                              edata.factor_id).as<binary_factor>();
  // Construct the conditional
  unary_factor conditional(parent, parent_vdata.potential.arity());
  conditional.condition(edge_factor, vdata.asg);
  conditional.normalize();
  // Construct the message
  unary_factor marginal = vdata.potential;
  graphlab::edge_list in_edges = scope.in_edge_ids();
  foreach(edge_id_t eid, in_edges) {
    vertex_id_t neighbor = scope.source(eid); 
    if(neighbor != parent) {
      const vertex_data& ndata = scope.neighbor_vertex_data(neighbor);
      const edge_data& edata = scope.edge_data(eid);
      const binary_factor& edge_factor =
        shared_data->get_constant(EDGE_FACTOR_ID + 
                                  edata.factor_id).as<binary_factor>();
      marginal.condition(edge_factor, ndata.asg);
    }
  }
  marginal.normalize();
  unary_factor message(parent, parent_vdata.potential.arity());
  message.convolve(edge_factor, marginal);
  message.normalize();
  // // Compute the l1 difference
  // double residual = message.residual(conditional);
  double residual = 0;
  for(size_t i =0; i < message.arity(); ++i) {
    residual += std::abs(message.logP(i) - conditional.logP(i));    
  }
  residual /= message.arity();
  residual = std::tanh(residual);
  //std::cout << residual << std::endl;
  return residual;
}








vertex_id_t get_next_root(size_t current_root, 
                          size_t nroots,
                          size_t nvertices) {
  vertex_id_t block = current_root % nroots;
  vertex_id_t next_root = 
    block +  
    graphlab::random::rand_int(nvertices/nroots - 1) * nroots;
  // Next root is a valid vertex
  assert(next_root < nvertices);
  // Next root is in the same block
  assert(next_root % nroots == block);
  return next_root;
}









void grow_root_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data) {
  
  vertex_data& vdata = scope.vertex_data();  
  graphlab::edge_list out_edges = scope.out_edge_ids();

  // Test to see if this is a valid root
  bool root_candidate = vdata.state == vertex_data::AVAILABLE;

    


  // If this vertex is still a root candidate
  if(root_candidate) {
    // check the neighbors to see if they are in a tree
    foreach(edge_id_t eid, out_edges) {
      vertex_id_t neighbor = scope.target(eid);
      const vertex_data& neighbor_vdata =
        scope.neighbor_vertex_data(neighbor);
      if(!(neighbor_vdata.state == vertex_data::AVAILABLE ||
           neighbor_vdata.state == vertex_data::BOUNDARY)   ) {
        root_candidate = false; break;
      }
    } // End of foreach
  } // end of if

  // If it is still a root candidate then we have succeeded and
  // can grow a tree here
  if(root_candidate) {
    assert(vdata.state  == vertex_data::AVAILABLE);
    assert(vdata.height == 0);
    
    // Make this the root
    vdata.parent = scope.vertex();
    vdata.state  = vertex_data::TREE_NODE;
    const size_t& max_height =
      shared_data->get_constant(MAX_HEIGHT_ID).as<size_t>();
    vdata.height = max_height;
    

    if(vdata.height == 0) {
      // Actually start the upward pass
      up_tree_update(scope, scheduler, shared_data);
    } else { 
      // Count the number of unvisited neighbors
      size_t unvisited_neighbors = 0;
      // Add all the children
      foreach(edge_id_t eid, out_edges) {
        const vertex_id_t neighbor = scope.target(eid);
        const vertex_data& neighbor_vdata =
          scope.neighbor_vertex_data(neighbor);
        // If there are any unexplored neighbors
        if(neighbor_vdata.state == vertex_data::AVAILABLE) {
          edge_data& edata = scope.edge_data(eid);
          unvisited_neighbors++;
          const gl::update_task task(neighbor, grow_tree_update);
          const double residual = grow_root_residual;
          vdata.child_candidates.inc();
          edata.exploring = true;
          scheduler.add_task(task, residual);
        } else {
          // The neighbor must be a boundary neighbor
          assert(neighbor_vdata.state == vertex_data::BOUNDARY);
        }
      }
      // If all the neighbors are visited this is a leaf node and we 
      // start the upward pass immediately
      if( unvisited_neighbors == 0 ) {
        // Actually start the upward pass
        up_tree_update(scope, scheduler, shared_data);
      }
    }
  } else {
    // Try to start a root elsewhere
    const size_t& nroots = shared_data->get_constant(NROOTS_ID).as<size_t>();
    vertex_id_t next_root = get_next_root(scope.vertex(),
                                          nroots,
                                          scope.num_vertices());
    gl::update_task task(next_root, grow_root_update);
    double residual = grow_root_residual;
    scheduler.add_task(task, residual);
  }
} // End of grow root update




















void grow_tree_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data) { 
  assert(shared_data != NULL);

  // Get the vertex data
  // const vertex_id_t vid = scope.vertex();
  vertex_data& vdata = scope.vertex_data();

  // This vertex must be available to grow
  if(vdata.state != vertex_data::AVAILABLE) return;

  // Get neighbors
  graphlab::edge_list out_edges = scope.out_edge_ids();


  // Find a candidate parent and count the number of candidate parents
  vertex_id_t parent             = NULL_VID;
  edge_id_t   parent_edge        = NULL_VID;
  size_t      parent_height      = 0;
  bool        parent_found       = false;
  bool        is_tree_candidate  = true;

  // Ensure that there is only one path to the tree
  foreach(edge_id_t eid, out_edges) {
    vertex_id_t neighbor = scope.target(eid);
    const vertex_data& neighbor_vdata =
      scope.neighbor_vertex_data(neighbor);
    if(neighbor_vdata.state == vertex_data::TREE_NODE) {
      if(!parent_found){
        parent = neighbor;
        parent_edge = eid;
        parent_height = neighbor_vdata.height;
        parent_found = true;
      } else is_tree_candidate = false;  // Parent already found      
    } else if(! (neighbor_vdata.state == vertex_data::AVAILABLE ||
                 neighbor_vdata.state == vertex_data::BOUNDARY)  ) {
      is_tree_candidate = false;
    }
    if(!is_tree_candidate) break;
  }
  
  // If there are no parents then this must become a root. However only
  // grow_root is allowed to create a new roots so terminate
  if(is_tree_candidate && !parent_found ) return;

  // This vertex is a tree candidate if it has one parent and that
  // parent has a positive height
  is_tree_candidate =
    is_tree_candidate
    && parent_found
    && parent_height > 0;



  double pruning = shared_data->get_constant(PRUNING_ID).as<double>();
  if( pruning > 0 && is_tree_candidate ) {
    assert(parent != NULL_VID);
    assert(parent_edge != NULL_VID);
    assert(scope.target(parent_edge) == parent);
    // compute the edge weight
    double edge_weight =
      compute_edge_weight(scope, shared_data,
                          parent, parent_edge);
    // Get the edge to the parent
    edge_data& edata = scope.edge_data(parent_edge);
    edata.weight = edge_weight;
    // Parent remains parent depending on edge weight value
    //   double unif = gl::random::rand01();
    is_tree_candidate = pruning  < edge_weight;
    //    is_tree_candidate = edge_weight > 1.0E-3 ;
  }
 

  // If this is a tree vertex
  if( is_tree_candidate ) {
    assert(parent != NULL_VID);
    vdata.state  = vertex_data::TREE_NODE;
    vdata.parent = parent;
    vdata.height = parent_height - 1;

    
    if(vdata.height == 0) {
      // Actually start the upward pass
      up_tree_update(scope, scheduler, shared_data);
    } else {

      // Determine if edge weights should be computed
      const bool use_weights =
      shared_data->get_constant(USE_WEIGHTS_ID).as<bool>();
        
      // Count the number of unvisited neighbors
      size_t unvisited_neighbors = 0;
      
      // Add all the children
      foreach(edge_id_t eid, out_edges) {
        vertex_id_t neighbor = scope.target(eid);
        const vertex_data& neighbor_vdata =
          scope.neighbor_vertex_data(neighbor);
        edge_data& edata = scope.edge_data(eid);
                
        // If there are any unexplored neighbors
        if(neighbor_vdata.state == vertex_data::AVAILABLE) {
          unvisited_neighbors++;
          const gl::update_task task(neighbor, grow_tree_update);

          double residual = grow_tree_residual;
          if(use_weights) {
            graphlab::general_scope<graph_type>
              neighbor_scope(global_graph, neighbor, NULL);            
            // compute the edge weight
            double edge_weight =
              compute_edge_weight(neighbor_scope,
                                  shared_data,
                                  scope.vertex(),
                                  eid);
            residual += edge_weight;
            edge_data& edata = scope.edge_data(eid);
            edata.weight = residual;
          }
          edata.exploring = true;
          vdata.child_candidates.inc();
          scheduler.add_task(task, residual);
        }
      }
      // If all the neighbors are visited this is a leaf node and we 
      // start the upward pass immediately
      if( unvisited_neighbors == 0 ) {
        // Actually start the upward pass
        up_tree_update(scope, scheduler, shared_data);
      }
    }
  } else {
    assert(vdata.state == vertex_data::AVAILABLE);
    assert(vdata.parent == NULL_VID);
    assert(vdata.height == 0);

    // This vertex becomes a boundary vertex.
    vdata.state = vertex_data::BOUNDARY;
    graphlab::edge_list in_edges = scope.in_edge_ids();
    // Schedule all uncallibrated tree neighbors for upward pass 
    foreach(edge_id_t eid, in_edges) {
      vertex_id_t neighbor = scope.source(eid);
      vertex_data& neighbor_vdata =
        scope.neighbor_vertex_data(neighbor);
      edge_data& edata = scope.edge_data(eid);
      const bool old_exploring_value = edata.exploring; 
      if(edata.exploring) {
        assert(neighbor_vdata.child_candidates.value > 0);
        neighbor_vdata.child_candidates.dec();
        edata.exploring = false;
      }
      // If the neighbor is in an uncallibrated tree neighbor 
      if(neighbor_vdata.state == vertex_data::TREE_NODE &&
         old_exploring_value == true &&
         neighbor_vdata.child_candidates.value == 0) { 
        gl::update_task task(neighbor, up_tree_update);
        double residual = up_tree_residual;
        scheduler.add_task(task, residual);
      }
    }
  }
} // end of grow tree
















void up_tree_update(gl::iscope& scope, 
                    gl::icallback& scheduler,
                    gl::ishared_data* shared_data) {
  assert(shared_data != NULL);

  // Get the vertex data
  vertex_id_t vid = scope.vertex();
  vertex_data& vdata = scope.vertex_data();

  // If this vertex is not visited, is not a tree node, or has already
  // been calibrated simply return
  if( vdata.state != vertex_data::TREE_NODE ) return;

  assert(vdata.parent != NULL_VID);
  
  // Get neighbors
  graphlab::edge_list in_edges = scope.in_edge_ids();

  // See if this vertex is ready to be callibrated.  A vertex is ready
  // to be callibrated if all its neighbors are visited and all its
  // children are callibrated. or its height = 0
  if(vdata.height != 0) {
    foreach(edge_id_t eid, in_edges) {
      vertex_id_t neighbor = scope.source(eid);
      const vertex_data& neighbor_vdata =
        scope.neighbor_vertex_data(neighbor);
      const bool is_child = neighbor_vdata.parent == vid;
      // If the neighbors state has not yet been determined or the
      // neighbor is an uncalibrated child then we are not ready to go
      // up so return
      if( neighbor_vdata.state == vertex_data::AVAILABLE ||
          ( is_child &&
           neighbor_vdata.state != vertex_data::CALIBRATED) ) return;
    }
  }


    /* Check the state of the neighbors (either:)
     1) Child && Calibrated
     2) Parent && Tree Node
     3) Boundary
   */
  foreach(edge_id_t eid, in_edges) {
    const vertex_id_t neighbor = scope.source(eid);
    const vertex_data& neighbor_vdata = scope.neighbor_vertex_data(neighbor);
    const bool is_child  = neighbor_vdata.parent == vid;
    const bool is_parent = vdata.parent == neighbor;
    // is_child IMPLIES calibrated
    assert(!is_child || neighbor_vdata.state == vertex_data::CALIBRATED);
    // is_parent IMPLIES available
    assert(!is_parent || neighbor_vdata.state == vertex_data::TREE_NODE);
    // boundary or height = 0
    assert(is_child ||
           is_parent ||
           neighbor_vdata.state == vertex_data::BOUNDARY ||
           (neighbor_vdata.state == vertex_data::AVAILABLE && vdata.height == 0));    
    // is_parent IMPLIES available
    assert(!is_parent || neighbor_vdata.state == vertex_data::TREE_NODE);
  }

  // This vertex is ready for moving up Mark self as calibrated
  vdata.state = vertex_data::CALIBRATED;  




  vdata.bp_marginal = vdata.potential;
  // Collect all in messages from children and condition on
  // assignments to non-children neighbors EXCEPT the parent.
  foreach(edge_id_t eid, in_edges) {
    vertex_id_t  neighbor = scope.source(eid);
    // if the neighbor is not the parent
    if(neighbor != vdata.parent) {
      const edge_data&   edata  = scope.edge_data(eid);
      const vertex_data& ndata  = scope.neighbor_vertex_data(neighbor);
      // Determine if the neighbor is a child of this vertex
      if(vid == ndata.parent) {
        // Receive the message 
        vdata.bp_marginal.times( edata.message );
      } else {
        // Get the edge factor
        const binary_factor& edge_factor =
          shared_data->get_constant(EDGE_FACTOR_ID + edata.factor_id).as<binary_factor>();
        // standard gibbs conditional
        vdata.bp_marginal.condition(edge_factor, ndata.asg);
      }
    }
  }
  vdata.bp_marginal.normalize();    
  


  
  // If this is the root we begin the downward pass
  if(vdata.parent == vid) {
    // Start down pass on root
    down_tree_update(scope, scheduler, shared_data);
    // Respawn another root if necessary
    // Schedule the start of a new root
    const size_t& nroots = shared_data->get_constant(NROOTS_ID).as<size_t>();
    vertex_id_t next_root = get_next_root(scope.vertex(),
                                          nroots,
                                          scope.num_vertices());
    gl::update_task task(next_root, grow_root_update);
    double residual = grow_root_residual;
    scheduler.add_task(task, residual);   
  } else {   
    // Send the message upward to this vertex parents
    edge_data& edata = scope.edge_data(scope.edge(vid, vdata.parent));
 
    // Get the edge factor
    const binary_factor& edge_factor =
      shared_data->get_constant(EDGE_FACTOR_ID + edata.factor_id).as<binary_factor>();
    edata.message.convolve(edge_factor, vdata.bp_marginal);
    edata.message.normalize();

    vertex_data& parent_vdata = scope.neighbor_vertex_data(vdata.parent);
    const edge_id_t rev_parent_eid     = scope.edge(vdata.parent, vid);
    edge_data& rev_parent_edata  = scope.edge_data(rev_parent_eid);

    if(rev_parent_edata.exploring) {
      assert(parent_vdata.child_candidates.value > 0);
      parent_vdata.child_candidates.dec();
      rev_parent_edata.exploring = false;
    }

    if(parent_vdata.child_candidates.value == 0) {
      // Schedule the parent
      const gl::update_task task(vdata.parent, up_tree_update);
      const double residual = up_tree_residual;
      scheduler.add_task(task, residual);
    }
  }
}


























void down_tree_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data) {
  // Get the vertex id
  vertex_id_t vid = scope.vertex();
  vertex_data& vdata = scope.vertex_data();   

  // If the vertex is no longer visited then the down was spurious and
  // can be skipped?
  if( vdata.state == vertex_data::AVAILABLE) return;


  // If the vertex is not in the tree then we need to determine if it
  // needs to remain visited
  if( vdata.state == vertex_data::BOUNDARY ) {
    assert(vdata.height == 0);
    assert(vdata.parent == NULL_VID);
    // Count the tree (or calibrated) neighbors
    bool safe_to_remove = true;
    graphlab::edge_list in_edges = scope.in_edge_ids();
    foreach(edge_id_t eid, in_edges) {
      const vertex_id_t neighbor_vid = scope.source(eid);
      const vertex_data& ndata = scope.neighbor_vertex_data(neighbor_vid);
      if( ndata.state == vertex_data::TREE_NODE && ndata.height > 0) {
        safe_to_remove = false;
        break;
      }
    }
    if(safe_to_remove) {
      vdata.state = vertex_data::AVAILABLE;
    }
    return;
  }
    
  if (vdata.state != vertex_data::CALIBRATED) return;

  // This vertex must be calibrated
  assert(vdata.state == vertex_data::CALIBRATED);
  assert(vdata.parent != NULL_VID);

  
  // If a premature down is invoked on a branch of the tree then
  // ignore
  if(vdata.parent != vid) {
    const vertex_data& parent_vdata =
      scope.neighbor_vertex_data(vdata.parent);
    if(parent_vdata.parent != NULL_VID) return;
  }

  
  // Get the in and out edges for the vertex
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();    

  // Verify the tree invariant (Should be removed in final version)
  /* Check the state of the neighbors (either:)
     1) Child and calibrated
     2) Boundary
     3) Parent and available
   */
  foreach(edge_id_t eid, in_edges) {
    vertex_id_t neighbor = scope.source(eid);
    const vertex_data& neighbor_vdata =
      scope.neighbor_vertex_data(neighbor);
    const bool is_child  = neighbor_vdata.parent == vid;
    const bool is_parent = vdata.parent == neighbor;
    // is_child IMPLIES calibrated
    assert(!is_child || neighbor_vdata.state == vertex_data::CALIBRATED);
    // is_parent IMPLIES available
    assert(!is_parent || neighbor_vdata.parent == NULL_VID);
    // boundary or height = 0
    assert(is_child ||
           is_parent ||
           neighbor_vdata.state == vertex_data::BOUNDARY ||
           neighbor_vdata.state == vertex_data::AVAILABLE);
  }


  // If not root receive parent message
  if(vdata.parent != vid) {
    const edge_id_t parent_eid = scope.edge(vdata.parent, vid);
    const edge_data&   parent_edata   =
      scope.edge_data(parent_eid);
    vdata.bp_marginal.times( parent_edata.message );
  }
  
  // // Collect ALL in messages (and condition on all assignments to
  // // neighbors that are not in the tree)
  // unary_factor marginal = vdata.potential; 
 
  // foreach(edge_id_t eid, in_edges) {
  //   vertex_id_t  neighborid     = scope.source(eid);
  //   const vertex_data& ndata    = scope.neighbor_vertex_data(neighborid);
  //   const edge_data&   edata    = scope.edge_data(eid);
  //   const bool is_child = vid == ndata.parent;
  //   const bool is_parent = neighborid == vdata.parent;
  //   // Determine if the neighbor is either a
  //   //  1) child:    vid == ndata.parent        \ __ BP UPdate
  //   //  2) parent:   vdata.parent = neighborid  /
  //   //  3) External:  ----> Gibbs conditional
  //   if(is_child || is_parent) {
  //     // BP conditional
  //     marginal.times( edata.message );
  //   } else {
  //     // Get the edge factor
  //     const binary_factor& edge_factor =
  //       shared_data->get_constant(EDGE_FACTOR_ID + edata.factor_id).as<binary_factor>();
  //     // standard gibbs conditional
  //     marginal.condition(edge_factor, ndata.asg);
  //   }      
  // }
  // marginal.normalize();

  unary_factor cavity;  
  // Send messages to all children   
  foreach(edge_id_t eid, out_edges) {
    vertex_id_t  neighborid  = scope.target(eid);
    const vertex_data& ndata = scope.neighbor_vertex_data(neighborid);
    // If this edge is to a child in the tree
    if(vid == ndata.parent) {
      edge_data& edata       = scope.edge_data(eid);      
      // Construct the cavity
      cavity = vdata.bp_marginal;
      // Get the reverse edge data
      // OPTIMIZE
      edge_id_t rev_eid = scope.reverse_edge(eid);
      const edge_data& rev_edata = scope.edge_data(rev_eid);
      // Divide out the upward message
      cavity.divide(rev_edata.message);
      // Get the edge factor
      const binary_factor& edge_factor =
        shared_data->get_constant(EDGE_FACTOR_ID + edata.factor_id).as<binary_factor>(); 
      // Convolve the cavity with the edge factor     
      edata.message.convolve(edge_factor, cavity);
      edata.message.normalize();
    }
  }

  // Save the marginal as the new Rao-Blackwellized belief estimate
  vdata.belief.plus(vdata.bp_marginal);

  // If this is not the root vertex then when drawing a sample we
  // need to condition on the parents assignment so we must divide
  // out the bp message and replace it with the conditional
  if( vdata.parent != vid ) { // if not the root vertex
    // Get the bp message from the parent
    // OPTIMIZE
    const edge_data& edata = 
      scope.edge_data(scope.edge(vdata.parent, vid));
    // Divide out the message
    vdata.bp_marginal.divide(edata.message);
    // Get the edge factor
    const binary_factor& edge_factor =
      shared_data->get_constant(EDGE_FACTOR_ID + edata.factor_id).as<binary_factor>();
    // Multiply in the conditional
    const vertex_data& pdata = 
      scope.neighbor_vertex_data(vdata.parent);
    vdata.bp_marginal.condition(edge_factor, pdata.asg);
    vdata.bp_marginal.normalize();
  }
  
  // Now we have the correct marginal so sample from it:
  vdata.asg = vdata.bp_marginal.sample();
  vdata.updates++;

  // Invoke down tree on all children and boundary
  if(vdata.height > 0) {
    foreach(edge_id_t eid, out_edges) {    
      const vertex_id_t neighbor_vid    = scope.target(eid);
      const vertex_data& ndata = scope.neighbor_vertex_data(neighbor_vid);
      if( neighbor_vid != vdata.parent &&
          ndata.state != vertex_data::AVAILABLE ) {      
        gl::update_task task(neighbor_vid, down_tree_update);
        double residual = down_tree_residual;
        scheduler.add_task(task, residual);
      }
    }
  }

  // This vertex becomes available
  vdata.parent = NULL_VID;
  vdata.state  = vertex_data::AVAILABLE;
  vdata.height = 0;


} // end of down tree update






















 


template<bool UseCallback>
void single_sample_update(gl::iscope& scope, 
                          gl::icallback& scheduler,
                          gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  vertex_data& vdata = scope.vertex_data();
  unary_factor conditional = vdata.potential;
  foreach(edge_id_t eid, scope.in_edge_ids()) {
    vertex_id_t source = scope.source(eid);
    const edge_data& edata = scope.edge_data(eid);
    const vertex_data& neighbor_vdata =
      scope.neighbor_vertex_data(source);
    // Get the edge factor
    const binary_factor& edge_factor =
      shared_data->get_constant(EDGE_FACTOR_ID + edata.factor_id).as<binary_factor>();
    conditional.condition(edge_factor, neighbor_vdata.asg);
  }
  conditional.normalize();
  size_t sample = conditional.sample();
  vdata.asg = sample;
  vdata.updates++;
  vdata.belief.plus(conditional);

  // Reschedule self
  if(UseCallback) {
    gl::update_task task(scope.vertex(), single_sample_update<UseCallback> );
    double residual = 1.0;
    scheduler.add_task(task, residual);
  }
}












// /** Color selection function */
// bool select_color(uint32_t color,
//                   gl::vertex_id_t vertex,
//                   const vertex_data& vdata) {
//   return vdata.color == color;
// }


// // Global count of the number colors
// uint16_t total_colors = 0;

// /**
//  * The scheduling function used by the set scheduler to execute in
//  * colored order
//  */
// void planned_color_schedule(gl::set_scheduler &sched) {
//   // All sets must be created before scheduling calls
//   std::vector<gl::ivertex_set*> colorsets;
//   // Construct color sets
//   for (size_t i = 0; i < total_colors; ++i){
//     colorsets.push_back(&sched.attach(gl::rvset(gl::selector_function(boost::bind(select_color, i, _1, _2)), true), sched.root_set()));
//   }
  
//   gl::execution_plan eplan;

//   const bool use_callback = false;
//   gl::update_function update_function =
//     single_sample_update<use_callback>;

  
//   // Actually construct color sets
//   sched.init();
  
//   for (size_t i = 0; i < total_colors; ++i) {
//     std::cout <<"Set " << i <<" has " <<colorsets[i]->size() << " vertices\n";
//     eplan.execute(*colorsets[i], update_function);
//   }
  
//   graphlab::timer timer;
//   timer.start();
//   eplan.generate_plan(sched.get_graph(), sched.num_cpus());
//   std::cout << "Execution Plan took " << timer.current_time() << " to compile\n";

//   // Actually execut the schedule generating all the samples
//   while(!sched.completed()) { 
//     sched.execute_plan(eplan);
//   }

//   std::cout << "Finished" << std::endl;
// }





gl::iengine* make_colored_engine(gl::graph& graph,
                                 gl::thread_shared_data& sdm,
                                 const std::vector<binary_factor>& edge_factors,
                                 size_t ncpus) {

  gl::iengine* engine =
    graphlab::engine_factory::new_engine("async", "colored", "null",
                                         graph,
                                         ncpus);

  assert(engine != NULL);  
  // Set the shared data 
  engine->set_shared_data_manager(&sdm);
  std::cout << "Using set schedule planner." << std::endl;

  const bool use_callback = false;
  gl::update_function update_function =
    single_sample_update<use_callback>;

  engine->get_scheduler().set_option(gl::scheduler_options::UPDATE_FUNCTION, 
                                     (void*) update_function);
  // Add all the edge factors to the shared data
  for(size_t i = 0; i < edge_factors.size(); ++i) {
    sdm.set_constant(EDGE_FACTOR_ID + i, graphlab::any(edge_factors[i]) );
  }

  return engine;
} // end of make colored engine









gl::iengine* make_engine(gl::graph& graph,
                         gl::thread_shared_data& sdm,
                         const std::vector<binary_factor>& edge_factors,
                         size_t ncpus,
                         const std::string& scheduler_type) {
  // Set the global graph
  global_graph = &graph;

  // Create an engine
  gl::iengine* engine =
    graphlab::engine_factory::new_engine("async",
                                         scheduler_type,
                                         "edge",
                                         graph,
                                         ncpus);
  assert(engine != NULL);


  // Set the shared data manager
  engine->set_shared_data_manager(&sdm);


  if(scheduler_type == "clustered_priority") {
    size_t partition = graph.num_vertices() / 32;
    std::cout << "Partition: " << partition << std::endl;
    ((graphlab::clustered_priority_scheduler<graph_type>*)
     (&(engine->get_scheduler())))->cluster(partition, 
                                            graphlab::partition_method::PARTITION_METIS);
  }

  
  // Add all the edge factors to the shared data
  for(size_t i = 0; i < edge_factors.size(); ++i) {
    sdm.set_constant(EDGE_FACTOR_ID + i, graphlab::any(edge_factors[i]) );
  }
  return engine;
}










void set_tree_sampler_constants(gl::thread_shared_data& sdm,
                                size_t nroots,
                                size_t tree_height,
                                bool use_weights, 
                                double pruning) {
  // Create a shared data manager to pass the edge_factor
  sdm.set_constant(MAX_HEIGHT_ID,  graphlab::any(tree_height) );
  sdm.set_constant(USE_WEIGHTS_ID, graphlab::any(use_weights) );
  sdm.set_constant(PRUNING_ID, graphlab::any(pruning) );
  sdm.set_constant(NROOTS_ID, graphlab::any(nroots) );
}








/** Get the update counts for a vertex */
size_t get_nsamples(const vertex_data& vdata) {
  return vdata.updates;
}


bool nsamples_terminator(const gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  const size_t& max_nsamples =
    shared_data->get_constant(MAX_NSAMPLES_ID).as<size_t>();
  size_t nsamples = shared_data->get(NSAMPLES_ID).as<size_t>();
  bool terminate = nsamples >= max_nsamples;
  if(terminate) {
    //     std::cout << "Termination condition reached" << std::endl;
  }
  return terminate;
}








// struct pgibbs_results{
//   size_t actual_samples;
//   size_t update_counts;
//   size_t partition;
//   double runtime;
//   pgibbs_results() :
//     actual_samples(0),
//     update_counts(0),
//     partition(0),
//     runtime(0) { } 
// };                   










// pgibbs_results 
// parallel_tree_sample(gl::graph& graph,
//                      const std::vector<binary_factor>& edge_factors,
//                      size_t nsamples, 
//                      size_t ncpus,
//                      const std::string& scheduler_type,
//                      size_t nroots,
//                      size_t tree_height,
//                      bool use_weights, 
//                      double pruning) {

//   // Create a shared data manager to pass the edge_factor
//   gl::thread_shared_data sdm;
//   set_tree_sampler_constants(sdm,
//                              nroots,
//                              tree_height,
//                              use_weights,
//                              pruning);
//   gl::iengine* engine = make_engine(graph,
//                                     sdm,
//                                     edge_factors,
//                                     ncpus,
//                                     scheduler_type);
//   assert(engine != NULL);


  
//   // Set the terminator
//   sdm.set_constant(MAX_NSAMPLES_ID,  graphlab::any(nsamples) );
//   sdm.set_sync(NSAMPLES_ID,
//                gl::sync_ops::sum< size_t, get_nsamples >,
//                gl::apply_ops::identity_print< size_t >,
//                size_t(0),
//                1000);
//   engine->add_terminator(nsamples_terminator);
  

//   pgibbs_results result;
//   result.actual_samples = 0;

  
//   // Create 
//   for(size_t j = 0; j < nroots; ++j) {
//     // Add the root
//     vertex_id_t root = get_next_root(j, nroots, graph.num_vertices());
//     // Create a task to grow a tree from the root
//     double residual = 1.0;
//     gl::update_task task(root, grow_root_update);
//     engine->get_scheduler().add_task(task, residual);
//   }

//   // Run the engine
//   graphlab::timer time; time.start();  
//   engine->start();
//   result.runtime = time.current_time();

       
//   // Get the actual number of samples
//   sdm.sync(graph, NSAMPLES_ID);
//   result.actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();
//   result.update_counts = engine->last_update_count();

//   std::cout << "Update Count:    " << result.update_counts << std::endl
//             << "Runtime:         " << result.runtime << std::endl
//             << "Efficiency:      "
//             << (double(result.update_counts) / result.runtime)<< std::endl;

//   std::cout << "Samples collected: " 
//             << result.actual_samples
//             << std::endl;

//   // Ensure that enough samples are made
//   //  assert(result.actual_samples >= nsamples);
//   if(result.actual_samples < nsamples) {
//     std::cout << "Failed to generate sufficient samples"
//               << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
//               << std::endl;

//   }
//   // Delete the engine
//   delete engine;  
//   return result;
// }












// pgibbs_results 
// parallel_async_sample(gl::graph& graph,
//                       const std::vector<binary_factor>& edge_factors, 
//                       size_t nsamples, 
//                       size_t ncpus,
//                       const std::string& scheduler_type) {                    
//   // Create an engine
//   pgibbs_results result;


//   gl::thread_shared_data sdm;
//   gl::iengine* engine = make_engine(graph,
//                                     sdm,
//                                     edge_factors,
//                                     ncpus,
//                                     scheduler_type);
//   assert(engine != NULL); 
//   sdm.set_constant(MAX_NSAMPLES_ID,  graphlab::any(nsamples) );
//   sdm.set_sync(NSAMPLES_ID,
//                gl::sync_ops::sum< size_t, get_nsamples >,
//                gl::apply_ops::identity_print< size_t >,
//                size_t(0),
//                1000);
//   engine->add_terminator(nsamples_terminator);
  
//   // Set the shared data manager
//   engine->set_shared_data_manager(&sdm);

//   // Add the update function to all vertices
//   std::vector<vertex_id_t> vertices(graph.num_vertices(), 0);
//   for(vertex_id_t i = 0; i < vertices.size(); ++i) vertices[i] = i;
//   std::random_shuffle(vertices.begin(), vertices.end());
//   const double residual = 1.0;
//   const bool UseCallback = true;
//   for(size_t i = 0; i < vertices.size(); ++i) {
//     gl::update_task task(vertices[i], single_sample_update<UseCallback>);
//     engine->get_scheduler().add_task(task, residual);
//   }

//   // Run the engine
//   graphlab::timer time; time.start();
//   engine->start();
//   result.runtime = time.current_time();

//   // Track the number of samples that were actually taken
//   sdm.sync(graph, NSAMPLES_ID);
//   result.actual_samples = sdm.get(NSAMPLES_ID).as<size_t>();    

//   assert(result.actual_samples >= nsamples);
  



//   std::cout << "Samples collected: " 
//             << result.actual_samples
//             << std::endl;

//   delete engine;  
//   return result;

// }









#include <graphlab/macros_undef.hpp>
#endif
