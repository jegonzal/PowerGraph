#ifndef UPDATE_FUNCTIONS_HPP
#define UPDATE_FUNCTIONS_HPP

#include <boost/unordered_map.hpp>

#include "data_structures.hpp"

// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>


typedef std::pair<variable_t, uint32_t> mini_asg_t;

typedef std::map<vertex_id_t, mini_asg_t> asg_map_type;

/** Important!! edge_factor_id must be last */
enum constants {MAX_HEIGHT_ID, 
                NSAMPLES_ID,
                MAX_NSAMPLES_ID,
                USE_WEIGHTS_ID,
                PRUNING_ID,
                NROOTS_ID,
                FACTOR_ID};



graph_type* global_graph;


const double grow_root_residual = 4.0;
const double grow_tree_residual = 1.0; /* between + [0,1) */
const double up_tree_residual   = 2.0;
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



//! Compute the unormalized likelihood of the current assignment
double unnormalized_likelihood(const graph_type& graph,
                               const std::vector<factor_t>& factors) {
  double sum = 0;
  // Sum the logprob of each factor
  foreach(const factor_t& factor, factors) {
    // Accumulate the assignments 
    domain_t dom = factor.args();
    assignment_t asg;
    for(size_t i = 0; i < dom.num_vars(); ++i) {
      const vertex_id_t vid = dom.var(i).id;
      const vertex_data& vdata = graph.vertex_data(vid);
      assert(vdata.variable == dom.var(i));
      asg &= assignment_t(vdata.variable, vdata.asg);
    }
    sum += factor.logP(asg);
  }
  return sum;
}



//! Compute the edge weight
double compute_edge_weight(vertex_id_t parent,
                           vertex_id_t child,
                           gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  assert(global_graph != NULL);
  const vertex_data& parent_vdata = global_graph->vertex_data(parent);
  const vertex_data& child_vdata = global_graph->vertex_data(child);
  graphlab::edge_list child_in_edges = global_graph->in_edge_ids(child);
  
  // Collect assignment to all variables in the neighborhood of parent and child
  typedef asg_map_type::const_iterator map_iterator;
  asg_map_type large_assignment;

  // loop over edges collecting assignments
  foreach(edge_id_t eid, child_in_edges) {
    const vertex_id_t neighbor = global_graph->source(eid);
    const vertex_data& neighbor_vdata = global_graph->vertex_data(neighbor);
    large_assignment[neighbor_vdata.variable.id] = 
      std::make_pair(neighbor_vdata.variable,
		     neighbor_vdata.asg);
  }

  // Construct the child and edge belief
  factor_t child_potential(child_vdata.variable);
  child_potential.uniform();
  
  factor_t edge_factor(domain_t(child_vdata.variable, parent_vdata.variable));
  edge_factor.uniform();
  
  // Loop through all factors attached to the child
  foreach(const size_t factor_id, child_vdata.factor_ids) {
    const factor_t& factor =
      shared_data->get_constant(FACTOR_ID + factor_id).as<factor_t>();
    
    // Try to build a nearly complete assignment
    const domain_t& args = factor.args();
    assignment_t asg;
    bool is_edge_potential = false;
    for(size_t i = 0; i < args.num_vars(); ++i) {
      const variable_t& var = args.var(i);
      if(var == parent_vdata.variable)  is_edge_potential = true;
      // if the ith variable is not a the parent or child we condition
      // on it
      if(var != parent_vdata.variable &&
         var != child_vdata.variable) {
        map_iterator iter = large_assignment.find(var.id);
        assert(iter != large_assignment.end());
	assignment_t local_asg(iter->second.first, 
			       iter->second.second);

        asg &= local_asg;
      }
    }
    // If the factor can be conditioned away (node_potential)
    if(is_edge_potential) {
      assert(asg.num_vars() + 2 == factor.num_vars());
      edge_factor.times_condition(factor, asg);
    } else {     
      assert(asg.num_vars() + 1 == factor.num_vars());
      child_potential.times_condition(factor, asg);
    }
  } 
  edge_factor.normalize();
  child_potential.normalize();


  // Construct the message
  
  factor_t bp_message(parent_vdata.variable);
  bp_message.convolve(edge_factor, child_potential);
  bp_message.normalize();
  factor_t conditional(parent_vdata.variable);
  assignment_t child_asg(child_vdata.variable, child_vdata.asg);
  conditional.condition(edge_factor, child_asg);
  conditional.normalize();

  // // Compute the l1 difference
  // double residual = message.residual(conditional);
  double residual = 0;
  for(size_t i =0; i < bp_message.size(); ++i) {
    residual += std::abs(bp_message.logP(i) - conditional.logP(i));    
  }
  residual /= bp_message.size();
  residual = std::tanh(residual);


  //std::cout << residual << std::endl;
  return residual;
}





// vertex_id_t get_next_root(size_t current_root, 
//                           size_t nroots,
//                           size_t nvertices) {
//   const vertex_id_t block = current_root % nroots;
//   const size_t rndnum =
//     size_t(std::floor(graphlab::random::rand01() * (nvertices/nroots)));
//   const vertex_id_t next_root =  block + rndnum*nroots;
//   // Next root is a valid vertex
//   assert(next_root < nvertices);
//   // Next root is in the same block
//   assert(next_root % nroots == block);
//   return next_root;
// }


vertex_id_t get_next_root(size_t current_root, 
                          size_t nroots,
                          size_t nvertices) {
  const vertex_id_t block = current_root % nroots;
  vertex_id_t next_root = current_root + nroots;
  if(next_root >= nvertices) next_root = block;
  // Next root is a valid vertex
  assert(next_root < nvertices);
  // Next root is in the same block
  assert(next_root % nroots == block);
  //  std::cout << current_root << "->" << next_root << std::endl;
  return next_root;
}











void grow_root_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data) {
  
  vertex_data& vdata = scope.vertex_data();
  graphlab::edge_list out_edges = scope.out_edge_ids();

  // Test to see if this is a valid root
  bool root_candidate = vdata.state == vertex_data::AVAILABLE;

  // This vertex is still a root candidate as long as all of its
  // neighbors are either available or boundary.
  if(root_candidate) {
    foreach(edge_id_t eid, out_edges) {
      const vertex_data& neighbor_vdata =
        scope.neighbor_vertex_data(scope.target(eid));
      if(!(neighbor_vdata.state == vertex_data::AVAILABLE ||
           neighbor_vdata.state == vertex_data::BOUNDARY)   ) {
        root_candidate = false; break;
      }
    }
  } // end of if

  // If it is still a root candidate then we have succeeded and can
  // grow a tree here
  if(root_candidate) {
    assert(vdata.state  == vertex_data::AVAILABLE);
    assert(vdata.parent == NULL_VID);
    assert(vdata.height == 0);
    // Make this the root
    vdata.parent = scope.vertex();
    vdata.state  = vertex_data::TREE_NODE;
    // Set the height to max tree height
    const size_t& max_height =
      shared_data->get_constant(MAX_HEIGHT_ID).as<size_t>();
    vdata.height = max_height;

    if(vdata.height == 0) {
      // Actually start the upward pass
      up_tree_update(scope, scheduler, shared_data);
    } else {
      size_t available_neighbors = 0;
      // Add all the children
      foreach(edge_id_t eid, out_edges) {
        vertex_id_t neighbor = scope.target(eid);
        vertex_data& neighbor_vdata =
          scope.neighbor_vertex_data(neighbor);
        if(neighbor_vdata.state == vertex_data::AVAILABLE) {
          edge_data& edata = scope.edge_data(eid);
          available_neighbors++;
          gl::update_task task(neighbor, grow_tree_update);
          double residual = grow_root_residual;
          
          //vdata.child_candidates++;
          vdata.child_candidates.inc();
          edata.exploring = true;
          scheduler.add_task(task, residual);

          
        } else {
          // The neighbor must be a boundary neighbor
          assert(neighbor_vdata.state == vertex_data::BOUNDARY);
        }        
      } // end of loop over out edges
      // If no neighbors were available schedule the upward pass
      if(available_neighbors == 0)
        up_tree_update(scope,scheduler,shared_data);
    }
  } else {
    // Try to start a root elsewhere
    const size_t& nroots = shared_data->get_constant(NROOTS_ID).as<size_t>();
    const vertex_id_t next_root = get_next_root(scope.vertex(),
                                          nroots,
                                          scope.num_vertices());
    const gl::update_task task(next_root, grow_root_update);
    double residual = grow_root_residual;
    scheduler.add_task(task, residual);
  }
} // End of grow root update




















void grow_tree_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data) { 
  assert(shared_data != NULL);

  // Get the vertex data
  const vertex_id_t vid = scope.vertex();
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
      compute_edge_weight(parent, vid, shared_data);
    // Get the edge to the parent
    edge_data& edata = scope.edge_data(parent_edge);
    edata.weight = edge_weight;
    // Parent remains parent depending on edge weight value
    //   double unif = gl::random::rand01();
    is_tree_candidate = pruning < edge_weight;
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
      size_t available_neighbors = 0;
      // Add all the children
      foreach(edge_id_t eid, out_edges) {
        vertex_id_t neighbor = scope.target(eid);
        edge_data& edata = scope.edge_data(eid);
        const vertex_data& neighbor_vdata = scope.neighbor_vertex_data(neighbor);        
        // Skip unavailable neighbors
        if(neighbor_vdata.state == vertex_data::AVAILABLE) {
          available_neighbors++;
          gl::update_task task(neighbor, grow_tree_update);
          // Compute the residual
          double residual = grow_tree_residual;
          if(use_weights) {
            graphlab::general_scope<graph_type>
              neighbor_scope(global_graph, neighbor, NULL);            
            // compute the edge weight
            double edge_weight =
              compute_edge_weight(vid, neighbor, shared_data);
         
            residual += edge_weight;
            edge_data& edata = scope.edge_data(eid);
            edata.weight = residual;
          }
          edata.exploring = true;
          vdata.child_candidates.inc();
          scheduler.add_task(task, residual);
        }
      }
      // If there were no available neighbors then begin upward pass
      if(available_neighbors == 0) up_tree_update(scope, scheduler, shared_data);
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
  const vertex_id_t vid = scope.vertex();
  vertex_data& vdata = scope.vertex_data();

  // If this vertex is not visited, is not a tree node, or has already
  // been calibrated simply return
  if( vdata.state != vertex_data::TREE_NODE ) return;
  
  // This vertex must have a parent
  assert(vdata.parent != NULL_VID);

  // Get neighbors
  graphlab::edge_list in_edges = scope.in_edge_ids();

  // See if this vertex is ready to be callibrated.  A vertex is ready
  // to be callibrated if all its neighbors are visited and all its
  // children are callibrated. or its height = 0
  if(vdata.height > 0) {
    foreach(edge_id_t eid, in_edges) {
      const vertex_id_t neighbor = scope.source(eid);
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



  
  // Construct an assignment map to all the neighbors that are not in
  // the tree
  typedef asg_map_type::const_iterator map_iterator;
  asg_map_type large_assignment;

  // Prepare to compute the tmp belief 
  vdata.tmp_bp_belief.uniform();

  // Collect assignments to all non tree neighbors and multiply child
  // messages into temporary bp belief
  foreach(edge_id_t eid, in_edges) {
    const vertex_id_t neighbor        = scope.source(eid);
    const vertex_data& neighbor_vdata = scope.neighbor_vertex_data(neighbor);
    const bool is_child               = neighbor_vdata.parent == vid;
    const bool is_parent              = vdata.parent == neighbor;
    // if the neighbor is neither a child or a parent then add it to
    // the conditioning set
    if(!(is_child || is_parent)) {
      large_assignment[neighbor_vdata.variable.id] = 
      std::make_pair(neighbor_vdata.variable,
		     neighbor_vdata.asg);
    } else  if(is_child) {
      // Multiply in the message from the children
      const edge_data& neighbor_edata = scope.edge_data(eid);
      vdata.tmp_bp_belief *= neighbor_edata.message;
    }
  }

  // Go through all factors that contain this variable and find those
  // factors that are not associated with tree edges
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor =
      shared_data->get_constant(FACTOR_ID + factor_id).as<factor_t>();
    // Try to build a nearly complete assignment
    const domain_t& args = factor.args();
    assignment_t asg;
    bool node_potential = true;
    for(size_t i = 0; i < args.num_vars(); ++i) {
      const variable_t& var = args.var(i);
      // if the ith variable is not the same as this vertex variable
      if(vdata.variable != var) {
        map_iterator iter = large_assignment.find(var.id);
        // If the map has a value for this variable
        if(iter == large_assignment.end()) {
          node_potential = false; break;
        } else {
	  assignment_t local_asg(iter->second.first, 
				 iter->second.second);
	  asg &= local_asg;
	}
      }
    }
    // If the factor can be conditioned away (node_potential)
    if(node_potential) {
      assert(asg.num_vars() + 1 == factor.num_vars());
      vdata.tmp_bp_belief.times_condition(factor, asg);
    }
  } 
  vdata.tmp_bp_belief.normalize();
  
  // This vertex is ready for moving up Mark self as calibrated
  vdata.state = vertex_data::CALIBRATED;  

    
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
  } else { // Upward pass
    // We need to construct the edge factor to the parent and then the
    // message
    // Construct the edge factor to parent ------------------------------
    vertex_data& parent_vdata =
      scope.neighbor_vertex_data(vdata.parent);
    const edge_id_t parent_eid     = scope.edge(vid, vdata.parent);
    edge_data& parent_edata  = scope.edge_data(parent_eid);
    const edge_id_t rev_parent_eid     = scope.edge(vdata.parent, vid);
    edge_data& rev_parent_edata  = scope.edge_data(rev_parent_eid);

    // Clear out the parent edge factor
    parent_edata.edge_factor.uniform();
    
    // Go through all factors that contain this variable and find
    // those that also contain the parent
    foreach(size_t factor_id, vdata.factor_ids) {
      const factor_t& factor =
        shared_data->get_constant(FACTOR_ID + factor_id).as<factor_t>();
      // Try to build a nearly complete assignment
      const domain_t& args = factor.args();
      assignment_t asg;
      bool is_edge_factor = false;
      for(size_t i = 0; i < args.num_vars(); ++i) {
        const variable_t& var = args.var(i);
        if(var == parent_vdata.variable) {
          assert(is_edge_factor == false);
          is_edge_factor = true;
        } else if(vdata.variable != var) {
          map_iterator iter = large_assignment.find(var.id);
          // If the map does not have a value for this variable
          if(iter == large_assignment.end()) {
            is_edge_factor = false; break;
          } else {
	    assignment_t local_asg(iter->second.first, 
				   iter->second.second);
	    asg &= local_asg;
	  } 
        }
      }
      // If the only two non-conditional variables are this and its
      // parent
      if(is_edge_factor) {
        assert(asg.num_vars() + 2 == factor.num_vars());
        parent_edata.edge_factor.times_condition(factor, asg);
      }
    } // End of loop over factors
    parent_edata.edge_factor.normalize();

    // Construct message to parent --------------------------------------
    parent_edata.message.convolve(parent_edata.edge_factor, vdata.tmp_bp_belief);


    if(rev_parent_edata.exploring) {
      assert(parent_vdata.child_candidates.value > 0);
      parent_vdata.child_candidates.dec();
      rev_parent_edata.exploring = false;
    }

    if(parent_vdata.child_candidates.value == 0) {
      // Schedule the parent
      gl::update_task task(vdata.parent, up_tree_update);
      double residual = up_tree_residual;
      scheduler.add_task(task, residual);
    }
  }
}


























void down_tree_update(gl::iscope& scope, 
                      gl::icallback& scheduler,
                      gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  // Get the vertex id
  const vertex_id_t vid = scope.vertex();
  vertex_data& vdata = scope.vertex_data();   

  if( vdata.state == vertex_data::AVAILABLE) return;
  
  // If this is a boundary vertex and it has no tree neighbors then we
  // can clear it
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
  
  // If this is not the root then the parent must have been removed
  // form the tree to proceed
  if(vdata.parent != vid) {
    const vertex_data& parent_vdata =
      scope.neighbor_vertex_data(vdata.parent);
    if(parent_vdata.parent != NULL_VID) return;
    //    if(parent_vdata.state != vertex_data::AVAILABLE) return;
  }

  // Get the in and out edges for the vertex
  graphlab::edge_list in_edges = scope.in_edge_ids();
  graphlab::edge_list out_edges = scope.out_edge_ids();    

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


  // Draw the sample conditioned on the parents value and then receive
  // message from parent if there is a parent
  factor_t sample_dist = vdata.tmp_bp_belief;
  if(vdata.parent != vid) {
    const vertex_data& parent_vdata = scope.neighbor_vertex_data(vdata.parent);
    const edge_id_t to_parent_eid = scope.edge(vid, vdata.parent);
    const edge_data& to_parent_edata = scope.edge_data(to_parent_eid);
    assignment_t parent_asg(parent_vdata.variable, parent_vdata.asg);
    sample_dist.times_condition(to_parent_edata.edge_factor, parent_asg);
    sample_dist.normalize();
    // Update the belief with message from parent
    const edge_id_t from_parent_eid = scope.edge(vdata.parent, vid);
    const edge_data& from_parent_edata = scope.edge_data(from_parent_eid);
    vdata.tmp_bp_belief *= from_parent_edata.message;
  }

  // Actually Draw the sample
  vdata.asg = sample_dist.sample().asg_at(0);
  vdata.updates++;
  vdata.belief += sample_dist;
  

  // Send messages to all children by constructing the cavity
  factor_t cavity;
  foreach(edge_id_t out_eid, out_edges) {
    const vertex_id_t  neighborid  = scope.target(out_eid);
    const vertex_data& ndata = scope.neighbor_vertex_data(neighborid);
    // If this edge is to a child in the tree
    if(vid == ndata.parent) {
      edge_data& out_edata = scope.edge_data(out_eid);
      const edge_id_t in_eid = scope.reverse_edge(out_eid);
      const edge_data& in_edata = scope.edge_data(in_eid);

      // Construct the cavity
      cavity = vdata.tmp_bp_belief;
      cavity /= in_edata.message;
      out_edata.message.convolve(in_edata.edge_factor, cavity);
      out_edata.message.normalize();
    }
  }
   
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
void single_sample_update1(gl::iscope& scope, 
			   gl::icallback& scheduler,
			   gl::ishared_data* shared_data);

template<bool UseCallback>
void single_sample_update2(gl::iscope& scope, 
			   gl::icallback& scheduler,
			   gl::ishared_data* shared_data);



template<bool UseCallback>
void single_sample_update(gl::iscope& scope, 
			  gl::icallback& scheduler,
			  gl::ishared_data* shared_data) {
  //single_sample_update1<UseCallback>(scope,scheduler, shared_data);
  single_sample_update2<UseCallback>(scope,scheduler, shared_data);

}




template<bool UseCallback>
void single_sample_update1(gl::iscope& scope, 
			   gl::icallback& scheduler,
			   gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  vertex_data& vdata = scope.vertex_data();

    
  // Construct an assignment map to all the neighbors that are not in
  // the tree
  typedef asg_map_type::const_iterator map_iterator;
  asg_map_type large_assignment;

  // Collect assignments to all non tree neighbors and multiply child
  // messages into temporary bp belief
  foreach(edge_id_t eid, scope.in_edge_ids()) {
    const vertex_id_t neighbor        = scope.source(eid);
    const vertex_data& neighbor_vdata = 
      scope.neighbor_vertex_data(neighbor);
    large_assignment[neighbor_vdata.variable.id] = 
      std::make_pair(neighbor_vdata.variable,
		     neighbor_vdata.asg);
    
  }

  // Prepare to compute the tmp belief 
  vdata.tmp_bp_belief.uniform();
  
  // Multiply in all the facotrs
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor =
      shared_data->get_constant(FACTOR_ID + factor_id).as<factor_t>();
    // Try to build a nearly complete assignment
    const domain_t& args = factor.args();
    assignment_t asg;
    for(size_t i = 0; i < args.num_vars(); ++i) {
      const variable_t& var = args.var(i);
      if(var != vdata.variable) {
        map_iterator iter = large_assignment.find(var.id);
        assert(iter != large_assignment.end());
	assignment_t local_asg(iter->second.first, 
			       iter->second.second);
	asg &= local_asg;
      }
    }
    assert(asg.num_vars() + 1 == factor.num_vars());
    vdata.tmp_bp_belief.times_condition(factor, asg);    
  } 
  vdata.tmp_bp_belief.normalize(); 
  vdata.asg = vdata.tmp_bp_belief.sample().asg_at(0);
  vdata.belief += vdata.tmp_bp_belief;
  vdata.updates++;
  
  // Reschedule self
  if(UseCallback) {
    gl::update_task task(scope.vertex(), single_sample_update<UseCallback> );
    double residual = 1.0;
    scheduler.add_task(task, residual);
  }
}


template<bool UseCallback>
void single_sample_update2(gl::iscope& scope, 
                          gl::icallback& scheduler,
                          gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  vertex_data& vdata = scope.vertex_data();

  factor_t conditional;
  factor_t& belief(vdata.tmp_bp_belief);
  belief.uniform();
  conditional.set_args(vdata.variable);
  foreach(vertex_id_t factor_id, vdata.factor_ids) {
    const factor_t& factor =
      shared_data->get_constant(FACTOR_ID + factor_id).as<factor_t>();
    // build the conditional
    assignment_t conditional_asg = factor.args() - vdata.variable;
    for(size_t i = 0; i < conditional_asg.num_vars(); ++i) {
      const vertex_data& other_vdata = 
	scope.const_neighbor_vertex_data(conditional_asg.args().var(i).id);
      assert(conditional_asg.args().var(i) == other_vdata.variable);
      conditional_asg &= 
	assignment_t(other_vdata.variable, other_vdata.asg);
    }
    belief.times_condition(factor, conditional_asg);
  }
  vdata.tmp_bp_belief.normalize(); 
  vdata.asg = vdata.tmp_bp_belief.sample().asg_at(0);
  vdata.belief += vdata.tmp_bp_belief;
  vdata.updates++;
  
  // Reschedule self
  if(UseCallback) {
    gl::update_task task(scope.vertex(), single_sample_update<UseCallback> );
    double residual = 1.0;
    scheduler.add_task(task, residual);
  }
}











// gl::iengine* make_colored_engine(gl::graph& graph,
//                                  gl::thread_shared_data& sdm,
//                                  size_t ncpus) {

//   std::cout << "Computing coloring " << std::endl;
//   graph.compute_coloring();

//   gl::iengine* engine =
//     graphlab::engine_factory::new_engine("threaded",
//                                          "colored",
//                                          "null",
//                                          graph,
//                                          ncpus);

//   assert(engine != NULL);  
//   // Set the shared data 
//   engine->set_shared_data_manager(&sdm);
//   std::cout << "Using set schedule planner." << std::endl;
  
//   const bool use_callback = false;
//   gl::update_function update_function =
//     single_sample_update<use_callback>;

//   engine->get_scheduler().set_option(gl::scheduler_options::UPDATE_FUNCTION, 
//                                      (void*) update_function);
//   return engine;
// } // end of make colored engine









gl::iengine* make_engine(gl::graph& graph,
                         gl::thread_shared_data& sdm,
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

  return engine;
}



void add_factors_to_sdm(gl::ishared_data_manager& sdm,
                        const std::vector<factor_t>& factors) {
  // Add all the edge factors to the shared data
  for(size_t i = 0; i < factors.size(); ++i) {
    sdm.set_constant(FACTOR_ID + i, graphlab::any(factors[i]) );
  }
}





void set_tree_sampler_constants(gl::ishared_data_manager& sdm,
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















#include <graphlab/macros_undef.hpp>
#endif
