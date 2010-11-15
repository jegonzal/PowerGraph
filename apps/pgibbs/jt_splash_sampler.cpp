
#include "jt_splash_sampler.hpp"
#include "pgibbs_tls.hpp"






#include <graphlab/macros_def.hpp>











void jtree_sample_update(jtree_gl::iscope& scope,
                         jtree_gl::icallback& callback,
                         jtree_gl::ishared_data* shared_data) {
    typedef factorized_model::factor_map_t factor_map_t;
    
    // get the vertex data
    jtree_vertex_data& vdata = scope.vertex_data();
    
    // get thread local storage to reduce hit on allocator
    pgibbs_tls& tls = get_pgibbs_tls();

    //////////////////////////////////////////////////////////////////
    // Initialize factor    
    
    // If the factor args have not been set then we need to initialize
    // the local factor by setting the args and taking the product of
    // all factors associated with the clique.  Some of these factors
    // may depend on variables not in the clique and are therefore
    // sliced (conditioned) on the current assignment to those
    // variables.
    if(vdata.factor.args() != vdata.variables) {
      assert(shared_data != NULL);
      // Resize the factor for the variables in the clique
      vdata.factor.set_args(vdata.variables);
      vdata.factor.uniform();

      // We now build up the factor by iteratoring over the dependent
      // factors conditioning if necessary into the conditional_factor
      // and then multiplying.
      factor_t& conditional_factor(tls.conditional_factor);
      // Iterate over the factors and multiply each into this factor
      foreach(size_t factor_id, vdata.factor_ids) {
        const factor_t& factor = get_factor(*shared_data, factor_id);
        // Build up an assignment for the conditional
        domain_t conditional_args = factor.args() - vdata.variables;
        assignment_t conditional_asg;
        for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
          const mrf_vertex_data& mrf_vdata = 
            get_mrf_vdata(*shared_data, conditional_args.var(i).id);
          assert(mrf_vdata.tree_info.tree_id == NULL_VID);
          conditional_asg &= 
	    assignment_t(mrf_vdata.variable, mrf_vdata.asg);         
        }
        // set the factor arguments
        conditional_factor.set_args(factor.args() - conditional_args);
        conditional_factor.condition(factor, conditional_asg);        
        // Multiply the conditional factor in
        vdata.factor *= conditional_factor;
      }
      // Extra normalization for stability on the table factors
      vdata.factor.normalize();
      // vdata.belief = vdata.factor;
    }

    //////////////////////////////////////////////////////////////////
    // receive any unreceived messages
    size_t received_neighbors = 0;
    if(!vdata.calibrated) {
      foreach(edge_id_t in_eid, scope.in_edge_ids()) {
        jtree_edge_data& in_edata = scope.edge_data(in_eid);
        // if the message has been calibrated but not received
        if(in_edata.calibrated && !in_edata.received) {
          // receive message and mark as calibrated
          vdata.factor *= in_edata.message;
          vdata.factor.normalize();
          in_edata.received = true;
        }
        // track total received neighbors
        if(in_edata.received) received_neighbors++;
      } // end of receive all in messages
      // if all messages have been received then set as calibrated
      vdata.calibrated = 
        received_neighbors == scope.in_edge_ids().size();
    } else {
      received_neighbors = scope.in_edge_ids().size();
    }


    //////////////////////////////////////////////////////////////////
    // send any unset messages 
    // if we have recieve enough in messages
    if(received_neighbors + 1 >= scope.in_edge_ids().size()) {
      factor_t& cavity(tls.cavity);
      foreach(edge_id_t out_eid, scope.out_edge_ids()) {
        jtree_edge_data& out_edata = scope.edge_data(out_eid);
        edge_id_t rev_eid = scope.reverse_edge(out_eid);
        // if the out message is not calibrated try to calibrate it:
        if(!out_edata.calibrated) {
          bool ready_to_send = true;
          // check that all in messages (except the one we want to
          // send) have been recieved
          foreach(edge_id_t in_eid, scope.in_edge_ids()) {
            const jtree_edge_data& in_edata = scope.const_edge_data(in_eid);
            // if the in edge has not been received and is not from
            // the destination of the out edge then we cannot send
            if(!in_edata.received && rev_eid != in_eid) {
              ready_to_send = false;
              break;
            }
          } // check all neighbors are go for send

          // if we are ready to send then compute message
          if(ready_to_send) {
            cavity = vdata.factor;
            const jtree_edge_data& in_edata = scope.const_edge_data(rev_eid);
            // construct cavity if necessary
            if(in_edata.received) {
              cavity /= in_edata.message;
              cavity.normalize();
            }
            // compute actual message
            out_edata.message.set_args(out_edata.variables);
            out_edata.message.marginalize(cavity);
            out_edata.message.normalize();
            out_edata.calibrated = true;
            // schedule the reception of the message
            callback.add_task(scope.target(out_eid), jtree_sample_update, 1.0);      
          } // end of if ready to send
        } // end of if not calibrated
      } // end of loop over outbound messages
    } // of send all out messages


    //////////////////////////////////////////////////////////////////
    // Construct RB estimate and Sample
    // if calibrated but not yet sampled
    if(vdata.calibrated && !vdata.sampled) {
      // check that the parent is sampled and also determine which
      // variables are going to be sampled at this clique.  This is
      // done by finding the parent assignment if there is one
      assignment_t parent_asg;
      edge_id_t to_parent_eid = NULL_EID;

      // find the parent
      bool parent_found = false;
      foreach(edge_id_t out_eid, scope.out_edge_ids()) {       
        const jtree_vertex_data& parent_vdata = 
          scope.const_neighbor_vertex_data(scope.target(out_eid));
        if(parent_vdata.sampled) {
          assert(parent_vdata.calibrated);
          assert(!parent_found);
          parent_found = true;
          to_parent_eid = out_eid;
          const jtree_edge_data& parent_edata = scope.const_edge_data(to_parent_eid);
          parent_asg = 
            parent_vdata.asg.restrict(parent_edata.variables);
          assert(parent_asg.args() == parent_edata.variables);            
          // break;
        }
      }

      
      // Determine the remaining variables for which we will need to
      // sample and construct RB estimates
      domain_t unsampled_variables = 
        vdata.variables - parent_asg.args();
      vdata.asg = parent_asg;

      // if there are unsampled variables then sample them
      if(unsampled_variables.num_vars() > 0) {
        // First update all the RB estimates for the unsampled
        // variables in the mrf graph
        factor_t& tmp_belief(tls.tmp_belief);
        for(size_t i = 0; i < unsampled_variables.num_vars(); ++i) {
          variable_t var = unsampled_variables.var(i);
          // Construct the RB belief estimate
          tmp_belief.set_args(var);
          tmp_belief.marginalize(vdata.factor);
          tmp_belief.normalize();
          // Update the MRF
          mrf_vertex_data& mrf_vdata = get_mrf_vdata(*shared_data, var.id);
          mrf_vdata.belief += tmp_belief;
        } 

        // Condition the belief on the parent assignmnet
        tmp_belief.set_args(unsampled_variables);
        tmp_belief.condition(vdata.factor, parent_asg);
        tmp_belief.normalize();

        // Sample the remaining variables from the belief
        assignment_t sample_asg = tmp_belief.sample();
        // Set the local assignment
        vdata.asg = sample_asg & parent_asg;
        // the assignment should exacty cover the variables
        assert(vdata.asg.args() == vdata.variables);

      
        //// Fill out the MRF with the sampled variables
        for(size_t i = 0; i < sample_asg.num_vars(); ++i) {
          variable_t var = sample_asg.args().var(i);
          mrf_vertex_data& mrf_vdata = get_mrf_vdata(*shared_data, var.id);
          assignment_t local_asg = sample_asg.restrict(var);
          if(mrf_vdata.asg != local_asg.asg_at(0)) {
            mrf_vdata.nchanges++;
          }
          mrf_vdata.asg = local_asg.asg_at(0);
          mrf_vdata.nsamples++;
          // std::cout << graphlab::thread::thread_id()
          //           << ": sampling " << mrf_vdata.variable << std::endl;
          // remove the vertex from any trees
          mrf_vdata.tree_info.tree_id = NULL_VID;
          mrf_vdata.tree_info.height = 0;
          
          // double& logP = mrf_vdata.belief.logP(mrf_vdata.asg.asg_at(0));
          // logP = std::log( std::exp(logP) + 1.0 );          
        } 
      } // end of sampling unsampled variables

      // mark as sampled
      vdata.sampled = true;

      // Reschedule unsampled neighbors
      foreach(edge_id_t out_eid, scope.out_edge_ids()) {
        if(out_eid != to_parent_eid) {
          const vertex_id_t neighbor_vid = scope.target(out_eid);
          assert(neighbor_vid < scope.num_vertices());
          callback.add_task(neighbor_vid, 
                            jtree_sample_update, 
                            1.0);
        }
      }
    } // End of if(!sampled) sampling procedure

  } // End of update function




















void jt_worker::init(size_t wid,
                     scope_factory_type& sf, 
                     const factorized_model::factor_map_t& factors,
                     const std::vector<vertex_id_t>& root_perm, 
                     size_t ncpus,         
                     size_t treesize,
                     size_t treewidth,
                     size_t factorsize,
                     size_t max_height,
                     bool priorities,
                     size_t internal_threads) {
  // Initialize parameters
  worker_id        = wid;
  scope_factory    = &sf;  
  roots            = &root_perm;    
  root_index = root_perm.size();
  current_root = worker_id;


  worker_count     = ncpus;
  max_tree_size    = treesize;
  max_tree_width   = treewidth;
  max_factor_size  = factorsize;
  max_tree_height  = max_height;
  use_priorities   = priorities;

  
  // Initialize local jtcore
  if(internal_threads > 1) {
    jt_core.set_scheduler_type("multiqueue_fifo");
    jt_core.set_scope_type("edge");
    jt_core.set_ncpus(internal_threads);
    jt_core.set_engine_type("async");
  } else {
    jt_core.set_scheduler_type("fifo");
    jt_core.set_scope_type("none");
    jt_core.set_ncpus(1);
    jt_core.set_engine_type("async_sim");
  }
  
  // Initialize the shared data manager for the jt core
  set_factor_map(factors, jt_core.shared_data());
  set_mrf_graph(sf.get_graph(), jt_core.shared_data());
} // end of init













// get a root
void jt_worker::run() {   
  // looup until runtime is reached
  while(graphlab::lowres_time_seconds() < finish_time_seconds) {
    /////////////////////////////////////////////////////////
    // Construct one tree (we must succeed in order to count a tree
    size_t sampled_variables = 0;
    //      move_to_next_root();
    while(sampled_variables == 0 && 
          graphlab::lowres_time_seconds() < finish_time_seconds) {
      move_to_next_root();
      sampled_variables = sample_once();
      if(sampled_variables == 0) {
        ncollisions++;
        //  sched_yield();
      }
    }

    //       // Get a local copy of the graph
    //       mrf::graph_type& mrf(scope_factory->get_graph());
    //       if(worker_id == 0) {
    //         std::cout << "Saving sample: " << std::endl;
    //         size_t rows = std::sqrt(mrf.num_vertices());
    //         image img(rows, rows);
    //         for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {   
    //           img.pixel(vid) = mrf.vertex_data(vid).asg.asg_at(0);
    //         }
    //         img.pixel(0) = 0;
    //         img.pixel(1) = mrf.vertex_data(0).variable.arity -1;
    //         img.save(make_filename("sample", ".pgm", tree_count).c_str());
    //       }

    ntrees++;
  } 
} // end of run










/**
 * Grab this vertex into the tree owned by worker id
 */
bool jt_worker::quick_try_vertex(vertex_id_t vid) {
  const mrf_graph_type& mrf(scope_factory->get_graph());
  const mrf_vertex_data& vdata = mrf.vertex_data(vid);
  // Check that this vertex is not already in a tree
  bool in_tree = vdata.tree_info.tree_id != NULL_VID;
  if(in_tree) return false;
  // check that the neighbors are not in any other trees than this
  // one
  const graphlab::edge_list& in_eids = mrf.in_edge_ids(vid);
  foreach(edge_id_t in_eid, in_eids) {
    vertex_id_t neighbor_vid = mrf.source(in_eid);
    const mrf_vertex_data& vdata = mrf.vertex_data(neighbor_vid);
    bool in_tree = vdata.tree_info.tree_id != NULL_VID;
    // if the neighbor is in a tree other than this one quit
    if(in_tree && worker_id != vdata.tree_info.tree_id) return false;
  }
  return true;
} // end of try grab vertex















/**
 * Grab this vertex into the tree owned by worker id
 */
bool jt_worker::try_grab_vertex(iscope_type& scope) {
  // Check that this vertex is not already in a tree
  bool in_tree = scope.vertex_data().tree_info.tree_id != NULL_VID;
  if(in_tree) return false;

  // check that the neighbors are not in any other trees than this
  // one
  foreach(edge_id_t in_eid, scope.in_edge_ids()) {
    vertex_id_t neighbor_vid = scope.source(in_eid);
    const mrf_vertex_data& vdata = 
      scope.const_neighbor_vertex_data(neighbor_vid);
    bool in_tree = vdata.tree_info.tree_id != NULL_VID;
    // if the neighbor is in a tree other than this one quit
    if(in_tree && worker_id != vdata.tree_info.tree_id) return false;
  }
  // Assert that this vertex is not in a tree and that none of the
  // neighbors are in other trees
  // This vertex does not neighbor any other trees than this one
  scope.vertex_data().tree_info.tree_id = worker_id;
  return true;
} // end of try grab vertex












/**
 * Release the vertex
 */
void jt_worker::release_vertex(iscope_type& scope) {
  // This vertex does not neighbor any other trees than this one
  scope.vertex_data().tree_info.tree_id = NULL_VID;
} // release the vertex












double jt_worker::score_vertex(vertex_id_t vid) {

  mrf_graph_type& mrf(scope_factory->get_graph());
  mrf_vertex_data& vdata = mrf.vertex_data(vid);

  //    return graphlab::random::rand01();

  //    return score_vertex_log_odds(vid); 
  // return score_vertex_lik(vid);
  if (vdata.nsamples < 100 || vdata.tree_info.priority < 0) {
    vdata.tree_info.priority = score_vertex_log_odds(vid); 
  }
  return vdata.tree_info.priority;

  // +
  //    return double(vdata.changes + graphlab::random::rand01()) 
  //      / sqrt(double(vdata.updates + 1));

  // if(vdata.updates < 100) {


  //    return score_vertex_lik(vid); // +
  // }
  //	graphlab::random::rand01();
  // return score_vertex_l1_diff(vid);
  //}
  //    return (1.0 + graphlab::random::rand01() + score) / (vdata.updates + 1); 
  //return graphlab::random::rand01();;
}


















double jt_worker::score_vertex_l1_diff(vertex_id_t vid) {
  // Get the scope factory
  const mrf_graph_type& mrf(scope_factory->get_graph());
  const mrf_vertex_data& vdata = mrf.vertex_data(vid);

  // Construct the domain of neighbors that are already in the tree
  domain_t vars = vdata.variable;
  foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
    const vertex_id_t neighbor_vid = mrf.source(ineid);
    const mrf_vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
    // test to see if the neighbor is in the tree by checking the
    // elimination time map
    if(jt_list.is_eliminated(neighbor_vid)) {
      vars += neighbor.variable;
      // If this vertex has too many tree neighbor than the priority
      // is set to -1;
      if(vars.num_vars() > max_tree_width) return -1;
      if(vars.size() > max_factor_size) return -1;
    } 
  }

  // Compute the clique factor
  clique_factor.set_args(vars);
  clique_factor.uniform();
  // get all the factors
  const factorized_model::factor_map_t& 
    factors(get_factor_map(jt_core.shared_data()));
  // Iterate over the factors and multiply each into this factor
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor = factors[factor_id];      
    // Build up an assignment for the conditional
    domain_t conditional_args = factor.args() - vars;
    if(conditional_args.num_vars() > 0) {
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& neighbor_vdata = 
          mrf.vertex_data(conditional_args.var(i).id);
        conditional_asg &= 
          assignment_t(neighbor_vdata.variable, neighbor_vdata.asg);
      }
      // set the factor arguments
      conditional_factor.set_args(factor.args() - conditional_args);
      conditional_factor.condition(factor, conditional_asg);        
      // Multiply the conditional factor in
      clique_factor *= conditional_factor;
      //       clique_factor.normalize();
    } else {
      clique_factor *= factor;
    }
  } // end of loop over factors
  clique_factor.normalize();


  // Compute the product of marginals
  product_of_marginals_factor.set_args(vars);
  product_of_marginals_factor.uniform();
  for(size_t i = 0; i < vars.num_vars(); ++i) {
    marginal_factor.set_args(vars.var(i));
    marginal_factor.marginalize(clique_factor);
    marginal_factor.normalize();
    product_of_marginals_factor *= marginal_factor;
  }
  product_of_marginals_factor.normalize();

  // Compute the residual
  double residual = clique_factor.l1_diff(product_of_marginals_factor);


  assert( residual >= 0);
  assert( !std::isnan(residual) );
  assert( std::isfinite(residual) );

  // ensure score is bounded
  //    residual = std::tanh(residual);

  return residual;

} // end of score l1 diff












double jt_worker::score_vertex_log_odds(vertex_id_t vid) {
  // Get the scope factory
  const mrf_graph_type& mrf(scope_factory->get_graph());
  const mrf_vertex_data& vdata(mrf.vertex_data(vid));

  // Construct the domain of neighbors that are already in the tree
  domain_t vars = vdata.variable;
  foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
    const vertex_id_t neighbor_vid = mrf.source(ineid);
    const mrf_vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
    // test to see if the neighbor is in the tree by checking the
    // elimination time map
    if(jt_list.is_eliminated(neighbor_vid)) {
      vars += neighbor.variable;
      // If this vertex has too many tree neighbor than the priority
      // is set to 0;
      if(vars.num_vars() > max_tree_width) return -1;
      if(vars.size() > max_factor_size) return -1;
    } 
  }
    
  assert(vars.num_vars() == 2);


  // Compute the clique factor
  clique_factor.set_args(vars);
  clique_factor.uniform();
  // get all the factors
  const factorized_model::factor_map_t& 
    factors(get_factor_map(jt_core.shared_data()));

  // Iterate over the factors and multiply each into this factor
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor = factors[factor_id];      
    // Build up an assignment for the conditional
    domain_t conditional_args = factor.args() - vars;
    if(conditional_args.num_vars() > 0) {
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& neighbor_vdata = 
          mrf.vertex_data(conditional_args.var(i).id);
        conditional_asg &= 
          assignment_t(neighbor_vdata.variable, neighbor_vdata.asg);
      }
      // set the factor arguments
      conditional_factor.set_args(factor.args() - conditional_args);
      conditional_factor.condition(factor, conditional_asg);        
      // Multiply the conditional factor in
      clique_factor *= conditional_factor;
      //        clique_factor.normalize();
    } else {
      clique_factor *= factor;
    }
  } // end of loop over factors
    // Compute the conditional factor and marginal factors
  conditional_factor.set_args(vars - vdata.variable);
  conditional_factor.condition(clique_factor, 
                               assignment_t(vdata.variable, vdata.asg));  
  marginal_factor.set_args(vars - vdata.variable);
  marginal_factor.marginalize(clique_factor);
    
  // Compute metric
  conditional_factor.normalize();
  marginal_factor.normalize();
  // double residual = conditional_factor.l1_logdiff(marginal_factor);
  double residual = conditional_factor.l1_diff(marginal_factor);

  // rescale by updates
  //    residual = residual / (vdata.updates + 1);

  assert( residual >= 0);
  assert( !std::isnan(residual) );
  assert( std::isfinite(residual) );

  // ensure score is bounded
  //    residual = std::tanh(residual);

  return residual;
} // end of score vertex























double jt_worker::score_vertex_lik(vertex_id_t vid) {
  // Get the scope factory
  const mrf_graph_type& mrf(scope_factory->get_graph());
  const mrf_vertex_data& vdata(mrf.vertex_data(vid));

  // Construct the domain of neighbors that are already in the tree
  domain_t vars = vdata.variable;
  foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
    const vertex_id_t neighbor_vid = mrf.source(ineid);
    const mrf_vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
    // test to see if the neighbor is in the tree by checking the
    // elimination time map
    if(jt_list.is_eliminated(neighbor_vid)) {
      vars += neighbor.variable;
      // If this vertex has too many tree neighbor than the priority
      // is set to 0;
      if(vars.num_vars() > max_tree_width) return -1;
      if(vars.size() > max_factor_size) return -1;
    } 
  }
    
  // Compute the clique factor
  clique_factor.set_args(vars);
  clique_factor.uniform();
  // get all the factors
  const factorized_model::factor_map_t& 
    factors(get_factor_map(jt_core.shared_data()));

  // Iterate over the factors and multiply each into this factor
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor = factors[factor_id];      
    // Build up an assignment for the conditional
    domain_t conditional_args = factor.args() - vars;
    if(conditional_args.num_vars() > 0) {
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& neighbor_vdata = 
          mrf.vertex_data(conditional_args.var(i).id);
        conditional_asg &= 
          assignment_t(neighbor_vdata.variable, neighbor_vdata.asg);
      }
      // set the factor arguments
      conditional_factor.set_args(factor.args() - conditional_args);
      conditional_factor.condition(factor, conditional_asg);        
      // Multiply the conditional factor in
      clique_factor *= conditional_factor;
      //        clique_factor.normalize();
    } else {
      clique_factor *= factor;
    }
  } // end of loop over factors

    // Compute the conditional factor and marginal factors
  marginal_factor.set_args(vdata.variable);
  marginal_factor.marginalize(clique_factor);
  marginal_factor.normalize();
  double residual =  1.0 - exp(marginal_factor.logP(vdata.asg));

  assert( residual >= 0);
  assert( !std::isnan(residual) );
  assert( std::isfinite(residual) );

  // // ensure score is bounded
  // residual = std::tanh(residual);


  return residual;
} // end of max lik


















void jt_worker::grow_bfs_jtree() {
  assert(scope_factory != NULL);
  // Get the scope factory
  mrf_graph_type& mrf = scope_factory->get_graph();
  // Clear local data structures
  jt_list.clear();
  bfs_queue.clear();
  visited.clear();
     
  // add the root
  bfs_queue.push_back(current_root);
  visited.insert(current_root);

  while(!bfs_queue.empty()) {
    // Take the top element
    const vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop_front();

    // pretest that the vertex is available before trying to get it
    bool grabbed = quick_try_vertex(next_vertex);
    if(!grabbed) continue;

    // Maybe we can get the vertex so actually try to get it
    iscope_type* scope_ptr = 
      scope_factory->get_edge_scope(worker_id, next_vertex);
    assert(scope_ptr != NULL);
    iscope_type& scope(*scope_ptr);

    // See if we can get the vertex for this tree
    grabbed = try_grab_vertex(scope);

    // If we failed to grab the scope then skip this vertex
    if(grabbed) {

      // if this is the root vertex
      vertex_id_t min_height = 0;
      bool is_root = jt_list.cliques.empty();
      if(max_tree_height != 0 && !is_root) {
        min_height = max_tree_height;
        // find the closest vertex to the root
        foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
          vertex_id_t neighbor_vid = mrf.target(eid);
          // if the neighbor is already in the tree
          if(jt_list.is_eliminated(neighbor_vid)) {
            min_height = 
              std::min(min_height, 
                       mrf.vertex_data(neighbor_vid).tree_info.height + 1);
          } 
        }
      } 

      // Check if it is safe to extend to the tree to include next variable
      bool safe_extension =
        (is_root ||  
         (max_tree_height == 0) ||
         (min_height < max_tree_height))
        &&
        extend_jtree_list(next_vertex,
                          mrf,
                          max_tree_width,
                          max_factor_size,
                          jt_list);

      if(safe_extension) {   
        // set the height
        mrf.vertex_data(next_vertex).tree_info.height = min_height;

        // add the neighbors to the search queue
        foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
          vertex_id_t neighbor_vid = mrf.target(eid);
          if(visited.count(neighbor_vid) == 0) {
            bfs_queue.push_back(neighbor_vid);
            visited.insert(neighbor_vid);
          }
        }
      } else {
        // release the vertex since it could not be used in the tree
        release_vertex(scope);
      }
    } // end of grabbed
      // release the scope
    scope_factory->release_scope(&scope);        
    // Limit the number of variables
    if(jt_list.cliques.size() > max_tree_size) break;
  } // end of while loop
} // end grow_bfs_jtree



void jt_worker::grow_prioritized_jtree() {
  assert(scope_factory != NULL);
  // Get the scope factory
  mrf_graph_type& mrf = scope_factory->get_graph();
  // Clear local data structures
  jt_list.clear();
  priority_queue.clear();
  visited.clear();
     
  // add the root
  priority_queue.push(current_root, 1.0);
  visited.insert(current_root);

  while(!priority_queue.empty()) {
    // Take the top element
    const vertex_id_t next_vertex = priority_queue.pop().first;

    // pretest that the vertex is available before trying to get it
    bool grabbed = quick_try_vertex(next_vertex);
    if(!grabbed) continue;


    // Get the scope
    bool released_scope = false;
    iscope_type* scope_ptr = 
      scope_factory->get_edge_scope(worker_id, next_vertex);
    assert(scope_ptr != NULL);
    iscope_type& scope(*scope_ptr);

    // See if we can get the vertex for this tree
    grabbed = try_grab_vertex(scope);

    // If we failed to grab the scope then skip this vertex
    if(grabbed) {
      // compute the tree height of the new vertex
      vertex_id_t min_height = 0;        
      // if this is not the root and we care about tree height
      bool is_root = jt_list.cliques.empty();
      if(max_tree_height != 0 && !is_root) {
        min_height = max_tree_height;
        // find the closest vertex to the root
        foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
          vertex_id_t neighbor_vid = mrf.target(eid);
          // if the neighbor is already in the tree
          if(jt_list.is_eliminated(neighbor_vid)) {
            min_height = 
              std::min(min_height, 
                       mrf.vertex_data(neighbor_vid).tree_info.height + 1);
          } 
        }
      } // end of tree height check for non root vertex 
          
        // test the 
      bool safe_extension =
        (is_root ||
         (max_tree_height == 0) ||
         (min_height < max_tree_height))
        &&
        extend_jtree_list(next_vertex,
                          mrf,
                          max_tree_width,
                          max_factor_size,
                          jt_list);


      // If the extension was safe than the elim_time_map and
      // cliques data structure are automatically extended
      if(safe_extension) {
        // // set the height
        mrf.vertex_data(next_vertex).tree_info.height = min_height;

        // release the scope early to allow other processors to move
        // in
        released_scope = true;
        scope_factory->release_scope(&scope);        

        // add the neighbors to the search queue or update their priority
        foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
          vertex_id_t neighbor_vid = mrf.target(eid);          
          if(visited.count(neighbor_vid) == 0) {
            // Vertex has not yet been visited
            double score = score_vertex(neighbor_vid);
            // if the score is greater than zero then add the
            // neighbor to the priority queue.  The score is zero if
            // there is no advantage or the treewidth is already too
            // large
            if(score >= 0) priority_queue.push(neighbor_vid, score);
            visited.insert(neighbor_vid);

          } 
          // else if(priority_queue.contains(neighbor_vid)) {
          //   // vertex is still in queue we may need to recompute
          //   // score
          //   double score = score_vertex(neighbor_vid);
          //   if(score >= 0) {
          //     // update the priority queue with the new score
          //     priority_queue.update(neighbor_vid, score);
          //   } else {
          //     // The score computation revealed that the clique
          //     // would be too large so simply remove the vertex from
          //     // the priority queue
          //     priority_queue.remove(neighbor_vid);
          //   }
          // } // otherwise the vertex has been visited and processed

        }
      } else {
        // release the vertex since it could not be used in the tree
        release_vertex(scope);
      }
    }

    if(!released_scope) {
      // release the scope
      scope_factory->release_scope(&scope);        
    }
    // Limit the number of variables
    if(jt_list.cliques.size() > max_tree_size) break;
  } // end of while loop
    
} // end grow_prioritized_jtree























size_t jt_worker::sample_once() {
  if(use_priorities) {
    // grow the prioritized junction tree data structure
    grow_prioritized_jtree();
  } else {
    // grow the bfs junction tree data structure
    grow_bfs_jtree();
  }

  assert(scope_factory != NULL);
  // Get the scope factory
  mrf_graph_type& mrf = scope_factory->get_graph();
    
  // If we failed to build a tree return failure
  if(jt_list.cliques.empty()) return 0;

  //  std::cout << "Varcount: " << cliques.size() << std::endl;  

  // ///////////////////////////////////
  // // plot the graph
  // if(worker_id == 0) {
  //   std::cout << "Saving treeImage:" << std::endl;
  //   size_t rows = std::sqrt(mrf.num_vertices());
  //   image img(rows, rows);
  //   for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {
  //     vertex_id_t tree_id = mrf.vertex_data(vid).tree_id;
  //     img.pixel(vid) = 
  //         tree_id == vertex_id_t(-1)? 0 : tree_id + worker_count;
  //   }
  //   img.save(make_filename("tree", ".pgm", tree_count).c_str());
  // }



  // Build the junction tree and sample
  jt_core.graph().clear();
  size_t num_factors = get_num_factors(jt_core.shared_data());
  jtree_list_to_jtree_graph(jt_list, mrf, num_factors, jt_core.graph());


  // jtree_from_cliques(mrf,  
  //                    cliques.begin(), cliques.end(), 
  //                    jt_core.graph());

  // Rebuild the engine (clear the old scheduler)
  jt_core.rebuild_engine();
  // add tasks to all vertices
  jt_core.add_task_to_all(jtree_sample_update, 1.0);

  // Run the core
  jt_core.start();



  // Check that the junction tree is sampled

  size_t actual_tree_width = 0;
  for(vertex_id_t vid = 0; 
      vid < jt_core.graph().num_vertices(); ++vid) {
    const jtree_vertex_data& vdata = jt_core.graph().vertex_data(vid);
    assert(vdata.sampled);
    assert(vdata.calibrated);
    assert(vdata.variables.num_vars() > 0);
    actual_tree_width = 
      std::max(vdata.variables.num_vars() - 1, actual_tree_width); 
  } 
    
  // std::cout << "Treewidth: " << actual_tree_width << std::endl;

  // Return the number of variables in the tree
  return jt_list.elim_time.size();
} // end of sample once


















jt_splash_sampler::
jt_splash_sampler(mrf_gl::core& mrf_core,
                  const graphlab::engine_options& eopts,
                  size_t max_tree_size,
                  size_t max_tree_width,
                  size_t max_factor_size,
                  size_t max_tree_height,
                  size_t internal_threads,
                  bool use_priorities) :
  workers(eopts.ncpus),
  scope_factory(mrf_core.graph(), eopts.ncpus, 
                graphlab::scope_range::EDGE_CONSISTENCY),
  roots(mrf_core.graph().num_vertices()),
  use_cpu_affinity(eopts.enable_cpu_affinities) { 

  // Shuffle ther oot ordering 
  for(vertex_id_t vid = 0; vid < mrf_core.graph().num_vertices(); ++vid)
    roots[vid] = vid;
  std::random_shuffle(roots.begin(), roots.end());
       
  for(size_t i = 0; i < eopts.ncpus; ++i) {
    // Initialize the worker
    workers[i].init(i, 
                    scope_factory, 
                    get_factor_map(mrf_core.shared_data()),
                    roots,
                    eopts.ncpus,    
                    max_tree_size,
                    max_tree_width,
                    max_factor_size,
                    max_tree_height,
                    use_priorities,
                    internal_threads);    
  }
} // end of constructor





size_t jt_splash_sampler::total_collisions() const {
  size_t total_collisions = 0;
  foreach(const jt_worker& worker, workers) 
    total_collisions += worker.ncollisions;
  return total_collisions;
}

size_t jt_splash_sampler::total_trees() const {
  size_t total_trees = 0;
  foreach(const jt_worker& worker, workers) 
    total_trees += worker.ntrees;
  return total_trees;
}










  
void jt_splash_sampler::sample_once(float runtime_secs) {
  // create workers
  graphlab::thread_group threads;
  
  for(size_t i = 0; i < workers.size(); ++i) {
    workers[i].set_runtime(runtime_secs);
    // Launch the threads
    if(use_cpu_affinity) threads.launch(&(workers[i]), i);
    else threads.launch(&(workers[i]));            
  }
  
  // Wait for all threads to finish
  threads.join();
  
}                   










