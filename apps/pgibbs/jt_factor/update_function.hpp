#ifndef PGIBBS_JT_UPDATE_FUNCTION_HPP
#define PGIBBS_JT_UPDATE_FUNCTION_HPP

#include <graphlab.hpp>


#include "data_structures.hpp"




#include <graphlab/macros_def.hpp>



namespace junction_tree{
  enum SDT_KEYS {FACTOR_KEY, MRF_KEY};


  /** Slow update */
  void slow_update(gl::iscope& scope,
                   gl::icallback& callback,
                   gl::ishared_data* shared_data);
 
  void fast_update(gl::iscope& scope,
                   gl::icallback& callback,
                   gl::ishared_data* shared_data);


  void calibrate_update(gl::iscope& scope,
                        gl::icallback& callback,
                        gl::ishared_data* shared_data) {
    
    // slow_update(scope, callback, shared_data);

    fast_update(scope, callback, shared_data);

  } // End of update function




















  void slow_update(gl::iscope& scope,
                   gl::icallback& callback,
                   gl::ishared_data* shared_data) {
    typedef factorized_model::factor_map_t factor_map_t;
    
    // get the vertex data
    vertex_data& vdata = scope.vertex_data();
    
    
    // If the factor args have not been set then we need to initialize
    // the local factor by setting the args and taking the product of
    // all factors associated with the clique.  Some of these factors
    // may depend on variables not in the clique and are therefore
    // sliced (conditioned) on the current assignment to those
    // variables.
    if(vdata.factor.args() != vdata.variables) {
      assert(shared_data != NULL);
      // get all the factors
      const factor_map_t& factors = 
        *shared_data->get_constant(FACTOR_KEY).as<const factor_map_t*>();
      // Get the mrf (needed to get the assignment to variables not in
      // the clique)
      const mrf::graph_type& mrf = 
        *shared_data->get_constant(MRF_KEY).as<mrf::graph_type*>();
      // Resize the factor for the variables in the clique
      vdata.factor.set_args(vdata.variables);
      vdata.factor.uniform();

      // We now build up the factor by iteratoring over the dependent
      // factors conditioning if necessary into the conditional_factor
      // and then multiplying.
      factor_t conditional_factor;
      // Iterate over the factors and multiply each into this factor
      foreach(size_t factor_id, vdata.factor_ids) {
        const factor_t& factor = factors[factor_id];
        // Build up an assignment for the conditional
        domain_t conditional_args = factor.args() - vdata.variables;
        assignment_t conditional_asg;
        for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
          const mrf::vertex_data& mrf_vdata = 
            mrf.vertex_data(conditional_args.var(i).id);
          assert(mrf_vdata.tree_id == NULL_VID);
          conditional_asg &= mrf_vdata.asg;         
        }
        // set the factor arguments
        conditional_factor.set_args(factor.args() - conditional_args);
        conditional_factor.condition(factor, conditional_asg);        
        // Multiply the conditional factor in
        vdata.factor *= conditional_factor;
      }
      // Extra normalization for stability on the table factors
      vdata.factor.normalize();
    }




    // Determine if their are any edges for which messages can be
    // computed
    factor_t cavity; // preallocate cavity factor
    foreach(edge_id_t out_eid, scope.out_edge_ids()) {
      vertex_id_t target = scope.target(out_eid);
      edge_data& out_edata = scope.edge_data(out_eid);
      // If the message has already been computed than we don't need
      // to compute it again
      bool ready = !out_edata.calibrated;
      if(ready) {
        foreach(edge_id_t in_eid, scope.in_edge_ids()) {       
          const edge_data& in_edata = scope.const_edge_data(in_eid);
          if(scope.source(in_eid) != target && !in_edata.calibrated) {
            ready = false;
            break;
          }
        }
      } // end of if ready
      
      // Compute message if its still ready
      if(ready) {
        cavity = vdata.factor;
        foreach(edge_id_t in_eid, scope.in_edge_ids()) {
          if(scope.source(in_eid) != target) {
            const edge_data& in_edata = scope.const_edge_data(in_eid);
            cavity *= in_edata.message;    
          }
        }
        // Marginalize all variables not in outbound message
        out_edata.message.set_args(out_edata.variables);
        out_edata.message.marginalize(cavity);
        out_edata.message.normalize();
        out_edata.calibrated = true;
        // Schedule neighbor to receive message
        assert(target < scope.num_vertices());
        callback.add_task(target, calibrate_update, 1.0);
      } // if still ready
    }

    // Sample if necessary
    if(!vdata.sampled) {

      // Determine the state of the vertex and its neighbors
      bool is_calibrated = true;
      // The sampling procedure "directs" the tree and so the parent
      // is the clique that is sampled before the child.
      size_t parent_count = 0;
      vertex_id_t parent_vid = NULL_VID;
      edge_id_t parent_eid = NULL_EID;

      // Loop over the variables collecting information about the
      // neighbors
      foreach(edge_id_t in_eid, scope.in_edge_ids()) {
        const edge_data& in_edata = scope.const_edge_data(in_eid);
        // If the message has not been computed (calibrated) then we
        // are done and return
        if(!in_edata.calibrated) { 
          is_calibrated = false; 
          break; 
        } else { 
          // Otherwise proceed to test whether the neighbor has been
          // sampled
          vertex_id_t neighbor_vid = scope.source(in_eid);
          const vertex_data& neighbor_vdata = 
            scope.const_neighbor_vertex_data(neighbor_vid);
          // If the neighbor has been sampled then it is also a parent
          if(neighbor_vdata.sampled) {
            parent_count++;
            parent_eid = in_eid;
            parent_vid = neighbor_vid;
          }
        }
      } // end of loop over in edges

      vdata.calibrated = is_calibrated;


      // There can be at most one parent in the tree (the root has no
      // parents). 
      assert(parent_count < 2);

      // If this is either the "root" or has one sampled neighbor then
      // we are ready to sample
      bool is_root = scope.vertex() == 0;
      if(is_calibrated && (is_root || parent_count == 1) ) {       
        // We are ready to sample!!!  

        if(parent_count == 1) {
          assert(vdata.parent == parent_vid);
        }

        // First we determine which variables are going to be sampled
        // in this instance.  This is done by finding the parent
        // assignment if there is one
        assignment_t parent_asg;
        if(!is_root) {
          const edge_data& parent_edata = 
            scope.const_edge_data(parent_eid);
          const vertex_data& parent_vdata = 
            scope.const_neighbor_vertex_data(parent_vid);
          // Restricted the parents assignment to an assignment over
          // the edge variables
          parent_asg = parent_vdata.asg.restrict(parent_edata.variables);
        }

        // Determine the remaining variables for which we will need to
        // sample and construct RB estimates
        domain_t unsampled_variables = vdata.variables - parent_asg.args();

        // If there is nothing to sample just skip along
        if(unsampled_variables.num_vars() == 0) {
          // if there was nothing to sample just mark this clique as
          // sampled and pass along to child
          vdata.sampled = true;
          // Reschedule unsampled neighbors
          foreach(edge_id_t in_eid, scope.in_edge_ids()) {
            if(in_eid != parent_eid) {
              const vertex_id_t neighbor_vid = scope.source(in_eid);
              assert(neighbor_vid < scope.num_vertices());
              callback.add_task(neighbor_vid, 
                                calibrate_update, 
                                1.0);
            }
          }
          return;
        } // end of if variables empty

        // Now construct the full belief  
        // We start by taking the full factor
        factor_t belief = vdata.factor;
       
        // Multiply in all the messages to compute the full belief
        foreach(edge_id_t in_eid, scope.in_edge_ids()) {
          const edge_data& edata = scope.const_edge_data(in_eid);
          assert(edata.calibrated);
          belief *= edata.message;
        } // end of foreach
        belief.normalize();        

        // vdata.belief = belief;

        // Fill out the variables in the mrf
        mrf::graph_type& mrf_graph = 
          *shared_data->get_constant(MRF_KEY).as<mrf::graph_type*>();

        // First update all the RB estimates for the unsampled
        // variables in the mrf graph
        factor_t tmp_belief;
        for(size_t i = 0; i < unsampled_variables.num_vars(); ++i) {
          variable_t var = unsampled_variables.var(i);
          // Construct the RB belief estimate
          tmp_belief.set_args(var);
          tmp_belief.marginalize(belief);
          tmp_belief.normalize();
          // Update the MRF
          mrf::vertex_data& mrf_vdata = mrf_graph.vertex_data(var.id);
          mrf_vdata.belief += tmp_belief;
        } 

        // Condition the belief on the parent assignmnet
        tmp_belief.set_args(unsampled_variables);
        tmp_belief.condition(belief, parent_asg);
        tmp_belief.normalize();


        // Sample the remaining variables from the belief
        assignment_t sample_asg = tmp_belief.sample();
        // Set the local assignment
        vdata.asg = sample_asg & parent_asg;
        // the assignment should exacty cover the variables
        assert(vdata.asg.args() == vdata.variables);
        vdata.sampled = true;


        //// Fill out the MRF with the sampled variables
        for(size_t i = 0; i < sample_asg.num_vars(); ++i) {
          variable_t var = sample_asg.args().var(i);
          mrf::vertex_data& mrf_vdata = mrf_graph.vertex_data(var.id);
          mrf_vdata.asg = sample_asg.restrict(var);
          mrf_vdata.updates++;
          // std::cout << graphlab::thread::thread_id()
          //           << ": sampling " << mrf_vdata.variable << std::endl;
          // remove the vertex from any trees
          mrf_vdata.tree_id = NULL_VID;
          mrf_vdata.height = 0;
        } 

        // Reschedule unsampled neighbors
        foreach(edge_id_t in_eid, scope.in_edge_ids()) {
          if(in_eid != parent_eid) {
            const vertex_id_t neighbor_vid = scope.source(in_eid);
            assert(neighbor_vid < scope.num_vertices());
            callback.add_task(neighbor_vid, 
                              calibrate_update, 
                              1.0);
          }
        }

      } // end of if(is_calibrated and ready to sample
    } // End of if(!sampled) sampling procedure
  } // End of update function























  void fast_update(gl::iscope& scope,
                   gl::icallback& callback,
                   gl::ishared_data* shared_data) {
    typedef factorized_model::factor_map_t factor_map_t;
    
    // get the vertex data
    vertex_data& vdata = scope.vertex_data();

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
      // get all the factors
      const factor_map_t& factors = 
        *shared_data->get_constant(FACTOR_KEY).as<const factor_map_t*>();
      // Get the mrf (needed to get the assignment to variables not in
      // the clique)
      const mrf::graph_type& mrf = 
        *shared_data->get_constant(MRF_KEY).as<mrf::graph_type*>();
      // Resize the factor for the variables in the clique
      vdata.factor.set_args(vdata.variables);
      vdata.factor.uniform();

      // We now build up the factor by iteratoring over the dependent
      // factors conditioning if necessary into the conditional_factor
      // and then multiplying.
      factor_t conditional_factor;
      // Iterate over the factors and multiply each into this factor
      foreach(size_t factor_id, vdata.factor_ids) {
        const factor_t& factor = factors[factor_id];
        // Build up an assignment for the conditional
        domain_t conditional_args = factor.args() - vdata.variables;
        assignment_t conditional_asg;
        for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
          const mrf::vertex_data& mrf_vdata = 
            mrf.vertex_data(conditional_args.var(i).id);
          assert(mrf_vdata.tree_id == NULL_VID);
          conditional_asg &= mrf_vdata.asg;         
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
        edge_data& in_edata = scope.edge_data(in_eid);
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
      factor_t cavity;
      foreach(edge_id_t out_eid, scope.out_edge_ids()) {
        edge_data& out_edata = scope.edge_data(out_eid);
        edge_id_t rev_eid = scope.reverse_edge(out_eid);
        // if the out message is not calibrated try to calibrate it:
        if(!out_edata.calibrated) {
          bool ready_to_send = true;
          // check that all in messages (except the one we want to
          // send) have been recieved
          foreach(edge_id_t in_eid, scope.in_edge_ids()) {
            const edge_data& in_edata = scope.const_edge_data(in_eid);
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
            const edge_data& in_edata = scope.const_edge_data(rev_eid);
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
            callback.add_task(scope.target(out_eid), calibrate_update, 1.0);      
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
        const vertex_data& parent_vdata = 
          scope.const_neighbor_vertex_data(scope.target(out_eid));
        if(parent_vdata.sampled) {
          assert(parent_vdata.calibrated);
          assert(!parent_found);
          parent_found = true;
          to_parent_eid = out_eid;
          const edge_data& parent_edata = 
            scope.const_edge_data(to_parent_eid);
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
        // Fill out the variables in the mrf
        mrf::graph_type& mrf_graph = 
          *shared_data->get_constant(MRF_KEY).as<mrf::graph_type*>();
      
        // First update all the RB estimates for the unsampled
        // variables in the mrf graph
        factor_t tmp_belief;
        for(size_t i = 0; i < unsampled_variables.num_vars(); ++i) {
          variable_t var = unsampled_variables.var(i);
          // Construct the RB belief estimate
          tmp_belief.set_args(var);
          tmp_belief.marginalize(vdata.factor);
          tmp_belief.normalize();
          // Update the MRF
          mrf::vertex_data& mrf_vdata = mrf_graph.vertex_data(var.id);
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
          mrf::vertex_data& mrf_vdata = mrf_graph.vertex_data(var.id);
          assignment_t local_asg = sample_asg.restrict(var);
          if(mrf_vdata.asg != local_asg) {
            mrf_vdata.changes++;
            vdata.changes++;
          }
          mrf_vdata.asg = local_asg;
          mrf_vdata.updates++;
          // std::cout << graphlab::thread::thread_id()
          //           << ": sampling " << mrf_vdata.variable << std::endl;
          // remove the vertex from any trees
          mrf_vdata.tree_id = NULL_VID;
          mrf_vdata.height = 0;
          
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
                            calibrate_update, 
                            1.0);
        }
      }
    } // End of if(!sampled) sampling procedure




  } // End of update function







}; // end of namespace junction tree





#include <graphlab/macros_undef.hpp>
#endif
