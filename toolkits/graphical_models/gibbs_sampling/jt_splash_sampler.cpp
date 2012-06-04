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



#include "util.hpp"
#include "jt_splash_sampler.hpp"
#include "pgibbs_tls.hpp"
#include "run_statistics.hpp"
#include "global_variables.hpp"




#include <graphlab/macros_def.hpp>











void run_jtsplash_sampler(mrf_graph_type& mrf_graph,
                          const std::string& jtsplash_results_fn,
                          const std::vector<double>& runtimes,
                          const bool draw_images,
                          const splash_settings& settings) {

  //  size_t ncpus = core.engine().get_ncpus();
  // bool affinities = 
  //   core.get_engine_options().get_cpu_affinities();
  // Initialize the jtsplash sampler
  jt_splash_sampler jtsplash_sampler(mrf_graph, settings);

  double total_runtime = 0;
  double actual_total_runtime = 0;
  foreach(const double experiment_runtime, runtimes) {
    total_runtime += experiment_runtime;
    // get the experiment id
    size_t experiment_id = file_line_count(jtsplash_results_fn);
    std::cout << "Running JTSplash sampler experiment " << experiment_id
              << " for " << experiment_runtime << " seconds." << std::endl;

    std::cout << "Settings: ======================" << std::endl
              << "Experiment:    " << experiment_id << std::endl
              << "runtime:       " << experiment_runtime << std::endl
              << "treesize:      " << settings.max_tree_size << std::endl
              << "treewidth:     " << settings.max_tree_width << std::endl
              << "treeheight:    " << settings.max_tree_height << std::endl
              << "factorsize:    " << settings.max_factor_size << std::endl
              << "subthreads:    " << settings.subthreads << std::endl
              << "priorities:    " << settings.priorities << std::endl   
              << "vanish:        " << settings.vanish_updates << std::endl;   


    graphlab::timer timer;
    timer.start();

    // run the sampler once
    jtsplash_sampler.sample_seconds(experiment_runtime);

    double actual_experiment_runtime = timer.current_time();
    std::cout << "Actual Experiment Runtime:" 
              << actual_experiment_runtime << std::endl;        
    actual_total_runtime += actual_experiment_runtime;
    std::cout << "Total Experiment Runtime (actual): "
              << total_runtime 
              << "(" << actual_total_runtime << ")" 
              << std::endl;


    // check mrf graph
    for(size_t i = 0; i < mrf_graph.num_vertices(); ++i) {
      ASSERT_EQ(mrf_graph.vertex_data(i).tree_info.tree_id, NULL_VID);
    }



    /// ==================================================================
    // Compute final statistics of the mode
    run_statistics stats(mrf_graph);
    stats.print();
    // Save the beliefs
    save_beliefs(mrf_graph,
                 make_filename("jtsplash_blfs_", ".tsv", experiment_id));
    // // Save the current assignments
    save_asg(mrf_graph,
             make_filename("jtsplash_asg_", ".asg", experiment_id));
    // Save the experiment
    std::ofstream fout(jtsplash_results_fn.c_str(), std::ios::app);
    fout.precision(8);
    fout << experiment_id << '\t'
         << total_runtime << '\t'
         << actual_total_runtime << '\t'
         << settings.ntrees << '\t'
         << stats.nsamples << '\t'
         << stats.nchanges << '\t'
         << stats.loglik << '\t'
         << stats.min_samples << '\t'
         << stats.max_samples << '\t'
         << settings.max_tree_size << '\t'
         << settings.max_tree_width << '\t'
         << settings.max_factor_size << '\t'
         << settings.max_tree_height << '\t'
         << settings.subthreads << '\t'
         << settings.priorities << '\t'
         << jtsplash_sampler.total_trees() << '\t'  
         << jtsplash_sampler.total_collisions() << '\t'
         << std::endl;
    fout.close();

    // Plot images if desired
    if(draw_images) draw_mrf(experiment_id, "jtsplash", mrf_graph);

  } // end of for loop over runtimes

} // end of run_jtsplash_sampler





void jtree_update::operator()(base::icontext_type& context) {
  typedef factorized_model::factor_map_t factor_map_t;
  ASSERT_NE(mrf_ptr, NULL);
  mrf_graph_type& mrf = *mrf_ptr;
    
  // get the vertex data
  jtree_vertex_data& vdata = context.vertex_data();
    
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
    // Resize the factor for the variables in the clique
    vdata.factor.set_args(vdata.variables);
    vdata.factor.uniform();

    // We now build up the factor by iteratoring over the dependent
    // factors conditioning if necessary into the conditional_factor
    // and then multiplying.
    factor_t& conditional_factor(tls.conditional_factor);
    // Iterate over the factors and multiply each into this factor
    foreach(size_t factor_id, vdata.factor_ids) {
      //const factor_t& factor = SHARED_FACTORS.get()[factor_id];
      const factor_t& factor = (*SHARED_FACTORS_PTR)[factor_id];
      // Build up an assignment for the conditional
      domain_t conditional_args = factor.args() - vdata.variables;
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& mrf_vdata = 
          mrf.vertex_data(conditional_args.var(i).id());
        ASSERT_EQ(mrf_vdata.tree_info.tree_id, NULL_VID);
        //        ASSERT_FALSE(mrf_vdata.tree_info.in_tree);
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
    foreach(edge_id_t in_eid, context.in_edge_ids()) {
      jtree_edge_data& in_edata = context.edge_data(in_eid);
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
      received_neighbors == context.in_edge_ids().size();
  } else {
    received_neighbors = context.in_edge_ids().size();
  }


  //////////////////////////////////////////////////////////////////
  // send any unset messages 
  // if we have recieve enough in messages
  if(received_neighbors + 1 >= context.in_edge_ids().size()) {
    factor_t& cavity(tls.cavity);
    foreach(edge_id_t out_eid, context.out_edge_ids()) {
      jtree_edge_data& out_edata = context.edge_data(out_eid);
      edge_id_t rev_eid = context.reverse_edge(out_eid);
      // if the out message is not calibrated try to calibrate it:
      if(!out_edata.calibrated) {
        bool ready_to_send = true;
        // check that all in messages (except the one we want to
        // send) have been recieved
        foreach(edge_id_t in_eid, context.in_edge_ids()) {
          const jtree_edge_data& in_edata = context.const_edge_data(in_eid);
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
          const jtree_edge_data& in_edata = context.const_edge_data(rev_eid);
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
          callback.add_task(context.target(out_eid), jtree_sample_update, 1.0);      
        } // end of if ready to send
      } // end of if not calibrated
    } // end of loop over outbound messages
  } // of send all out messages


  //////////////////////////////////////////////////////////////////
  // Construct RB estimate and Sample if calibrated but not yet
  // sampled
  if(vdata.calibrated && !vdata.sampled) {
    // check that the parent is sampled and also determine which
    // variables are going to be sampled at this clique.  This is
    // done by finding the parent assignment if there is one
    assignment_t parent_asg;
    edge_id_t to_parent_eid = NULL_EID;

    // find the parent
    bool parent_found = false;
    foreach(edge_id_t out_eid, context.out_edge_ids()) {       
      const jtree_vertex_data& parent_vdata = 
        context.const_neighbor_vertex_data(context.target(out_eid));
      if(parent_vdata.sampled) {
        ASSERT_TRUE(parent_vdata.calibrated);
        ASSERT_FALSE(parent_found);
        parent_found = true;
        to_parent_eid = out_eid;
        const jtree_edge_data& parent_edata = 
          context.const_edge_data(to_parent_eid);
        parent_asg = 
          parent_vdata.asg.restrict(parent_edata.variables);
        ASSERT_EQ(parent_asg.args(), parent_edata.variables);            
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
        mrf_vertex_data& mrf_vdata = get_mrf_vdata(var.id());
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
      ASSERT_EQ(vdata.asg.args(), vdata.variables);

      
      //// Fill out the MRF with the sampled variables
      for(size_t i = 0; i < sample_asg.num_vars(); ++i) {
        variable_t var = sample_asg.args().var(i);
        mrf_vertex_data& mrf_vdata = get_mrf_vdata(var.id());
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
        //        mrf_vdata.tree_info.in_tree = false;
        mrf_vdata.tree_info.height = 0;
          
        // double& logP = mrf_vdata.belief.logP(mrf_vdata.asg.asg_at(0));
        // logP = std::log( std::exp(logP) + 1.0 );          
      } 
    } // end of sampling unsampled variables

      // mark as sampled
    vdata.sampled = true;

    // Reschedule unsampled neighbors
    foreach(edge_id_t out_eid, context.out_edge_ids()) {
      if(out_eid != to_parent_eid) {
        const vertex_id_t neighbor_vid = context.target(out_eid);
        ASSERT_LT(neighbor_vid, context.num_vertices());
        callback.add_task(neighbor_vid, 
                          jtree_sample_update, 
                          1.0);
      }
    }
  } // End of if(!sampled) sampling procedure

} // End of update function




termination_condition::termination_condition() : 
  error(false), finish_time_seconds(-1), target_nsamples(0), target_ntrees(0),
  atomic_nsamples(0), atomic_ntrees(0) { }

bool termination_condition::finished() const {
  return
    error ||
    (finish_time_seconds > 0 &&
     finish_time_seconds < graphlab::lowres_time_seconds()) ||
    (target_nsamples > 0 && atomic_nsamples.value > target_nsamples) ||
    (target_ntrees > 0 && atomic_ntrees.value > target_ntrees);
}

void termination_condition::reset() {
  error = false;
  finish_time_seconds = -1;
  target_nsamples = 0;
  atomic_nsamples.value = 0;
  target_ntrees = 0;
  atomic_ntrees.value = 0;
}










jt_worker::jt_worker(size_t worker_id, 
                     const splash_settings& settings,
                     scope_factory_type& scope_factory, 
                     const std::vector<vertex_id_t>& root_perm,
                     termination_condition& terminator) :
  worker_id(worker_id),
  settings(settings),
  scope_factory_ptr(&scope_factory),
  root_index(root_perm.size()),
  root_perm_ptr(&root_perm),
  current_root(root_perm.at(worker_id)),
  terminator_ptr(&terminator),
  ncollisions(0) {  
  // Initialize local jtcore
  if(settings.subthreads > 1) {
    jt_core.set_scheduler_type("multiqueue_fifo");
    jt_core.set_scope_type("edge");
    jt_core.set_ncpus(settings.subthreads);
    jt_core.set_engine_type("async");
  } else {
    jt_core.set_scheduler_type("fifo");
    jt_core.set_scope_type("none");
    jt_core.set_ncpus(1);
    jt_core.set_engine_type("async_sim");
  }
} // end of init








// get a root
void jt_worker::run() {   
  // looup until runtime is reached
  while(!terminator_ptr->finished()) {
    /////////////////////////////////////////////////////////
    // Construct one tree (we must succeed in order to count a tree)
    advance_root();
    //    std::cout << "Worker " << worker_id << ": " << current_root << std::endl;
    // here we loop until the current root is sampled
    size_t sampled_variables = 0;
    while(sampled_variables == 0 &&  !terminator_ptr->finished()) {
      sampled_variables = splash_once();
      // If sample once failed due to collision record a collision event
      if(sampled_variables == 0) {
        ncollisions++;
        // sched_yield();
      }
    }
    // if variables where sampled in the splash increment the atomic
    // counters.
    if(sampled_variables > 0) {
      terminator_ptr->atomic_nsamples.inc(sampled_variables);
      terminator_ptr->atomic_ntrees.inc();
    }
  } 
  //  std::cout << "N Collisions: " << ncollisions << std::endl;
} // end of run




void jt_worker::advance_root() {  
  root_index += settings.ntrees;
  if(root_index >= root_perm_ptr->size()) root_index = worker_id;
  current_root = root_perm_ptr->at(root_index);
}




///////////////////////////////////////////////////////////////////////
/// Markov Blanket Locking Helper functions


/**
 * See if the vertex can be grabbed into this workers tree. If true we
 * still need to actually grab the vertex (which could still
 * fail). However if the vertex is currently unavailable we could save
 * time by not even trying (although it may become available later).
 */
bool jt_worker::is_vertex_available(vertex_id_t vid) {
  ASSERT_NE(scope_factory_ptr, NULL);
  const mrf_graph_type& mrf(scope_factory_ptr->get_graph());
  const mrf_vertex_data& vdata = mrf.vertex_data(vid);
  // Check that this vertex is not already in a tree
  bool in_tree = vdata.tree_info.tree_id != NULL_VID;
  if(in_tree) return false;
  // check that the neighbors are not in any other trees than this
  // one
  const mrf_gl::edge_list& in_eids = mrf.in_edge_ids(vid);
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
 * Grab this vertex into the tree owned by worker id.  If this returns
 * true than the vertex is marked as grabbed. This must be called
 * within the context of an edge scope.
 */
bool jt_worker::try_grab_vertex(iscope_type& scope) {
  // Check that this vertex is not already in a tree
  bool in_tree = scope.vertex_data().tree_info.tree_id != NULL_VID;
  if(in_tree) return false;

  // check that the neighbors are not in any other trees than this
  // one
  vertex_id_t min_height(std::numeric_limits<vertex_id_t>::max());
  foreach(edge_id_t in_eid, scope.in_edge_ids()) {
    vertex_id_t neighbor_vid = scope.source(in_eid);
    const mrf_vertex_data& vdata = 
      scope.const_neighbor_vertex_data(neighbor_vid);
    bool in_tree = vdata.tree_info.tree_id != NULL_VID;
    // if the neighbor is in a tree other than this one quit
    if(in_tree && worker_id != vdata.tree_info.tree_id) return false;
    if(in_tree) min_height = std::min(min_height, vdata.tree_info.height);
  }
  // Assert that this vertex is not in a tree and that none of the
  // neighbors are in other trees This vertex does not neighbor any
  // other trees than this one
  scope.vertex_data().tree_info.tree_id = worker_id;
  //  scope.vertex_data().tree_info.in_tree = true;
  scope.vertex_data().tree_info.height = min_height + 1;
  return true;
} // end of try grab vertex




/**
 * Release the vertex
 */
void jt_worker::release_vertex(iscope_type& scope) {
  // This vertex does not neighbor any other trees than this one
  scope.vertex_data().tree_info.tree_id = NULL_VID;
  //  scope.vertex_data().tree_info.in_tree = false;
  scope.vertex_data().tree_info.height = 0;
} // release the vertex











///////////////////////////////////////////////////////////////////////
/// Scoring helper functions

/**
 * This function returns the priority of a particular vertex.
 */
double jt_worker::score_vertex(vertex_id_t vid) {
  ASSERT_NE(scope_factory_ptr, NULL);
  mrf_graph_type& mrf(scope_factory_ptr->get_graph());
  mrf_vertex_data& vdata = mrf.vertex_data(vid);

  if (vdata.nsamples < settings.vanish_updates || 
      vdata.tree_info.priority < 0) {
    vdata.tree_info.priority = score_vertex_log_odds(vid); 
  }
  return vdata.tree_info.priority;
}






double jt_worker::score_vertex_l1_diff(vertex_id_t vid) {
  // Get the scope factory
  ASSERT_NE(scope_factory_ptr, NULL);
  const mrf_graph_type& mrf(scope_factory_ptr->get_graph());
  const mrf_vertex_data& vdata = mrf.vertex_data(vid);

  // Construct the domain of neighbors that are already in the tree
  domain_t vars = vdata.variable;
  foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
    const vertex_id_t neighbor_vid = mrf.source(ineid);
    const mrf_vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
    // test to see if the neighbor is in the tree by checking the
    // elimination time map
    if(jt_list.contains(neighbor_vid)) {
      vars += neighbor.variable;
      // If this vertex has too many tree neighbor than the priority
      // is set to -1;
      if(vars.num_vars() > settings.max_tree_width + 1) return -1;
      if(vars.size() > settings.max_factor_size) return -1;
    } 
  }

  // Compute the clique factor
  clique_factor.set_args(vars);
  clique_factor.uniform();
  // get all the factors const factorized_model::factor_map_t&
  // factors(SHARED_FACTORS.get());
  const factorized_model::factor_map_t& factors(*SHARED_FACTORS_PTR);
  // Iterate over the factors and multiply each into this factor
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor = factors[factor_id];      
    // Build up an assignment for the conditional
    domain_t conditional_args = factor.args() - vars;
    if(conditional_args.num_vars() > 0) {
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& neighbor_vdata = 
          mrf.vertex_data(conditional_args.var(i).id());
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

  ASSERT_GE( residual, 0);
  ASSERT_FALSE( std::isnan(residual) );
  ASSERT_TRUE( std::isfinite(residual) );

  // ensure score is bounded
  //    residual = std::tanh(residual);

  return residual;

} // end of score l1 diff








double jt_worker::score_vertex_log_odds(vertex_id_t vid) {
  // Get the scope factory
  const mrf_graph_type& mrf(scope_factory_ptr->get_graph());
  const mrf_vertex_data& vdata(mrf.vertex_data(vid));

  // Construct the domain of neighbors that are already in the tree
  domain_t vars = vdata.variable;
  foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
    const vertex_id_t neighbor_vid = mrf.source(ineid);
    const mrf_vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
    // test to see if the neighbor is in the tree by checking the
    // elimination time map
    if(jt_list.contains(neighbor_vid)) {
      vars += neighbor.variable;
      // If this vertex has too many tree neighbor than the priority
      // is set to 0;
      if(vars.num_vars() > settings.max_tree_width + 1) return -1;
      if(vars.size() > settings.max_factor_size) return -1;
    } 
  }
    
  ASSERT_EQ(vars.num_vars(),  2);


  // Compute the clique factor
  clique_factor.set_args(vars);
  clique_factor.uniform();
  // get all the factors
  // const factorized_model::factor_map_t& factors(SHARED_FACTORS.get());
  const factorized_model::factor_map_t& factors(*SHARED_FACTORS_PTR);

  // Iterate over the factors and multiply each into this factor
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor = factors[factor_id];      
    // Build up an assignment for the conditional
    domain_t conditional_args = factor.args() - vars;
    if(conditional_args.num_vars() > 0) {
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& neighbor_vdata = 
          mrf.vertex_data(conditional_args.var(i).id());
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

  ASSERT_GE( residual, 0);
  ASSERT_FALSE( std::isnan(residual) );
  ASSERT_TRUE( std::isfinite(residual) );

  // ensure score is bounded
  //    residual = std::tanh(residual);

  return residual;
} // end of score vertex















double jt_worker::score_vertex_lik(vertex_id_t vid) {
  // Get the scope factory
  const mrf_graph_type& mrf(scope_factory_ptr->get_graph());
  const mrf_vertex_data& vdata(mrf.vertex_data(vid));

  // Construct the domain of neighbors that are already in the tree
  domain_t vars = vdata.variable;
  foreach(edge_id_t ineid, mrf.in_edge_ids(vid)) {
    const vertex_id_t neighbor_vid = mrf.source(ineid);
    const mrf_vertex_data& neighbor = mrf.vertex_data(neighbor_vid);
    // test to see if the neighbor is in the tree by checking the
    // elimination time map
    if(jt_list.contains(neighbor_vid)) {
      vars += neighbor.variable;
      // If this vertex has too many tree neighbor than the priority
      // is set to 0;
      if(vars.num_vars() > settings.max_tree_width + 1) return -1;
      if(vars.size() > settings.max_factor_size) return -1;
    } 
  }
    
  // Compute the clique factor
  clique_factor.set_args(vars);
  clique_factor.uniform();
  // get all the factors
  // const factorized_model::factor_map_t& factors(SHARED_FACTORS.get());
  const factorized_model::factor_map_t& factors(*SHARED_FACTORS_PTR);

  // Iterate over the factors and multiply each into this factor
  foreach(size_t factor_id, vdata.factor_ids) {
    const factor_t& factor = factors[factor_id];      
    // Build up an assignment for the conditional
    domain_t conditional_args = factor.args() - vars;
    if(conditional_args.num_vars() > 0) {
      assignment_t conditional_asg;
      for(size_t i = 0; i < conditional_args.num_vars(); ++i) {
        const mrf_vertex_data& neighbor_vdata = 
          mrf.vertex_data(conditional_args.var(i).id());
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

  ASSERT_GE( residual, 0);
  ASSERT_FALSE( std::isnan(residual) );
  ASSERT_TRUE(  std::isfinite(residual) );

  // // ensure score is bounded
  // residual = std::tanh(residual);


  return residual;
} // end of max lik





















///////////////////////////////////////////////////////////////////////
/// Tree Growing helper functions


void jt_worker::grow_bfs_jtree() {
  ASSERT_NE(scope_factory_ptr, NULL);
  // Get the scope factory
  mrf_graph_type& mrf = scope_factory_ptr->get_graph();
  // Clear local data structures
  jt_list.clear();
  bfs_queue.clear();
  visited.clear();
     
  // add the root
  bfs_queue.push_back(current_root);
  visited.insert(current_root);

  while(!bfs_queue.empty() && 
        jt_list.cliques.size() < settings.max_tree_size) {
    // Take the top element
    const vertex_id_t next_vertex = bfs_queue.front();
    bfs_queue.pop_front();

    // pretest that the vertex is available before trying to get it
    if(!is_vertex_available(next_vertex)) continue;

    // Maybe we can get the vertex so actually try to get it by first
    // grabbing the lock (scope) for the vertex
    iscope_type* scope_ptr = 
      scope_factory_ptr->get_edge_scope(worker_id, next_vertex);
    ASSERT_NE(scope_ptr, NULL);
    iscope_type& scope(*scope_ptr);

    // See if we can get the vertex for this tree
    if(!try_grab_vertex(scope)) {
      // release the scope and move on
      scope_factory_ptr->release_scope(&scope);        
      continue;
    }

    // Assert that we own the vertex at this point
    ASSERT_EQ(scope.vertex_data().tree_info.tree_id, worker_id);
    
    // Determine if this is the root (it is the root if there are no
    // cliques in the junction tree).
    bool is_root = jt_list.cliques.empty();   
    // Set the height of the root to be zero explicity
    if(is_root) scope.vertex_data().tree_info.height = 0;
 
    // Check if it is safe to extend to the tree to include next variable
    bool extended_jtree =
      scope.vertex_data().tree_info.height < settings.max_tree_height
      &&
      jt_list.extend(next_vertex,
                     mrf,
                     settings.max_tree_width,
                     settings.max_factor_size);

    // If we were unable to extend the tree then release the vertex
    if(!extended_jtree) {
      release_vertex(scope);     
      scope_factory_ptr->release_scope(&scope);        
      continue;
    } 


    // add the neighbors to the search queue
    foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
      vertex_id_t neighbor_vid = mrf.target(eid);
      if(visited.count(neighbor_vid) == 0) {
        bfs_queue.push_back(neighbor_vid);
        visited.insert(neighbor_vid);
      }
    }    

    // Release the scope and let neighbors start to run
    scope_factory_ptr->release_scope(&scope);        
  } // end of while loop

} // end grow_bfs_jtree







void jt_worker::grow_prioritized_jtree() {
  ASSERT_NE(scope_factory_ptr, NULL);
  // Get the scope factory
  mrf_graph_type& mrf = scope_factory_ptr->get_graph();
  // Clear local data structures
  jt_list.clear();
  priority_queue.clear();
  visited.clear();
     
  // add the root
  priority_queue.push(current_root, 1.0);
  visited.insert(current_root);

  while(!priority_queue.empty() && 
        jt_list.cliques.size() < settings.max_tree_size) {
    // Take the top element
    const vertex_id_t next_vertex = priority_queue.pop().first;

    // pretest that the vertex is available before trying to get it
    if(!is_vertex_available(next_vertex)) continue;

    // Maybe we can get the vertex so actually try to get it by first
    // grabbing the lock (scope) for the vertex
    iscope_type* scope_ptr = 
      scope_factory_ptr->get_edge_scope(worker_id, next_vertex);
    ASSERT_NE(scope_ptr, NULL);
    iscope_type& scope(*scope_ptr);

    // See if we can get the vertex for this tree
    if(!try_grab_vertex(scope)) {
      // release the scope and move on
      scope_factory_ptr->release_scope(&scope);        
      continue;
    }

    // Assert that we own the vertex at this point
    ASSERT_EQ(scope.vertex_data().tree_info.tree_id, worker_id);

    // Determine if this is the root (it is the root if there are no
    // cliques in the junction tree).
    bool is_root = jt_list.cliques.empty();   
    // Set the height of the root to be zero explicity
    if(is_root) scope.vertex_data().tree_info.height = 0;
          
    // test the 
    bool extended_jtree =
      scope.vertex_data().tree_info.height < settings.max_tree_height
      &&
      jt_list.extend(next_vertex,
                     mrf,
                     settings.max_tree_width,
                     settings.max_factor_size);                    

    // If we were unable to extend the tree then release the vertex
    if(!extended_jtree) {
      release_vertex(scope);
      scope_factory_ptr->release_scope(&scope);
      continue;
    }

    // If the tree was extended, extend the boundary by adding the
    // neighbors of the newly added vertex to the tree

    // add the neighbors to the search queue or update their priority
    foreach(edge_id_t eid, mrf.out_edge_ids(next_vertex)) {
      vertex_id_t neighbor_vid = mrf.target(eid);          
      if(visited.count(neighbor_vid) == 0) {
        visited.insert(neighbor_vid);
        // Vertex has not yet been visited
        double score = score_vertex(neighbor_vid);
        // if the score is greater than zero then add the neighbor
        // to the priority queue.  The score is zero if there is no
        // advantage or the treewidth is already too large
        if(score >= 0) priority_queue.push(neighbor_vid, score);       
      } else if(priority_queue.contains(neighbor_vid)) {
        // vertex is still in queue we may need to recompute
        // score
        double score = score_vertex(neighbor_vid);
        if(score >= 0) {
          // update the priority queue with the new score
          priority_queue.update(neighbor_vid, score);
        } else {
          // The score computation revealed that the clique would be
          // too large so simply remove the vertex from the priority
          // queue
          priority_queue.remove(neighbor_vid);
        }
      } // otherwise the vertex has been visited and processed
    } // end of foreach    

    // Release the scope and let neighbors start to run
    scope_factory_ptr->release_scope(&scope);        

  } // end of while loop
    
} // end grow_prioritized_jtree













size_t jt_worker::splash_once() {
  if(settings.priorities) {
    // grow the prioritized junction tree data structure
    grow_prioritized_jtree();
  } else {
    // grow the bfs junction tree data structure
    grow_bfs_jtree();
  }

  ASSERT_NE(scope_factory_ptr, NULL);
  // Get the scope factory
  mrf_graph_type& mrf = scope_factory_ptr->get_graph();
    
  // If we failed to build a tree return failure
  if(jt_list.cliques.empty()) return 0;

  //  std::cout << "Varcount: " << jt_list.cliques.size() << std::endl;  

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
  size_t num_factors = (*SHARED_FACTORS_PTR).size();

  // jt_list.validate();

  jt_list.load_graph(mrf, num_factors, jt_core.graph());

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
    ASSERT_TRUE(vdata.sampled);
    ASSERT_TRUE(vdata.calibrated);
    ASSERT_GT(vdata.variables.num_vars(), 0);
    actual_tree_width = 
      std::max(vdata.variables.num_vars() - 1, actual_tree_width); 
  }    
  //  std::cout << "Treewidth: " << actual_tree_width << std::endl;

  // Return the number of variables in the tree
  return jt_list.elim_time.size();
} // end of sample once











jt_splash_sampler::
jt_splash_sampler(mrf_graph_type& mrf, 
                  const splash_settings& settings) :
  workers(settings.ntrees, NULL),
  scope_factory(mrf, settings.ntrees, 
                graphlab::scope_range::EDGE_CONSISTENCY),
  root_perm(mrf.num_vertices()) { 
  ASSERT_LE(settings.ntrees, mrf.num_vertices());
  
  // Set the shared graph pointer
  shared_mrf_ptr = &mrf;

  // Shuffle the root ordering 
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid)
    root_perm[vid] = vid;
  graphlab::random::shuffle(root_perm.begin(), root_perm.end());
  
  // initialize the worker thread objects
  for(size_t i = 0; i < workers.size(); ++i) {
    workers[i] = 
      new jt_worker(i, settings, scope_factory, root_perm, terminator);
  }
} // end of constructor

jt_splash_sampler::~jt_splash_sampler() {
  for(size_t i = 0; i < workers.size(); ++i) {
    if(workers[i] != NULL) {
      delete workers[i];
      workers[i] = NULL;
    }
  }
}

size_t jt_splash_sampler::total_collisions() const {
  size_t total_collisions = 0;
  foreach(const jt_worker* worker, workers)  {
    ASSERT_NE(worker, NULL);
    total_collisions += worker->ncollisions;
  }
  return total_collisions;
}


size_t jt_splash_sampler::total_trees() const {
  return terminator.atomic_ntrees.value;
}

size_t jt_splash_sampler::total_samples() const {
  return terminator.atomic_nsamples.value;
}

  
void jt_splash_sampler::sample_seconds(float runtime_secs) {
  // Set the termination condition
  terminator.reset();
  terminator.finish_time_seconds = 
    graphlab::lowres_time_seconds() + runtime_secs;
  run();
}                   

void jt_splash_sampler::sample_trees(size_t total_trees) {
  // Set the termination condition
  terminator.reset();
  terminator.target_ntrees = total_trees;
  run();
}

void jt_splash_sampler::sample_updates(size_t total_updates) {
  // Set the termination condition
  terminator.reset();
  terminator.target_nsamples = total_updates;
  run();
}








void jt_splash_sampler::run() {
  // create worker threads
  graphlab::thread_group threads;
  if(workers.size() == 1) {
    ASSERT_NE(workers[0], NULL);
    workers[0]->run();
  } else {
    // Launch the threads
    for(size_t i = 0; i < workers.size(); ++i) {   
      ASSERT_NE(workers[i], NULL);
      // if(use_cpu_affinity) 
      //   threads.launch(boost::bind(&jt_worker::run, &(workers[i])), i);
      // else 
      threads.launch(boost::bind(&jt_worker::run, workers[i]));
    }
    const char* exception_message = "Exception!";
    // Wait for all threads to finish
    while (threads.running_threads() > 0) {
      try {
        threads.join();
      } catch(const char* error) {
        logstream(LOG_ERROR) << "Exception Caught:\n\t" << error << std::endl;
        exception_message = error;
        // killall the running threads
        terminator.error = true;
      }
    }
    if(terminator.error) {
      throw exception_message;
    }
  }
}







