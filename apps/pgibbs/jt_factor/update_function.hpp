#ifndef PGIBBS_JT_UPDATE_FUNCTION_HPP
#define PGIBBS_JT_UPDATE_FUNCTION_HPP

#include <graphlab.hpp>


#include "data_structures.hpp"




#include <graphlab/macros_def.hpp>



namespace junction_tree{
  enum SDT_KEYS {FACTOR_KEY, MRF_KEY};


  void calibrate_update(gl::iscope& scope,
                        gl::icallback& callback,
                        gl::ishared_data* shared_data) {
    typedef factorized_model::factor_map_t factor_map_t;
    
    // get the vertex data
    vertex_data vdata = scope.vertex_data();
        

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
        for(size_t i = 0; i < conditional_args.num_vars(); ++i)
          conditional_asg &= 
            mrf.vertex_data(conditional_args.var(i).id).asg;
        // set the factor arguments
        conditional_factor.set_args(factor.args() - conditional_args);
        conditional_factor.condition(factor, conditional_asg);        
        // Multiply the conditional factor in
        vdata.factor *= conditional_factor;
      }
    }

    // Determine if their are any edges for which messages can be
    // computed
    foreach(edge_id_t out_eid, scope.out_edge_ids()) {
      vertex_id_t target = scope.target(out_eid);
      edge_data& out_edata = scope.edge_data(out_eid);
      // If the message has already been computed than we don't need
      // to compute it again
      bool ready = (out_edata.message.args() != out_edata.variables);      
      if(ready) {
        foreach(edge_id_t in_eid, scope.in_edge_ids()) {       
          const edge_data& in_edata = scope.const_edge_data(in_eid);
          ready =  // target edge
            (scope.source(in_eid) == target) ||
            // calibrated
            (in_edata.message.args() == in_edata.variables);
          if(!ready) break;
        }
      } // end of if ready
      
      // Compute message if its still ready
      if(ready) {
        factor_t belief_factor = vdata.factor;
        foreach(edge_id_t in_eid, scope.in_edge_ids()) {
          if(scope.source(in_eid) != target) {
            const edge_data& in_edata = scope.const_edge_data(in_eid);
            belief_factor *= in_edata.message;    
          }
        }
        // Marginalize all variables not in outbound message
        out_edata.message.set_args(out_edata.variables);
        out_edata.message.marginalize(belief_factor);
        callback.add_task(target, calibrate_update, 1.0);
      } // if still ready
    }
  }



}; // end of namespace junction tree





#include <graphlab/macros_undef.hpp>
#endif
