#ifndef SINGLE_GIBBS_UPDATE_FUNCTION_HPP
#define SINGLE_GIBBS_UPDATE_FUNCTION_HPP



#include "factorized_model.hpp"
#include "mrf.hpp"

// Include the macro for the foreach operation
#include <graphlab/macros_def.hpp>





void single_gibbs_update(mrf_gl::iscope& scope, 
                         mrf_gl::icallback& scheduler,
                         mrf_gl::ishared_data* shared_data) {
  assert(shared_data != NULL);
  mrf_vertex_data& vdata = scope.vertex_data();

  //TODO: switch to use tls
  factor_t belief(vdata.variable);
  belief.uniform();
  foreach(const factor_id_t factor_id, vdata.factor_ids) {
    const factor_t& factor =
      shared_data->get_constant(FACTOR_ID + factor_id).as<factor_t>();
    // build the conditional
    assignment_t conditional_asg = factor.args() - vdata.variable;
    for(size_t i = 0; i < conditional_asg.num_vars(); ++i) {
      const mrf_vertex_data& other_vdata = 
	scope.const_neighbor_vertex_data(conditional_asg.args().var(i).id);
      assert(conditional_asg.args().var(i) == other_vdata.variable);
      conditional_asg &= 
	assignment_t(other_vdata.variable, other_vdata.asg);
    }
    belief.times_condition(factor, conditional_asg);
    belief.normalize();
  }
  size_t new_asg = belief.sample().asg_at(0);
  vdata.nchanges += (new_asg != vdata.asg);
  vdata.asg = new_asg;
  vdata.belief += belief;
  vdata.nsamples++;
  
}




/** Get the update counts for a vertex */
size_t get_nsamples(const mrf_vertex_data& vdata) {
  return vdata.nsamples;
}



bool nsamples_terminator(const mrf_gl::ishared_data* shared_data) {
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
