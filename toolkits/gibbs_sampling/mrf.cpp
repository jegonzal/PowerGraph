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



#include "mrf.hpp"

#include "util.hpp"
#include "image.hpp"

#include "global_variables.hpp"


#include <graphlab/macros_def.hpp>





/** Save the beliefs stored in the graph */
void save_beliefs(const mrf_graph_type& mrf,
                  const std::string& filename) {
  std::ofstream fout(filename.c_str());
  fout.precision(16);
  factor_t marginal;
  for(size_t v = 0; v < mrf.num_vertices(); ++v) {
    const mrf_vertex_data& vdata = mrf.vertex_data(v);
    marginal = vdata.belief;
    marginal.normalize();
    fout << vdata.nsamples << '\t';
    size_t arity = marginal.args().var(0).size();
    for(size_t asg = 0; asg < arity; ++asg) {
      fout << std::exp( marginal.logP(asg) );
      if(asg + 1 < arity ) fout << '\t';      
    }
    fout << '\n';
  } 
  fout.close();
} // End of save beliefs










void save_asg(const mrf_graph_type& mrf,
              const std::string& filename) {
  std::ofstream fout(filename.c_str());
  for(size_t v = 0; v < mrf.num_vertices(); ++v) 
    fout << mrf.vertex_data(v).asg << '\n';
  fout.close();
} // End of save beliefs









/** Construct an MRF from the factorized model */
void mrf_from_factorized_model(const factorized_model& model,
                               mrf_graph_type& mrf) {
  typedef mrf_graph_type::vertex_id_type vertex_id_type;
  typedef mrf_graph_type::edge_id_type   edge_id_type;
  ///======================================================================
  // Add all the variables
  factor_t conditional, belief;
  foreach(variable_t variable, model.variables()) {
    mrf_vertex_data vdata(variable, model.factor_ids(variable));
    { // Construct a uniformly random initial assignment
      assignment_t asg(vdata.variable);
      asg.uniform_sample();
      vdata.asg = asg.asg_at(0);
      double& logP = vdata.belief.logP(vdata.asg);
      logP = log(exp(logP) + 1.0);
    }
    // { // construct mode center initial assignment
    //   belief.set_args(variable);
    //   belief.uniform();
    //   conditional.set_args(variable);
    //   const std::set<vertex_id_t>& factor_ids = model.factor_ids(variable);
    //   foreach(vertex_id_t fid, factor_ids) {
    //     conditional.marginalize(model.factors()[fid]);
    // 	   belief *= conditional;
    //   }
    //   belief.normalize();
    //   assignment_t asg = belief.sample();
    //   vdata.asg = asg.asg_at(0);
    //   double& logP = vdata.belief.logP(vdata.asg);
    //   logP = log(exp(logP) + 1.0);
    // }
    const vertex_id_type vid = mrf.add_vertex(vdata);
    // We require variable ids to match vertex id (this simplifies a
    // lot of stuff).
    ASSERT_EQ(vid, variable.id());
  }  
  ASSERT_EQ(mrf.num_vertices(), model.variables().size());

  ///======================================================================
  // Add all the edges
  const factorized_model::factor_map_t& factors(model.factors());
  for(vertex_id_type vid = 0; vid < mrf.num_vertices(); ++vid) {
    const mrf_vertex_data& vdata = mrf.vertex_data(vid);
    // Compute all the neighbors of this vertex by looping over all
    // the variables in all the factors that contain this vertex
    std::set<variable_t> neighbors;
    foreach(const factor_id_t fid, vdata.factor_ids) {
      const domain_t& args = factors[fid].args();
      for(size_t n = 0; n < args.num_vars(); ++n) {
        variable_t neighbor_var = args.var(n);
        if(vdata.variable != neighbor_var )
          neighbors.insert(neighbor_var);
      }
    }
    // For each of those variables add an edge from this varaible to
    // that variable
    foreach(const variable_t neighbor_variable, neighbors) {
      const vertex_id_type neighbor_vid = neighbor_variable.id();
      mrf_edge_data edata;
      mrf.add_edge(vid, neighbor_vid, edata);      
    }
  } // loop over factors
  mrf.finalize();
} // End of construct_mrf













//! Compute the unormalized likelihood of the current assignment
double unnormalized_loglikelihood(const mrf_graph_type& mrf) {
  typedef mrf_graph_type::vertex_id_type vertex_id_type;
  double sum = 0;
  //  size_t num_factors = SHARED_FACTORS.get().size();
  size_t num_factors = SHARED_FACTORS_PTR->size();
  // Sum the logprob of each factor
  for(factor_id_t fid = 0; fid < num_factors; ++fid) {
    // const factor_t& factor(SHARED_FACTORS.get()[fid]);
    const factor_t& factor((*SHARED_FACTORS_PTR)[fid]);
    // Accumulate the assignments 
    domain_t dom = factor.args();
    assignment_t asg;
    for(size_t i = 0; i < dom.num_vars(); ++i) {
      const vertex_id_type vid = dom.var(i).id();
      const mrf_vertex_data& vdata = mrf.vertex_data(vid);
      ASSERT_EQ(vdata.variable, dom.var(i));
      asg &= assignment_t(vdata.variable, vdata.asg);
    }
    sum += factor.logP(asg);
  }
  return sum;
}













void draw_mrf(const size_t experiment_id,
              const std::string& base_name, 
              const mrf_graph_type& mrf) {
  typedef mrf_graph_type::vertex_id_type vertex_id_type;
  size_t rows = std::sqrt(mrf.num_vertices());
  std::cout << "Rows: " << rows << std::endl;
  image img(rows, rows);
  std::vector<double> values(1);
  factor_t belief;
  for(vertex_id_type vid = 0; vid < mrf.num_vertices(); ++vid) {
    const mrf_vertex_data& vdata = mrf.vertex_data(vid);
    belief = vdata.belief;
    belief.normalize();
    belief.expectation(values);
    img.pixel(vid) = values[0];
  }
  img.pixel(0) = 0;
  img.pixel(1) = mrf.vertex_data(0).variable.size()-1;
  img.save(make_filename(base_name + "_pred_", ".pgm", experiment_id).c_str());
  
  for(vertex_id_type vid = 0; vid < mrf.num_vertices(); ++vid) {   
    img.pixel(vid) = mrf.vertex_data(vid).nsamples;
  }
  img.save(make_filename(base_name + "_updates_", ".pgm", experiment_id).c_str());
  
  for(vertex_id_type vid = 0; vid < mrf.num_vertices(); ++vid) {   
    img.pixel(vid) = mrf.vertex_data(vid).nsamples == 0;
  }
  img.save(make_filename(base_name + "_unsampled_", ".pgm", experiment_id).c_str());
  
  for(vertex_id_type vid = 0; vid < mrf.num_vertices(); ++vid) {   
    img.pixel(vid) = mrf.vertex_data(vid).asg;
  }
  img.pixel(0) = 0;
  img.pixel(1) = mrf.vertex_data(0).variable.size()-1;
  img.save(make_filename(base_name + "_final_sample_", ".pgm", experiment_id).c_str());
} // end of draw_mrf

#include <graphlab/macros_undef.hpp>
