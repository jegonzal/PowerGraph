#ifndef PGIBBS_MRF_HPP
#define PGIBBS_MRF_HPP




/**
 *
 * This code is ued to represent a markov random field
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>



#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cassert>


#include <graphlab.hpp>

#include "factorized_model.hpp"




#include <graphlab/macros_def.hpp>


struct mrf_vertex_data {
  //! Problem specific variables
  variable_t               variable;
  //! current assignment
  size_t                   asg;
  //! The vector of factor_ids associated with this vertex
  std::vector<factor_id_t> factor_ids;

  //! Current belief estimate
  factor_t                 belief;
  //! The number of times this vertex has been sampled
  size_t                   nsamples;
  //! The number of itmes this vertex has changed its value
  size_t                   nchanges;

  //! Properties associated with the tree
  struct tree_info_type {
    bool           in_tree; 
    vertex_id_t    tree_id;
    vertex_id_t    height;
    double         priority;
    tree_info_type () : 
      in_tree(false),
      tree_id(NULL_VID),
      height(0),
      priority(-1) { }
    void save(graphlab::oarchive& arc) const {
      arc << in_tree 
          << tree_id
          << height
          << priority;
    }
    void load(graphlab::iarchive& arc) {
      arc >> in_tree
          >> tree_id
          >> height
          >> priority;
    }
  };

  //! tree info
  tree_info_type tree_info;
  
  
  
  mrf_vertex_data() :
    asg(0),
    nsamples(0),
    nchagnes(0) { }

  mrf_vertex_data(const variable_t& variable,
                  const std::set<factor_id_t>& factor_ids_set) :
    variable(variable),
    asg(0),
    factor_ids(factor_ids_set.begin(), factor_ids_set.end()),
    belief(domain_t(variable)),
    nsamples(0),
    nchanges(0) {
    // Initialize the belief to "0"
    belief.uniform(-std::numeric_limits<double>::max());
    // Require that factor ids be non empty
    assert(!factor_ids.empty());
  }

  void save(graphlab::oarchive& arc) const {
    arc << variable
        << asg
        << factor_ids
        << belief
        << nsamples
        << nchanges
        << tree_info;
  }

  void load(graphlab::iarchive& arc) {
    arc >> variable
        >> asg
        >> factor_ids
        >> belief
        >> nsamples
        >> nchanges
        >> tree_info;
  }
  
}; // End of mrf vertex data

/**
 * The data associated with each directed edge in the pairwise markov
 * random field
 */
struct mrf_edge_data { 
  // Currently empty
  void save(graphlab::oarchive &arc) const {  }
  void load(graphlab::iarchive &arc) { }
};

// define the graph type:
typedef graphlab::graph< mrf_vertex_data, mrf_edge_data> mrf_graph_type;

typedef graplab::types<mrf_graph_type> mrg_gl;











/** Save the beliefs stored in the graph */
void save_beliefs(const mrf_graph_type& mrf,
                  const std::string& filename) {
  std::ofstream fout(filename.c_str());
  fout.precision(16);
  factor_t marginal;
  for(size_t v = 0; v < mrf.num_vertices(); ++v) {
    const vertex_data& vdata = mrf.vertex_data(v);
    marginal = vdata.belief;
    marginal.normalize();
    fout << vdata.updates << '\t';
    size_t arity = marginal.args().var(0).arity;
    for(size_t asg = 0; asg < arity; ++asg) {
      fout << std::exp( marginal.logP(asg) );
      if(asg + 1 < arity ) fout << '\t';      
    }
    fout << '\n';
  } 
  fout.close();
} // End of save beliefs












//! Compute the unormalized likelihood of the current assignment
double unnormalized_loglikelihood(const mrf_graph_type& graph,
                                  const std::vector<factor_t>& factors) {
  double sum = 0;
  // Sum the logprob of each factor
  foreach(const factor_t& factor, factors) {
    // Accumulate the assignments 
    domain_t dom = factor.args();
    assignment_t asg;
    for(size_t i = 0; i < dom.num_vars(); ++i) {
      const vertex_id_t vid = dom.var(i).id;
      const mrf::vertex_data& vdata = graph.vertex_data(vid);
      assert(vdata.variable == dom.var(i));
      asg &= assignment_t(vdata.variable, vdata.asg);
    }
    sum += factor.logP(asg);
  }
  return sum;
}












/** Construct an MRF from the factorized model */
void mrf_from_factorized_model(const factorized_model& model,
                               mrf_graph_type& mrf) {
  
  ///======================================================================
  // Add all the variables
  factor_t conditional, belief;
  foreach(variable_t variable, model.variables()) {
    mrf::vertex_data vdata(variable, model.factor_ids(variable));
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
    const vertex_id_t vid = mrf.add_vertex(vdata);
    // We require variable ids to match vertex id (this simplifies a
    // lot of stuff).
    assert(vid == variable.id);
  }  
  assert(mrf.num_vertices() == model.variables().size());

  ///======================================================================
  // Add all the edges
  const std::vector<factor_t>& factors = model.factors();
  for(vertex_id_t vid = 0; vid < mrf.num_vertices(); ++vid) {
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
      const vertex_id_t neighbor_vid = neighbor_variable.id;
      mrf::edge_data edata;
      mrf.add_edge(vid, neighbor_vid, edata);      
    }
  } // loop over factors
  mrf.finalize();
} // End of construct_mrf



#include <graphlab/macros_undef.hpp>
#endif
