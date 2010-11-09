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

enum mrf_sdt_keys {
  NSAMPLES_ID = 0,
  MAX_NSAMPLES_ID = 1,
  NUM_FACTORS_KEY = 99,
  FACTOR_OFFSET = 100,
};


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
    nchanges(0) { }

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

typedef graphlab::types<mrf_graph_type> mrf_gl;



/** Save the beliefs stored in the graph */
void save_beliefs(const mrf_graph_type& mrf,
                  const std::string& filename);



void save_asg(const mrf_graph_type& mrf,
              const std::string& filename);



//! look up a factor from the shared data table
inline const factor_t& get_factor(const mrf_gl::ishared_data& shared_data,
                                  const factor_id_t factor_id) {
  return shared_data.get_constant(FACTOR_OFFSET + factor_id).as<factor_t>();
}


//! look up a factor from the shared data table
inline size_t get_num_factors(const mrf_gl::ishared_data& shared_data) {
  return shared_data.get_constant(NUM_FACTORS_KEY).as<size_t>();
}





/** Construct an MRF from the factorized model */
void mrf_core_from_factorized_model(const factorized_model& model,
                                    mrf_gl::core& core);





//! Compute the unormalized likelihood of the current assignment
double unnormalized_loglikelihood(const mrf_gl::core& core);



void draw_mrf(const size_t experiment_id,
              const std::string& base_name, 
              const mrf_graph_type& mrf);


#endif
