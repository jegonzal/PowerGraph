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
    double         priority;
    size_t         tree_id;
    size_t         height;
    bool           in_tree; 
    tree_info_type () : 
      priority(-1), tree_id(-1), height(0), in_tree(false) { }
    void save(graphlab::oarchive& arc) const {
      arc << in_tree << tree_id << height << priority;
    }
    void load(graphlab::iarchive& arc) {
      arc >> in_tree >> tree_id >> height >> priority;
    }
  };
  //! tree info
  tree_info_type tree_info;
  mrf_vertex_data() : asg(0), nsamples(0), nchanges(0) { }
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
    ASSERT_FALSE(factor_ids.empty());
  }
  void save(graphlab::oarchive& arc) const {
    arc << variable << asg << factor_ids << belief << nsamples
        << nchanges << tree_info;
  }
  void load(graphlab::iarchive& arc) {
    arc >> variable >> asg >> factor_ids >> belief >> nsamples
        >> nchanges >> tree_info;
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

typedef graphlab::graph< mrf_vertex_data, mrf_edge_data> mrf_graph_type;

/** Save the beliefs stored in the graph */
void save_beliefs(const mrf_graph_type& mrf,
                  const std::string& filename);

void save_asg(const mrf_graph_type& mrf,
              const std::string& filename);

/** Construct an MRF from the factorized model */
void mrf_from_factorized_model(const factorized_model& model,
                               mrf_graph_type& mrf);

//! Compute the unormalized likelihood of the current assignment
double unnormalized_loglikelihood(const mrf_graph_type& mrf);

void draw_mrf(const size_t experiment_id,
              const std::string& base_name, 
              const mrf_graph_type& mrf);


#endif
