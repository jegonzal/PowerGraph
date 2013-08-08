/**  
 *  Software submitted by 
 *  Systems & Technology Research / Vision Systems Inc., 2013
 *
 *  Approved for public release; distribution is unlimited. [DISTAR Case #21428]
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


#ifndef VSI_BP_GRAPH_DATA_H
#define VSI_BP_GRAPH_DATA_H

#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include "table_factor.hpp"


namespace belief_prop {


/**
 * The type of the distributed graph representing the Factor Graph.
 */
template<size_t MAX_DIM> class vertex_data;
template<size_t MAX_DIM> class edge_data;

template<size_t MAX_DIM>
struct graph_type {
  typedef graphlab::distributed_graph<vertex_data<MAX_DIM>, edge_data<MAX_DIM> > type;
};


// Edge and Vertex data =============================================>

/**
 * The data associated with each variable in the factor graph
 */
template<size_t MAX_DIM>
class vertex_data {
  typedef graphlab::table_factor<MAX_DIM> factor_type;

public:
  double DAMPING;
  double BOUND;
  double REGULARIZATION;

  std::string name;
  // might be nice to be able to use different datatypes, however any
  // such class would be required to implement the same protocols as 
  // table_factor (and table_base if it was to use bp_vertex_program) 
  factor_type potential;
  factor_type belief;
  bool        isVariable; // is the vertex a variable or a factor 

  vertex_data() : 
      DAMPING(0.0),
      BOUND(0.0),
      REGULARIZATION(0.0),
      name(""),
      potential(factor_type::nil), 
      belief(factor_type::nil),
      isVariable(false) { } 

  vertex_data(const factor_type& potential_, 
              const factor_type& belief_, 
              bool isVariable_,
              std::string name_ = "") : 
      DAMPING(0.0), 
      BOUND(0.0), 
      REGULARIZATION(0.0),
      name(name_), 
      potential(potential_), 
      belief(belief_),
      isVariable(isVariable_) { }
  
  void load(graphlab::iarchive& arc) {
    arc >> DAMPING;
    arc >> BOUND;
    arc >> REGULARIZATION;
    arc >> name;
    arc >> potential;
    arc >> belief;
    arc >> isVariable;
  }
  void save(graphlab::oarchive& arc) const {
    arc << DAMPING;
    arc << BOUND;
    arc << REGULARIZATION;
    arc << name;
    arc << potential;
    arc << belief;
    arc << isVariable;
  }
}; // End of vertex data


/**
 * The data associated with each edge in the factor graph
 */
template<size_t MAX_DIM>
class edge_data {
  // REVIEW this could be a dense_table<1>, but not sure how operations
  // on a dense_table<16> would work
  typedef graphlab::dense_table<MAX_DIM> msg_type;

  msg_type messages[4];

  size_t message_idx(size_t source_id, size_t target_id, bool is_new) {
    return size_t(source_id < target_id) + 2 * size_t(is_new);
  }
public:
  edge_data() { 
    for(size_t i = 0; i < 4; ++i) {
      messages[i] = msg_type();
    }
  } // end of constructor

  edge_data(const msg_type& msg) {
    for(size_t i = 0; i < 4; ++i) {
      messages[i] = msg;
    }
  } // end of constructor

  msg_type& message(size_t source_id, size_t target_id) { 
    return messages[message_idx(source_id, target_id, true)];
  }
  msg_type& old_message(size_t source_id, size_t target_id) { 
    return messages[message_idx(source_id, target_id, false)];
  }
  void update_old(size_t source_id, size_t target_id) { 
    old_message(source_id, target_id) = message(source_id, target_id);
  }
  void save(graphlab::oarchive& arc) const {
    for(size_t i = 0; i < 4; ++i) arc << messages[i];
  }
  void load(graphlab::iarchive& arc) {
    for(size_t i = 0; i < 4; ++i) arc >> messages[i];
  }
}; // End of edge data


} // end of namespace belief_prop

#endif // VSI_BP_GRAPH_DATA_H
