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


#ifndef PGIBBS_JUNCTION_TREE_HPP
#define PGIBBS_JUNCTION_TREE_HPP


/**
 *
 * Represents a junction tree
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


#include <boost/unordered_map.hpp>


#include <graphlab/logger/assertions.hpp>
#include <graphlab.hpp>


#include "factorized_model.hpp"
#include "mrf.hpp"




struct jtree_vertex_data; 
struct jtree_edge_data;

typedef graphlab::graph< jtree_vertex_data, jtree_edge_data> jtree_graph_type;


struct jtree_vertex_data {
  jtree_graph_type::vertex_id_type parent;
  domain_t variables;
  bool calibrated;
  bool sampled;
  std::vector<factor_id_t> factor_ids;
  factor_t factor;
  assignment_t asg;
  
  jtree_vertex_data() : parent(-1), calibrated(false), sampled(false) { }
}; // End of vertex data


struct jtree_edge_data {
  domain_t variables;
  factor_t message;
  bool calibrated;
  bool received;
  jtree_edge_data() : 
    calibrated(false), 
    received(false) { }
}; // End of edge data



//// Junction tree construction code
//// =====================================================================>

//! The fast set used in junction tree construction
typedef graphlab::small_set<2*MAX_DIM, jtree_graph_type::vertex_id_type> 
vertex_set;




struct jtree_list {
  struct elim_clique {
    //! The parent of this elim clique in the jtree_list
    jtree_graph_type::vertex_id_type parent;
    //! The vertex eliminated when this clique was created 
    mrf_graph_type::vertex_id_type elim_vertex;
    //! The vertices created in this clique EXCLUDING elim_vertex
    vertex_set vertices; 
    elim_clique() : parent(-1) { }
  };
  typedef std::vector<elim_clique> clique_list_type;
  typedef boost::unordered_map<mrf_graph_type::vertex_id_type,
                               mrf_graph_type::vertex_id_type> 
  elim_time_type;
  
  //! The collection of cliques
  clique_list_type cliques;
  //! the time variable i was eliminated
  elim_time_type   elim_time;

  inline bool contains(const mrf_graph_type::vertex_id_type vid) const {
    return elim_time.find(vid) != elim_time.end();
  }

  inline mrf_graph_type::vertex_id_type
  elim_time_lookup(const mrf_graph_type::vertex_id_type vid) const {
    elim_time_type::const_iterator iter(elim_time.find(vid));
    ASSERT_TRUE(iter != elim_time.end());
    return iter->second;
  }


  inline void clear() {
    cliques.clear();
    elim_time.clear();
  }

  /**
   * Extend the jtree_list data structure by eliminating the vertex.
   * If the jtree list can be extended then it is extended and this
   * function returns true.
   *
   **/
  bool extend(const mrf_graph_type::vertex_id_type elim_vertex,
              const mrf_graph_type& mrf,
              const size_t max_tree_width,
              const size_t max_factor_size);

  /**
   * Convert a jtree_list into a jtree_graph 
   */
  void load_graph(const mrf_graph_type& mrf,
                  const size_t num_factors,
                  jtree_graph_type& jt_graph) const;

  /**
   * Check internal data structures
   */
  void validate() const;

};






// /**
//  * Extend the jtree_list data structure by eliminating the vertex.  If
//  * the jtree list can be extended then it is extended and this
//  * function returns true.
//  *
//  **/
// bool extend_jtree_list(const vertex_id_t elim_vertex,
//                        const mrf_graph_type& mrf,
//                        const size_t max_tree_width,
//                        const size_t max_factor_size,
//                        jtree_list& jt_list);







// /**
//  * Convert a jtree_list into a jtree_graph 
//  */
// void jtree_list_to_jtree_graph(const jtree_list& jt_list,
//                                const mrf_graph_type& mrf,
//                                const size_t num_factors,
//                                jtree_graph_type& jt_graph);



// /**
//  * Scan the junction tree list to ensure that all invariants hold.
//  */
// bool validate_jtree_list(const jtree_list& jt_list);


#endif
