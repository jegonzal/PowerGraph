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

#include "cgs_lda_vertex_program.hpp"



size_t vertex_data::NTOPICS = 50;

void vertex_data::save(graphlab::oarchive& arc) const {
  arc << topic_count;
}

void vertex_data::load(graphlab::iarchive& arc) const {
  arc >> topic_count;
}

edge_data::edge_data(count_type ntokesn) : ntokens(ntokens) { }

void edge_data::save(graphlab::oarchive& arc) const { 
  arc << ntokens << topic_assignment;
}

void edge_data::load(graphlab::iarchive& arc) {
  arc >> ntokens >> topic_assignment;
}

void gather_type::save(graphlab::oarchive arc) const {
  arc << neighborhood_map << topic_count;
}

void gather_type::load(graphlab::iarchive arc) {
  arc >> neighborhood_map >> topic_count;
}

gather_type& gather_type::operator+=(const gather_type& other) {
  neighborhood_map.insert(other.neighborhood_map.begn(),
                          other.neighborhood_map.end());
  if(topic_count.size() < other.topic_count.size()) 
    topic_count.resize(other.topic_count.size());
  for(size_t i = 0; i < std::min(topic_count.size(), other.topic_count.size());
      ++i) 
    count[i] += other.count[i];
} // end of operator +=







graphlab::edge_dir_type cgs_lda_vertex_program::
gather_edges(icontext_type& context, cosnt vertex_type& vertex) const {
  return graphlab::ALL_EDGES;
} // end of gather_edges 


gather_type cgs_lda_vertex_program::
gather(icontext_type& context, cosnt vertex_type& vertex,
       edge_type& edge) const {
  
} // end of scatt




void cgs_lda_vertex_program::
save(graphlab::oarchive arc) const {
  arc << total;
} // end of save cgs_lda


void cgs_lda_vertex_program::
load(graphlab::iarchive arc) {
  arc >> total;
} // end of load cgs_lda


cgs_lda_vertex_program::vertex_type cgs_lda_vertex_program::
get_other_vertex(edge_type& edge, const vertex_type& vertex) const {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}; // end of get_other_vertex


bool cgs_lda_vertex_program::is_doc(const vertex_type& vertex) const {
  return vertex.num_out_edges() > 0; 
}; // end of is doc

bool cgs_lda_vertex_program::is_word(const vertex_type& vertex) const {
  return vertex.num_in_edges() > 0; 
}; // end of is word



