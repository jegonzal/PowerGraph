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



#include <graphlab/macros_def.hpp>

const topic_id_type NULL_TOPIC(-1);


void vertex_data::save(graphlab::oarchive& arc) const {
  arc << factor << nupdates;
}


void vertex_data::load(graphlab::iarchive& arc) {
  arc >> factor >> nupdates;
}


gather_type::gather_type(const factor_type& topic_count) :
  topic_count(topic_count) { } // end of gather_type constructor

gather_type::gather_type(graphlab::vertex_id_type vid,
                         const factor_type& factor,
                         const assignment_type& assignment) {
  neighborhood_map[vid] = edge_pair_type(factor, assignment);
}


void gather_type::save(graphlab::oarchive& arc) const {
  arc << neighborhood_map << topic_count;
}

void gather_type::load(graphlab::iarchive& arc) {
  arc >> neighborhood_map >> topic_count;
}

gather_type& gather_type::operator+=(const gather_type& other) {
  neighborhood_map.insert(other.neighborhood_map.begin(),
                          other.neighborhood_map.end());
  if(topic_count.size() < other.topic_count.size()) 
    topic_count.resize(other.topic_count.size());
  for(size_t i = 0; i < std::min(topic_count.size(), other.topic_count.size());
      ++i) 
    topic_count[i] += other.topic_count[i];
  return *this;
} // end of operator +=




double cgs_lda_vertex_program::ALPHA   = 0.1;
double cgs_lda_vertex_program::BETA    = 0.1;
size_t cgs_lda_vertex_program::NITERS  = -1;
size_t cgs_lda_vertex_program::NTOPICS = 50;
size_t cgs_lda_vertex_program::NWORDS  = 0;
factor_type cgs_lda_vertex_program::global_topic_count;

void cgs_lda_vertex_program::set_ntopics(size_t ntopics) {
  NTOPICS = ntopics;
  global_topic_count.resize(ntopics);
}

graphlab::edge_dir_type cgs_lda_vertex_program::
gather_edges(icontext_type& context, const vertex_type& vertex) const {
  return graphlab::ALL_EDGES;
} // end of gather_edges 


gather_type cgs_lda_vertex_program::
gather(icontext_type& context, const vertex_type& vertex,
       edge_type& edge) const {
  gather_type ret_val; ret_val.topic_count.resize(NTOPICS);
  vertex_type other_vertex = get_other_vertex(edge, vertex);
  // VIOLATING THE ABSTRACTION!
  vertex_data& vdata = graph_type::vertex_type(vertex).data();
  // VIOLATING THE ABSTRACTION!
  vertex_data& other_vdata = other_vertex.data();
  factor_type& doc_topic_count = 
    is_doc(vertex) ? vdata.factor : other_vdata.factor;
  factor_type& word_topic_count = 
    is_word(vertex) ? vdata.factor : other_vdata.factor;
  ASSERT_EQ(doc_topic_count.size(), NTOPICS);
  ASSERT_EQ(word_topic_count.size(), NTOPICS);
  // run the actual gibbs sampling 
  std::vector<double> prob(NTOPICS);
  assignment_type& assignment = edge.data();
  // Resample the topics
  foreach(topic_id_type& asg, assignment) {
    if(asg != NULL_TOPIC) { // construct the cavity
      --doc_topic_count[asg];
      --word_topic_count[asg];
      --global_topic_count[asg];
    } 
    for(size_t t = 0; t < NTOPICS; ++t) {
      const double n_dt = 
        std::max(count_type(doc_topic_count[t]), count_type(0));
      ASSERT_GE(n_dt, 0);
      const double n_wt = 
        std::max(count_type(word_topic_count[t]), count_type(0)); 
      ASSERT_GE(n_wt, 0);
      const double n_t  = global_topic_count[t]; ASSERT_GE(n_t, 0);
      prob[t] = (ALPHA + n_dt) * (BETA + n_wt) / (BETA * NWORDS + n_t);
    }
    asg = graphlab::random::multinomial(prob);
    ++doc_topic_count[asg];
    ++word_topic_count[asg];                    
    ++global_topic_count[asg];
    ++ret_val.topic_count[asg];
  } // End of loop over each token
  return ret_val;
} // end of scatter



void cgs_lda_vertex_program::
apply(icontext_type& context, vertex_type& vertex, const gather_type& sum) {
  const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
  ASSERT_GT(num_neighbors, 0);
  // There should be no new edge data since the vertex program has been cleared
  vertex_data& vdata = vertex.data();
  ASSERT_EQ(sum.topic_count.size(), NTOPICS);
  ASSERT_EQ(vdata.factor.size(), NTOPICS);
  vdata.nupdates++; vdata.nchanges = 0; 
  for(size_t t = 0; t < vdata.factor.size(); ++t) {
    vdata.nchanges += std::abs(vdata.factor[t] - sum.topic_count[t]);
    vdata.factor[t] = sum.topic_count[t];
  }
} // end of apply



graphlab::edge_dir_type cgs_lda_vertex_program::
scatter_edges(icontext_type& context, const vertex_type& vertex) const { 
  return graphlab::ALL_EDGES; 
}; // end of scatter edges


void cgs_lda_vertex_program::
scatter(icontext_type& context, const vertex_type& vertex, 
        edge_type& edge) const {
  const vertex_type other_vertex = get_other_vertex(edge, vertex);
  context.signal(other_vertex);
} // end of scatter function



void cgs_lda_vertex_program::save(graphlab::oarchive& arc) const {
  arc << new_edge_data;
} // end of save cgs_lda


void cgs_lda_vertex_program::load(graphlab::iarchive& arc) {
  arc >> new_edge_data;
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






