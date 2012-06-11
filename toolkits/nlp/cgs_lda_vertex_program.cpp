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




/**
 * \brief The gather type used to accumulate information about the
 * words in a document
 */
struct gather_type {   
  typedef std::pair<factor_type, assignment_type> edge_pair_type;
  typedef std::map<graphlab::vertex_id_type, edge_pair_type> 
  neighborhood_map_type;
  
  neighborhood_map_type neighborhood_map;
  factor_type topic_count;

  gather_type() { }
  gather_type(const factor_type& topic_count);
  gather_type(graphlab::vertex_id_type vid,
              const factor_type& factor,
              const assignment_type& assignment);
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);

  gather_type& operator+=(const gather_type& other);
}; // end of gather type


class cgs_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, gather_type> {
private:
  typedef std::map<graphlab::vertex_id_type, assignment_type> edge_data_map_type;
  edge_data_map_type new_edge_data;
public:
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);

  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const;
  gather_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const;
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& sum);
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const;
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const;
}; // end of cgs_lda_vertex_program



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



graphlab::edge_dir_type cgs_lda_vertex_program::
gather_edges(icontext_type& context, const vertex_type& vertex) const {
  return graphlab::ALL_EDGES;
} // end of gather_edges 


gather_type cgs_lda_vertex_program::
gather(icontext_type& context, const vertex_type& vertex,
       edge_type& edge) const {
  const vertex_type other_vertex = get_other_vertex(edge, vertex);
  if(is_doc(other_vertex)) {
    gather_type ret; ret.topic_count.resize(NTOPICS);
    const assignment_type& assignment = edge.data();
    foreach(topic_id_type asg, assignment) ++ret.topic_count[asg];
    return ret;
  } else {
    return gather_type(other_vertex.id(), other_vertex.data().factor,
                       edge.data());
  }
} // end of scatter



void cgs_lda_vertex_program::
apply(icontext_type& context, vertex_type& vertex, const gather_type& sum) {
  const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
  ASSERT_GT(num_neighbors, 0);
  ASSERT_EQ(new_edge_data.size(), 0); 
  ASSERT_EQ(vertex.data().factor.size(), NTOPICS);
  if(is_word(vertex)) {
    vertex.data().nupdates++; 
    vertex.data().factor = sum.topic_count;
  } else { ASSERT_TRUE(is_doc(vertex));
    vertex_data& vdata = vertex.data();
    vdata.nupdates++; vdata.nchanges = 0;
    factor_type& doc_topic_count = vdata.factor;
    // run the actual gibbs sampling 
    std::vector<double> prob(NTOPICS);
    typedef gather_type::neighborhood_map_type::value_type pair_type; 
    foreach(const pair_type& nbr_pair, sum.neighborhood_map) {
      const graphlab::vertex_id_type wordid = nbr_pair.first;
      factor_type word_topic_count = nbr_pair.second.first;
      assignment_type assignment = nbr_pair.second.second;
      ASSERT_EQ(word_topic_count.size(), NTOPICS);
      // Resample the topics
      foreach(topic_id_type& asg, assignment) {
        const topic_id_type old_asg = asg;
        if(asg != NULL_TOPIC) { // construct the cavity
          --doc_topic_count[asg];
          if(word_topic_count[asg] > 0) --word_topic_count[asg];
          --GLOBAL_TOPIC_COUNT[asg];
        }
        for(size_t t = 0; t < NTOPICS; ++t) {
          const double n_dt = doc_topic_count[t]; ASSERT_GE(n_dt, 0);
          const double n_wt = word_topic_count[t]; ASSERT_GE(n_wt, 0);
          const double n_t  = GLOBAL_TOPIC_COUNT[t]; ASSERT_GE(n_t, 0);
          prob[t] = (ALPHA + n_dt) * (BETA + n_wt) / (BETA * NWORDS + n_t);
        }
        asg = graphlab::random::multinomial(prob);
        ++doc_topic_count[asg];
        ++word_topic_count[asg];                    
        ++GLOBAL_TOPIC_COUNT[asg];
        // record a change if one occurs
        if(old_asg != asg) vdata.nchanges++;
      } // End of loop over each token
      // test to see if the topic assignments have change
      // sort the topic assignment to be in a "canonical order" 
      std::sort(assignment.begin(), assignment.end());
      const assignment_type& old_assignment = nbr_pair.second.second;
      bool is_same = (old_assignment.size() == assignment.size());  
      for(size_t i = 0; i < assignment.size() && is_same; ++i)
        is_same = (assignment[i] == old_assignment[i]);
      if(!is_same) new_edge_data[wordid] = assignment;
    } // end of loop over neighbors
  } // end of else document
} // end of apply



graphlab::edge_dir_type cgs_lda_vertex_program::
scatter_edges(icontext_type& context, const vertex_type& vertex) const { 
  return graphlab::ALL_EDGES; 
}; // end of scatter edges


void cgs_lda_vertex_program::
scatter(icontext_type& context, const vertex_type& vertex, 
        edge_type& edge) const {
  if(is_doc(vertex)) { 
    const vertex_type word_vertex = get_other_vertex(edge, vertex);
    ASSERT_TRUE(is_word(word_vertex));
    // if this is a document then update the topic assignment along the edge
    edge_data_map_type::const_iterator iter = 
      new_edge_data.find(word_vertex.id());
    // If there is an assignment then something changed so update and
    // signal
    if(iter != new_edge_data.end()) {
      const assignment_type& new_topic_assignment = iter->second;
      ASSERT_EQ(new_topic_assignment.size(), edge.data().size());
      edge.data() = new_topic_assignment;
    }
  } 
  context.signal(get_other_vertex(edge, vertex));
} // end of scatter function



void cgs_lda_vertex_program::save(graphlab::oarchive& arc) const {
  arc << new_edge_data;
} // end of save cgs_lda


void cgs_lda_vertex_program::load(graphlab::iarchive& arc) {
  arc >> new_edge_data;
} // end of load cgs_lda




void run_cgs_lda(graphlab::distributed_control& dc, 
                 graph_type& graph,
                 graphlab::command_line_options& clopts) {
  typedef graphlab::omni_engine<cgs_lda_vertex_program> engine_type;
  typedef cgs_lda_vertex_program::icontext_type icontext_type;
  typedef topk_aggregator<icontext_type> topk_type;
  engine_type engine(dc, graph, clopts, "asynchronous");
  ///! Add an aggregator
  bool success = false;
  success = engine.add_vertex_aggregator<topk_type>("topk",
                                                    topk_type::map, 
                                                    topk_type::finalize);
  ASSERT_TRUE(success);
  success = engine.aggregate_periodic("topk", INTERVAL);
  ASSERT_TRUE(success);
  ///! schedule only documents
  engine.map_reduce_vertices<graphlab::empty>(signal_docs<icontext_type>);
  std::cout << "Running The Collapsed Gibbs Sampler" << std::endl;
  graphlab::timer timer;
  engine.start();  
  const double runtime = timer.current_time();
  if(dc.procid() == 0) {
    std::cout << "----------------------------------------------------------"
              << std::endl;
    std::cout << "Final Runtime (seconds):   " << runtime 
              << std::endl
              << "Updates executed: " << engine.num_updates() << std::endl
              << "Update Rate (updates/second): " 
              << engine.num_updates() / runtime << std::endl;
  }
} // end of run_cgs_lda





