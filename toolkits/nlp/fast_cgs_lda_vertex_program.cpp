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

#include "fast_cgs_lda_vertex_program.hpp"



#include <graphlab/macros_def.hpp>


class fast_cgs_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, factor_type>,
  public graphlab::IS_POD_TYPE {
public:
  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const;
  factor_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const;
  void apply(icontext_type& context, vertex_type& vertex,
             const factor_type& sum);
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const;
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const;
}; // end of fast_cgs_lda_vertex_program




graphlab::edge_dir_type fast_cgs_lda_vertex_program::
gather_edges(icontext_type& context, const vertex_type& vertex) const {
  return graphlab::ALL_EDGES;
} // end of gather_edges 


factor_type fast_cgs_lda_vertex_program::
gather(icontext_type& context, const vertex_type& vertex,
       edge_type& edge) const {
  factor_type topic_count; topic_count.resize(NTOPICS);
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
      --GLOBAL_TOPIC_COUNT[asg];
    } 
    for(size_t t = 0; t < NTOPICS; ++t) {
      const double n_dt = 
        std::max(count_type(doc_topic_count[t]), count_type(0));
      ASSERT_GE(n_dt, 0);
      const double n_wt = 
        std::max(count_type(word_topic_count[t]), count_type(0)); 
      ASSERT_GE(n_wt, 0);
      const double n_t  = GLOBAL_TOPIC_COUNT[t]; ASSERT_GE(n_t, 0);
      prob[t] = (ALPHA + n_dt) * (BETA + n_wt) / (BETA * NWORDS + n_t);
    }
    asg = graphlab::random::multinomial(prob);
    ++doc_topic_count[asg];
    ++word_topic_count[asg];                    
    ++GLOBAL_TOPIC_COUNT[asg];
    ++topic_count[asg];
  } // End of loop over each token
  return topic_count;
} // end of gather



void fast_cgs_lda_vertex_program::
apply(icontext_type& context, vertex_type& vertex, const factor_type& sum) {
  const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
  ASSERT_GT(num_neighbors, 0);
  // There should be no new edge data since the vertex program has been cleared
  vertex_data& vdata = vertex.data();
  ASSERT_EQ(sum.size(), NTOPICS);
  ASSERT_EQ(vdata.factor.size(), NTOPICS);
  vdata.nupdates++; vdata.nchanges = 0; 
  for(size_t t = 0; t < vdata.factor.size(); ++t) {
    vdata.nchanges += std::abs(vdata.factor[t] - sum[t]);
    vdata.factor[t] = sum[t];
  }
} // end of apply



graphlab::edge_dir_type fast_cgs_lda_vertex_program::
scatter_edges(icontext_type& context, const vertex_type& vertex) const { 
  return graphlab::ALL_EDGES; 
}; // end of scatter edges


void fast_cgs_lda_vertex_program::
scatter(icontext_type& context, const vertex_type& vertex, 
        edge_type& edge) const {
  const vertex_type other_vertex = get_other_vertex(edge, vertex);
  context.signal(other_vertex);
} // end of scatter function



void run_fast_cgs_lda(graphlab::distributed_control& dc, 
                      graph_type& graph,
                      graphlab::command_line_options& clopts) {
  typedef graphlab::omni_engine<fast_cgs_lda_vertex_program> engine_type;
  typedef fast_cgs_lda_vertex_program::icontext_type icontext_type;
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
} // end of run_fast_cgs_lda






