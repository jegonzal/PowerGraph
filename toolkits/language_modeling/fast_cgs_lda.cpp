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


#include <vector>
#include <algorithm>
#include <graphlab.hpp>


#include "cgs_lda_common.hpp"


#include <graphlab/macros_def.hpp>


// struct gather_type {
//   factor_type factor;
//   size_t nchanges;
//   gather_type() : nchanges(0) { };
//   void save(graphlab::oarchive& arc) const { arc << factor << nchanges; }
//   void load(graphlab::iarchive& arc) { arc >> factor >> nchanges; }
//   gather_type& operator+=(const gather_type& other) {
//     factor += other.factor;
//     nchanges += other.nchanges;
//     return *this;
//   }
// }; // end of gather type


class cgs_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, factor_type>,
  public graphlab::IS_POD_TYPE {
public:

  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  } // end of gather_edges 

  factor_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const {
    factor_type factor(NTOPICS);
    const assignment_type& assignment = edge.data();
    foreach(topic_id_type asg, assignment) {
      if(asg != NULL_TOPIC) ++factor[asg];
    }
    return factor;
  } // end of gather

  void apply(icontext_type& context, vertex_type& vertex,
             const factor_type& sum_factor) {
    const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
    ASSERT_GT(num_neighbors, 0);
    // There should be no new edge data since the vertex program has been cleared
    vertex_data& vdata = vertex.data();
    ASSERT_EQ(sum_factor.size(), NTOPICS);
    ASSERT_EQ(vdata.factor.size(), NTOPICS);
    vdata.nupdates++; 
    vdata.factor = sum_factor;
  } // end of apply

  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {
    factor_type& doc_topic_count =  is_doc(edge.source()) ? 
      edge.source().data().factor : edge.target().data().factor;
    factor_type& word_topic_count = is_word(edge.source()) ? 
      edge.source().data().factor : edge.target().data().factor;
    ASSERT_EQ(doc_topic_count.size(), NTOPICS);
    ASSERT_EQ(word_topic_count.size(), NTOPICS);
    // run the actual gibbs sampling 
    std::vector<double> prob(NTOPICS);
    assignment_type& assignment = edge.data();
    // Resample the topics
    size_t nchanges = 0;
    foreach(topic_id_type& asg, assignment) {
      const topic_id_type old_asg = asg;
      if(asg != NULL_TOPIC) { // construct the cavity
        --doc_topic_count[asg];
        --word_topic_count[asg];
        --GLOBAL_TOPIC_COUNT[asg];
      } 
      for(size_t t = 0; t < NTOPICS; ++t) {
        const double n_dt = 
          std::max(count_type(doc_topic_count[t]), count_type(0));
        const double n_wt = 
          std::max(count_type(word_topic_count[t]), count_type(0)); 
        const double n_t  = 
          std::max(count_type(GLOBAL_TOPIC_COUNT[t]), count_type(0)); 
        prob[t] = (ALPHA + n_dt) * (BETA + n_wt) / (BETA * NWORDS + n_t);
      }
      asg = graphlab::random::multinomial(prob);
      // asg = std::max_element(prob.begin(), prob.end()) - prob.begin();
      ++doc_topic_count[asg];
      ++word_topic_count[asg];                    
      ++GLOBAL_TOPIC_COUNT[asg];
      if(asg != old_asg) ++nchanges;
    } // End of loop over each token
    // singla the other vertex
    context.signal(get_other_vertex(edge, vertex));
  } // end of scatter function

}; // end of cgs_lda_vertex_program





typedef graphlab::omni_engine<cgs_lda_vertex_program> engine_type;
typedef cgs_lda_vertex_program::icontext_type icontext_type;
typedef topk_aggregator<icontext_type> topk_type;
typedef selective_signal<icontext_type> signal_only;
typedef global_counts_aggregator<icontext_type> global_counts_agg;





int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::distributed_control dc;

  // Parse command line options -----------------------------------------------
  const std::string description = 
    "\n=========================================================================\n"
    "The fast Collapsed Gibbs Sampler for the LDA model implements\n" 
    "a highly asynchronous version of parallel LDA in which document\n"
    "and word counts are maintained in an eventually consistent\n"
    "manner.\n"
    "\n"
    "The standard usage is: \n"
    "\t./fast_cgs_lda --dictionary dictionary.txt --matrix doc_word_count.tsv\n"
    "where dictionary.txt contains: \n"
    "\taaa \n\taaai \n\tabalone \n\t   ... \n\n"
    "and doc_word_count.tsv is formatted <docid> <wordid> <count>:\n"
    "\t0\t0\t3\n"
    "\t0\t5\t1\n"
    "\t ...\n\n"
    "To learn more about the NLP package and its applications visit\n\n"
    "\t\t http://graphlab.org \n\n"
    "Additional Options";
  graphlab::command_line_options clopts(description);
  std::string matrix_dir; 
  std::string dictionary_fname;
  clopts.attach_option("dictionary", &dictionary_fname, dictionary_fname,
                       "The file containing the list of unique words");
  clopts.add_positional("dictionary");
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory or file containing the matrix data.");
  clopts.add_positional("matrix");
  clopts.attach_option("ntopics", &NTOPICS, NTOPICS,
                       "Number of topics to use.");
  clopts.attach_option("alpha", &ALPHA, ALPHA,
                       "The document hyper-prior");
  clopts.attach_option("beta", &BETA, BETA,                       
                       "The word hyper-prior");
  clopts.attach_option("topk", &TOPK, TOPK,
                       "The number of words to report");
  clopts.attach_option("interval", &INTERVAL, INTERVAL,
                       "statistics reporting interval");
  clopts.attach_option("max_count", &MAX_COUNT, MAX_COUNT,
                       "The maximum number of occurences of a word in a document.");

  if(!clopts.parse(argc, argv)) {
    graphlab::mpi_tools::finalize();
    return clopts.is_set("help")? EXIT_SUCCESS : EXIT_FAILURE;
  }

  if(dictionary_fname.empty()) {
    logstream(LOG_ERROR) << "No dictionary file was provided." << std::endl;
    return EXIT_FAILURE;
  }

  if(matrix_dir.empty()) {
    logstream(LOG_ERROR) << "No matrix file was provided." << std::endl;
    return EXIT_FAILURE;
  }


  ///! Initialize global variables
  GLOBAL_TOPIC_COUNT.resize(NTOPICS);
  bool success = load_dictionary(dictionary_fname); 
  if(!success) {
    logstream(LOG_ERROR) << "Error loading dictionary." << std::endl;
    return EXIT_FAILURE;
  }
  
  ///! load the graph
  graph_type graph(dc, clopts);  
  success = load_and_initialize_graph(dc, graph, matrix_dir);
  if(!success) {
    logstream(LOG_ERROR) << "Error loading graph." << std::endl;
    return EXIT_FAILURE;
  }


  const size_t ntokens = graph.map_reduce_edges<size_t>(count_tokens);
  dc.cout() << "Total tokens: " << ntokens << std::endl;


  
  engine_type engine(dc, graph, clopts, "asynchronous");
  ///! Add an aggregator
  success = 
    engine.add_vertex_aggregator<topk_type>
    ("topk", topk_type::map, topk_type::finalize) &&
    engine.aggregate_periodic("topk", INTERVAL);
  ASSERT_TRUE(success);
  success = 
    engine.add_vertex_aggregator<factor_type>
    ("global_counts", global_counts_agg::map, global_counts_agg::finalize) &&
    engine.aggregate_periodic("global_counts", 5);
  ASSERT_TRUE(success);


  ///! schedule only documents
  dc.cout() << "Running The Collapsed Gibbs Sampler" << std::endl;
  engine.map_reduce_vertices<graphlab::empty>(signal_only::docs);
  graphlab::timer timer;
  engine.start();  
  const double runtime = timer.current_time();
    dc.cout() 
    << "----------------------------------------------------------" << std::endl
    << "Final Runtime (seconds):   " << runtime 
    << std::endl
    << "Updates executed: " << engine.num_updates() << std::endl
    << "Update Rate (updates/second): " 
    << engine.num_updates() / runtime << std::endl;



  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;


} // end of main

























