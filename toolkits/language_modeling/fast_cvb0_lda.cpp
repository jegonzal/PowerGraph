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


#include "cvb0_lda_common.hpp"


#include <graphlab/macros_def.hpp>


class cvb0_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, factor_type>,
  public graphlab::IS_POD_TYPE {
public:

  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  } // end of gather_edges 

  factor_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const {
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
    factor_type& belief = edge.data().belief;
    const uint32_t count = edge.data().count;
    // Resample the topics
    double sum = 0, old_sum = 0;
    for(size_t t = 0; t < NTOPICS; ++t) {
      old_sum += belief[t];
      doc_topic_count[t] -= belief[t];
      word_topic_count[t] -= belief[t];
      GLOBAL_TOPIC_COUNT[t] -= belief[t];
      const double n_dt = 
        std::max(count_type(doc_topic_count[t]), count_type(0));
      ASSERT_GE(n_dt, 0);
      const double n_wt = 
        std::max(count_type(word_topic_count[t]), count_type(0)); 
      ASSERT_GE(n_wt, 0);
      const double n_t  = 
        std::max(count_type(GLOBAL_TOPIC_COUNT[t]), count_type(0)); 
      ASSERT_GE(n_t, 0);
      belief[t] = (ALPHA + n_dt) * (BETA + n_wt) / (BETA * NWORDS + n_t);
      sum += belief[t];
    } // End of loop over each token
    ASSERT_GT(sum, 0);
    if(old_sum == 0) {
      size_t asg = graphlab::random::multinomial(belief);
      for(size_t i = 0; i < NTOPICS; ++i) belief[i] = 0;
      belief[asg] = count;
      return belief;
    }
    for(size_t t = 0; t < NTOPICS; ++t) {
      belief[t] = count * (belief[t]/sum);
      doc_topic_count[t] += belief[t];
      word_topic_count[t] += belief[t];
      GLOBAL_TOPIC_COUNT[t] += belief[t];
    }
    return belief;
  } // end of gather


  void apply(icontext_type& context, vertex_type& vertex,
             const factor_type& sum) {
    const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
    ASSERT_GT(num_neighbors, 0);
    // There should be no new edge data since the vertex program has been cleared
    vertex_data& vdata = vertex.data();
    ASSERT_EQ(sum.size(), NTOPICS);
    ASSERT_EQ(vdata.factor.size(), NTOPICS);
    vdata.nupdates++; vdata.nchanges = 0; 
    vdata.factor = sum;
  } // end of apply

  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const { 
    return graphlab::ALL_EDGES; 
  }; // end of scatter edges

  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const {
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    context.signal(other_vertex);
  } // end of scatter function

}; // end of cvb0_lda_vertex_program





typedef graphlab::omni_engine<cvb0_lda_vertex_program> engine_type;
typedef cvb0_lda_vertex_program::icontext_type icontext_type;
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
    "The fast Collapsed Variational Bayes Alg for the LDA model implements\n" 
    "a highly asynchronous version of parallel CVB0 in which document\n"
    "and word counts are maintained in an eventually consistent\n"
    "manner.\n"
    "\n"
    "The standard usage is: \n"
    "\t./fast_cvb0_lda --dictionary dictionary.txt --matrix doc_word_count.tsv\n"
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

























