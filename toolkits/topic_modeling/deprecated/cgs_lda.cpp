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

#include "cgs_lda_common.hpp"

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

  gather_type() { }

  gather_type(const graph_type::edge_type& edge,
              const graph_type::vertex_type& vertex) {
    const graph_type::vertex_type other_vertex =
      get_other_vertex(edge, vertex);
    neighborhood_map[other_vertex.id()] =
      edge_pair_type(other_vertex.data().factor, edge.data().assignment);
  }

  void save(graphlab::oarchive& arc) const { arc << neighborhood_map; }

  void load(graphlab::iarchive& arc) { arc >> neighborhood_map; }

  gather_type& operator+=(const gather_type& other) {
    neighborhood_map.insert(other.neighborhood_map.begin(),
                            other.neighborhood_map.end());
    return *this;
  } // end of operator +=

}; // end of gather type




class cgs_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, gather_type> {
private:
  typedef std::map<graphlab::vertex_id_type, assignment_type> edge_data_map_type;
  edge_data_map_type new_edge_data;
public:
  void save(graphlab::oarchive& arc) const {
    arc << new_edge_data;
  } // end of save cgs_lda

  void load(graphlab::iarchive& arc) {
    arc >> new_edge_data;
  } // end of load cgs_lda


  edge_dir_type gather_edges(icontext_type& context,
                             const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  } // end of gather_edges

  gather_type gather(icontext_type& context, const vertex_type& vertex,
                     edge_type& edge) const {
    return gather_type(edge, vertex);
  } // end of gather


  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& sum) {
    const size_t num_neighbors = vertex.num_in_edges() + vertex.num_out_edges();
    ASSERT_GT(num_neighbors, 0);
    ASSERT_EQ(new_edge_data.size(), 0);
    vertex_data& vdata = vertex.data();
    factor_type& factor = vdata.factor;
    ASSERT_EQ(factor.size(), NTOPICS);

    // first update the factor count for this vertex
    typedef gather_type::neighborhood_map_type::value_type pair_type;
    foreach(const pair_type& nbr_pair, sum.neighborhood_map) {
      const assignment_type& assignment = nbr_pair.second.second;
      foreach(const topic_id_type& asg, assignment) {
        if(asg != NULL_TOPIC) ++factor[asg];
      } // end of loop over assignments
    } // end of loop over neighborhood

    // Resample the vertex
    vdata.nchanges = 0;
    // run the actual gibbs sampling
    std::vector<double> prob(NTOPICS);
    typedef gather_type::neighborhood_map_type::value_type pair_type;
    foreach(const pair_type& nbr_pair, sum.neighborhood_map) {
      const graphlab::vertex_id_type other_id = nbr_pair.first;
      factor_type other_factor = nbr_pair.second.first;
      assignment_type assignment = nbr_pair.second.second;
      factor_type& doc_topic_count  = is_doc(vertex)? factor : other_factor;
      factor_type& word_topic_count = is_word(vertex)? factor : other_factor;
      ASSERT_EQ(word_topic_count.size(), NTOPICS);
      ASSERT_EQ(doc_topic_count.size(), NTOPICS);
      // Resample the topics
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
          ASSERT_GE(n_dt, 0);
          const double n_wt =
            std::max(count_type(word_topic_count[t]), count_type(0));
          ASSERT_GE(n_wt, 0);
          const double n_t  =
            std::max(count_type(GLOBAL_TOPIC_COUNT[t]), count_type(0));
          ASSERT_GE(n_t, 0);
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
      if(!is_same) new_edge_data[other_id] = assignment;
    } // end of loop over neighbors
  } // end of apply


  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
    return graphlab::ALL_EDGES;
  }; // end of scatter edges

  void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const  {
    const vertex_type other_vertex = get_other_vertex(edge, vertex);
    edge_data_map_type::const_iterator iter =
      new_edge_data.find(other_vertex.id());
    // If there is an assignment then something changed
    if(iter != new_edge_data.end()) {
      const assignment_type& new_topic_assignment = iter->second;
      ASSERT_EQ(new_topic_assignment.size(),
                edge.data().assignment.size());
      edge.data().assignment = new_topic_assignment;
    }
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
    "\taaa \n\taaai \n\tabalone \n\t   ... \n"
    "each line number corresponds to wordid (i.e aaa has wordid=0)\n\n"
    "and doc_word_count.tsv is formatted <docid> <wordid> <count>:\n"
    "(where wordid is indexed starting from zero and docid are positive integers)\n"
    "\t0\t0\t3\n"
    "\t0\t5\t1\n"
    "\t ...\n\n"
    "For JSON format, make sure docid are negative integers index starting from -2 \n\n"
    "To learn more about the NLP package and its applications visit\n\n"
    "\t\t http://graphlab.org \n\n"
    "Additional Options";
  graphlab::command_line_options clopts(description);
  std::string matrix_dir;
  std::string dictionary_fname;
  bool loadjson = false;
  clopts.attach_option("dictionary", &dictionary_fname, dictionary_fname,
                       "The file containing the list of unique words");
  clopts.add_positional("dictionary");
  clopts.attach_option("matrix", &matrix_dir, matrix_dir,
                       "The directory containing the matrix file");
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
  clopts.attach_option("loadjson",&loadjson,loadjson,
                        "Boolean if in json format (matrix arg is dir or gzip file)");

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
  success = load_and_initialize_graph(dc, graph, matrix_dir,loadjson);
  if(!success) {
    logstream(LOG_ERROR) << "Error loading graph." << std::endl;
    return EXIT_FAILURE;
  }



  engine_type engine(dc, graph, clopts, "synchronous");
  ///! Add an aggregator
  success =
    engine.add_vertex_aggregator<topk_type>
    ("topk", topk_type::map, topk_type::finalize) &&
    engine.aggregate_periodic("topk", INTERVAL);
  ASSERT_TRUE(success);
  // success =
  //   engine.add_vertex_aggregator<factor_type>
  //   ("global_counts", global_counts_agg::map, global_counts_agg::finalize) &&
  //   engine.aggregate_periodic("global_counts", 5);
  // ASSERT_TRUE(success);


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

























