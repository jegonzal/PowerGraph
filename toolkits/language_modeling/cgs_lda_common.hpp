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

#ifndef CGS_LDA_HPP
#define CGS_LDA_HPP
 
#include <vector>
#include <algorithm>
#include <graphlab/parallel/atomic.hpp>

// Global Types
// ============================================================================

typedef int count_type;

typedef uint16_t topic_id_type;

#define NULL_TOPIC (topic_id_type(-1))

/**
 * \brief The factor type is used to store the counts of tokens in
 * each topic for words, documents, and assignments.
 *
 * Atomic counts are used because we violate the abstraction by
 * modifying adjacent vertex data on scatter.  As a consequence
 * multiple threads on the same machine may try to update the same
 * vertex data at the same time.  The graphlab::atomic type ensures
 * that multiple increments are serially consistent.
 */
typedef std::vector< graphlab::atomic<count_type> > factor_type;


/**
 * \brief We use the factor type in accumulators and so we define an
 * operator+=
 */
inline factor_type& operator+=(factor_type& lvalue, 
                               const factor_type& rvalue) {
  if(!rvalue.empty()) {
    if(lvalue.empty()) lvalue = rvalue;
    else {
      for(size_t t = 0; t < lvalue.size(); ++t) lvalue[t] += rvalue[t];
    }
  }
  return lvalue;
} // end of operator +=





#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>


/**
 * \brief The assignment type is used on each edge to store the
 * assignments of each token.  There can be several occurrences of the
 * same word in a given document and so a vector is used to store the
 * assignments of each occurrence.
 */
typedef std::vector< topic_id_type > assignment_type;


// Global Variables
// ============================================================================

/**
 * \brief The alpha parameter determines the sparsity of topics for
 * each document.
 */
double ALPHA = 0.1;

/**
 * \brief the Beta parameter determines the sparsity of words in each
 * document.
 */
double BETA = 0.1;

/**
 * \brief the total number of topics to uses
 */
size_t NTOPICS = 50;

/**
 * \brief The total number of words in the dataset.
 */
size_t NWORDS = 0;

/**
 * \brief The total number of docs in the dataset.
 */
size_t NDOCS = 0;


/**
 * \brief The number of top words to display during execution (from
 * each topic).
 */
size_t TOPK = 5;

/**
 * \brief The interval to display topics during execution.
 */
size_t INTERVAL = 10;

/**
 * \brief The global variable storing the global topic count across
 * all machines.  This is maintained periodically using aggregation.
 */
factor_type GLOBAL_TOPIC_COUNT;

/**
 * \brief A dictionary of words used to print the top words during
 * execution.
 */
std::vector<std::string> DICTIONARY;

/**
 * \brief The maximum occurences allowed for an individual term-doc
 * pair. (edge data)
 */
size_t MAX_COUNT = 100;






// Graph Types
// ============================================================================

/**
 * \brief The vertex data represents each term and document in the
 * corpus and contains the counts of tokens in each topic.
 */
struct vertex_data {
  ///! The total number of updates
  uint32_t nupdates;
  ///! The total number of changes to adjacent tokens
  uint32_t nchanges;
  ///! The count of tokens in each topic
  factor_type factor;
  vertex_data() : nupdates(0), nchanges(0), factor(NTOPICS) { }
  void save(graphlab::oarchive& arc) const {
    arc << nupdates << nchanges << factor; 
  }
  void load(graphlab::iarchive& arc) { 
    arc >> nupdates >> nchanges >> factor; 
  }
}; // end of vertex_data


/**
 * \brief The edge data represents the individual tokens (word,doc)
 * pairs and their assignment to topics. 
 */
struct edge_data {
  ///! The number of changes on the last update
  uint16_t nchanges;
  ///! The assignment of all tokens
  assignment_type assignment;
  edge_data(size_t ntokens = 0) : nchanges(0), assignment(ntokens, NULL_TOPIC) { }
  void save(graphlab::oarchive& arc) const { arc << nchanges << assignment; }
  void load(graphlab::iarchive& arc) { arc >> nchanges >> assignment; }
}; // end of edge_daga


/**
 * \brief The LDA graph is a bipartite graph with docs connected to
 * terms if the term occurs in the document.  
 *
 * The edges store the number of occurrences of the term in the
 * document as a vector of the assignments of that term in that
 * document to topics.
 *
 * The vertices store the total topic counts.
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;


/**
 * \brief The graph loader is used by graph.load to parse lines of the
 * text data file.
 *
 * The global variable MAX_COUNT limits the number of tokens that can
 * be constructed on a particular edge.  
 */
bool graph_loader(graph_type& graph, const std::string& fname, 
                  const std::string& line) {
  ASSERT_FALSE(line.empty()); 
  const int BASE = 10;
  char* next_char_ptr = NULL;
  graph_type::vertex_id_type doc_id = 
    strtoul(line.c_str(), &next_char_ptr, BASE);
  if(next_char_ptr == NULL)  return false;
  const graph_type::vertex_id_type word_id = 
    strtoul(next_char_ptr, &next_char_ptr, BASE);
  if(next_char_ptr == NULL) return false;
  size_t count = 
    strtoul(next_char_ptr, &next_char_ptr, BASE);
  if(next_char_ptr == NULL) return false;
  // Threshold the count
  count = std::min(count, MAX_COUNT);
  // since this is a bipartite graph I need a method to number the
  // left and right vertices differently.  To accomplish I make sure
  // all vertices have non-zero ids and then negate the right vertex.
  // Unfortunatley graphlab reserves -1 and so we add 2 and negate.
  doc_id += 2; 
  ASSERT_GT(doc_id, 1); 
  doc_id = -doc_id;
  ASSERT_NE(doc_id, word_id);
  // Create an edge and add it to the graph
  graph.add_edge(doc_id, word_id, edge_data(count));
  return true; // successful load
}; // end of graph loader





/**
 * \brief Determine if the given vertex is a word vertex or a doc
 * vertex.
 *
 * For simplicity we connect docs --> words and therefore if a vertex
 * has in edges then it is a word.
 */
inline bool is_word(const graph_type::vertex_type& vertex) {
  return vertex.num_in_edges() > 0 ? 1 : 0;
}


/**
 * \brief Determine if the given vertex is a doc vertex
 *
 * For simplicity we connect docs --> words and therefore if a vertex
 * has out edges then it is a doc
 */
inline bool is_doc(const graph_type::vertex_type& vertex) {
  return vertex.num_out_edges() > 0 ? 1 : 0;
}

/**
 * \brief return the number of tokens on a particular edge.
 */
inline size_t count_tokens(const graph_type::edge_type& edge) {
  return edge.data().assignment.size();
}


/**
 * \brief Get the other vertex in the edge.
 */
inline graph_type::vertex_type 
get_other_vertex(const graph_type::edge_type& edge, 
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}



/**
 * \brief The topk aggregator is used to periodically compute and
 * display the topk most common words in each topic.
 *
 * The number of words is determined by the global variable \ref TOPK
 * and the interval is determined by the global variable \ref INTERVAL.
 *
 */
template<typename IContext>
class topk_aggregator {
  typedef IContext icontext_type;
  typedef std::pair<float, graphlab::vertex_id_type> cw_pair_type;
private:
  std::vector< std::set<cw_pair_type> > top_words;
  size_t nchanges, nupdates;
public:
  topk_aggregator(size_t nchanges = 0, size_t nupdates = 0) : 
    nchanges(nchanges), nupdates(nupdates) { }

  void save(graphlab::oarchive& arc) const { arc << top_words << nchanges; }
  void load(graphlab::iarchive& arc) { arc >> top_words >> nchanges; }


  topk_aggregator& operator+=(const topk_aggregator& other) {
    nchanges += other.nchanges;
    nupdates += other.nupdates;
    if(other.top_words.empty()) return *this;
    if(top_words.empty()) top_words.resize(NTOPICS);
    for(size_t i = 0; i < top_words.size(); ++i) {
      // Merge the topk
      top_words[i].insert(other.top_words[i].begin(), 
                          other.top_words[i].end());
      // Remove excess elements        
      while(top_words[i].size() > TOPK) 
        top_words[i].erase(top_words[i].begin());
    }
    return *this;
  } // end of operator +=
  
  static topk_aggregator map(icontext_type& context, 
                             const graph_type::vertex_type& vertex) {
    topk_aggregator ret_value;
    const vertex_data& vdata = vertex.data();
    ret_value.nchanges = vdata.nchanges;
    ret_value.nupdates = vdata.nupdates;
    if(is_word(vertex)) {
      const graphlab::vertex_id_type wordid = vertex.id();    
      ret_value.top_words.resize(vdata.factor.size());
      for(size_t i = 0; i < vdata.factor.size(); ++i) {
        const cw_pair_type pair(vdata.factor[i], wordid);
        ret_value.top_words[i].insert(pair);
      }
    } 
    return ret_value;   
  } // end of map function

  static void finalize(icontext_type& context,
                       const topk_aggregator& total) {
    if(context.procid() != 0) return;
    for(size_t i = 0; i < total.top_words.size(); ++i) {
      std::cout << "Topic " << i << ": ";
      rev_foreach(cw_pair_type pair, total.top_words[i])  {
        ASSERT_LT(pair.second, DICTIONARY.size());
         std::cout << DICTIONARY[pair.second] 
                   << "(" << pair.first << ")" << ", "; 
        // std::cout << DICTIONARY[pair.second] << ",  ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nNumber of token changes: " << total.nchanges << std::endl;
    std::cout << "\nNumber of updates:       " << total.nupdates << std::endl;
  } // end of finalize
}; // end of topk_aggregator struct



/**
 * \brief The global counts aggregator computes the total number of
 * tokens in each topic across all words and documents and then
 * updates the \ref GLOBAL_TOPIC_COUNT variable.
 */
template<typename IContext>
struct global_counts_aggregator {
  typedef graph_type::vertex_type vertex_type;
  static factor_type map(IContext& context, const vertex_type& vertex) {
    return vertex.data().factor;
  } // end of map function

  static void finalize(IContext& context, const factor_type& total) {
    size_t sum = 0;
    for(size_t t = 0; t < total.size(); ++t) {
      GLOBAL_TOPIC_COUNT[t] =
        std::max(count_type(total[t]/2), count_type(0));
      sum += GLOBAL_TOPIC_COUNT[t];
    }
    context.cout() << "Total Tokens: " << sum << std::endl;
  } // end of finalize
}; // end of global_counts_aggregator struct



/**
 * \brief The Likelihood aggregators maintains the current estimate of
 * the log-likelihood of the current token assignments.
 * 
 *  llik_words_given_topics = ...
 *    ntopics * (gammaln(nwords * beta) - nwords * gammaln(beta)) - ...
 *    sum_t(gammaln( n_t + nwords * beta)) +
 *    sum_w(sum_t(gammaln(n_wt + beta)));
 *
 *  llik_topics = ...
 *    ndocs * (gammaln(ntopics * alpha) - ntopics * gammaln(alpha)) + ...
 *    sum_d(sum_t(gammaln(n_td + alpha)) - gammaln(sum_t(n_td) + ntopics * alpha));
 */
template<typename IContext>
class likelihood_aggregator : public graphlab::IS_POD_TYPE {
  typedef graph_type::vertex_type vertex_type;
  double lik_words_given_topics;
  double lik_topics;
public:
  likelihood_aggregator() : lik_words_given_topics(0), lik_topics(0) { }
  
  likelihood_aggregator& operator+=(const likelihood_aggregator& other) {
    lik_words_given_topics += other.lik_words_given_topics;
    lik_topics += other.lik_topics;
    return *this;
  } // end of operator +=

  static likelihood_aggregator 
  map(IContext& context, const vertex_type& vertex) {
    using boost::math::lgamma;
    const factor_type& factor = vertex.data().factor;
    ASSERT_EQ(factor.size(), NTOPICS);
   likelihood_aggregator ret;
    if(is_word(vertex)) {
      for(size_t t = 0; t < NTOPICS; ++t) 
        ret.lik_words_given_topics += lgamma(factor[t] + BETA);
    } else {  ASSERT_TRUE(is_doc(vertex));
      size_t ntokens_in_doc = 0;
      for(size_t t = 0; t < NTOPICS; ++t) {
        ret.lik_topics += lgamma(factor[t] + ALPHA);
        ntokens_in_doc += factor[t];
      }
      ret.lik_topics -= lgamma(ntokens_in_doc + NTOPICS * ALPHA);
    }
    return ret;
  } // end of map function

  static void finalize(IContext& context, const likelihood_aggregator& total) {
    using boost::math::lgamma;
    // Address the global sum terms
    double denominator = 0;
    for(size_t t = 0; t < NTOPICS; ++t) {
      denominator += lgamma(GLOBAL_TOPIC_COUNT[t] + NWORDS * BETA);
    } // end of for loop

    const double lik_words_given_topics = 
      NTOPICS * (lgamma(NWORDS * BETA) - NWORDS * lgamma(BETA)) -
      denominator + total.lik_words_given_topics;

    const double lik_topics =
      NDOCS * (lgamma(NTOPICS * ALPHA) - NTOPICS * lgamma(ALPHA)) +
      total.lik_topics;
    
    const double lik = lik_words_given_topics + lik_topics;
    context.cout() << "Likelihood: " << lik << std::endl;
  } // end of finalize
}; // end of likelihood_aggregator struct




template<typename IContext>
struct selective_signal {
 static graphlab::empty 
 docs(IContext& context, const graph_type::vertex_type& vertex) {
   if(is_doc(vertex)) context.signal(vertex);
   return graphlab::empty();
 } // end of signal_docs
 static graphlab::empty 
 words(IContext& context, const graph_type::vertex_type& vertex) {
    if(is_word(vertex)) context.signal(vertex);
    return graphlab::empty();
  } // end of signal_words
}; // end of selective_signal




bool load_and_initialize_graph(graphlab::distributed_control& dc,
                               graph_type& graph,
                               const std::string& matrix_dir) {  
  dc.cout() << "Loading graph." << std::endl;
  graphlab::timer timer; timer.start();
  graph.load(matrix_dir, graph_loader);

  dc.cout() << ": Loading graph. Finished in " 
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Finalizing graph." << std::endl;
  timer.start();
  graph.finalize();
  dc.cout() << "Finalizing graph. Finished in " 
            << timer.current_time() << " seconds." << std::endl;

  dc.cout() << "Computing number of words and documents." << std::endl;
  NWORDS = graph.map_reduce_vertices<size_t>(is_word);
  NDOCS = graph.map_reduce_vertices<size_t>(is_doc);
  dc.cout() << "Number of words:   " << NWORDS;
  dc.cout() << "Number of docs:    " << NDOCS;
  //ASSERT_LT(NWORDS, DICTIONARY.size());
  return true;
} // end of load and initialize graph




/** populate the global dictionary */
bool load_dictionary(const std::string& fname)  {
  // std::cout << "staring load on: " 
  //           << graphlab::get_local_ip_as_str() << std::endl;
  const bool gzip = boost::ends_with(fname, ".gz");
  // test to see if the graph_dir is an hadoop path
  if(boost::starts_with(fname, "hdfs://")) {
    graphlab::hdfs hdfs;
    graphlab::hdfs::fstream in_file(hdfs, fname);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    fin.set_auto_close(false);
    if(gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good()) {
      logstream(LOG_ERROR) << "Error loading dictionary: "
                           << fname << std::endl;
      return false;
    }
    std::string term;
    while(std::getline(fin,term).good()) DICTIONARY.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } else {
    std::cout << "opening: " << fname << std::endl;
    std::ifstream in_file(fname.c_str(), 
                          std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_stream<boost::iostreams::input> fin;  
    if (gzip) fin.push(boost::iostreams::gzip_decompressor());
    fin.push(in_file);
    if(!fin.good() || !fin.good()) {
      logstream(LOG_ERROR) << "Error loading dictionary: "
                           << fname << std::endl;
      return false;
    }
    std::string term;
    std::cout << "Loooping" << std::endl;
    while(std::getline(fin, term).good()) DICTIONARY.push_back(term);
    if (gzip) fin.pop();
    fin.pop();
    in_file.close();
  } // end of else
  // std::cout << "Finished load on: " 
  //           << graphlab::get_local_ip_as_str() << std::endl;
  std::cout << "Dictionary Size: " << DICTIONARY.size() << std::endl;
  return true;
} // end of load dictionary




#include <graphlab/macros_undef.hpp>
#endif
