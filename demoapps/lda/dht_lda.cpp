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


/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */

#include <iostream>
#include <vector>

#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/delta_dht.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <graphlab.hpp>





#include "corpus.hpp"





#include <graphlab/macros_def.hpp>




struct topic_vector : public std::vector<topic_id_type> {
  typedef std::vector<topic_id_type> base; 
  topic_vector& operator+=(const topic_vector& other) {
    const size_t max_size = std::max(base::size(), other.size());
    base::resize(max_size);
    for(size_t i = 0; i < max_size; ++i) 
      (*this)[i] += (i < other.size()? other[i] : 0);
    return *this;
  }
  topic_vector& operator-=(const topic_vector& other) {
    const size_t max_size = std::max(size(), other.size());
    base::resize(max_size);
    for(size_t i = 0; i < max_size; ++i) 
      (*this)[i] -= (i < other.size()? other[i] : 0);
    return *this;
  } 
  topic_vector operator+(const topic_vector& other) const {
    topic_vector result(*this); result += other;
    return result;
  } 
  topic_vector operator-(const topic_vector& other) const {
    topic_vector result(*this); result -= other;
    return result;
  } 
  int l1norm() const { 
    int sum = 0;
    for(size_t i = 0; i < base::size(); ++i)
      sum += std::abs((*this)[i]);
    return sum;
  }
  int l1diff(const topic_vector& other) const {
    const size_t max_size = std::max(size(), other.size());
    int sum = 0;
    for(size_t i = 0; i < max_size; ++i) {
      sum += std::abs
        ((i < other.size()? other[i]:0) - (i < base::size()? (*this)[i]:0));
    }
    return sum;
  }
  void load(graphlab::iarchive& arc) { 
    arc >> *dynamic_cast<base*>(this); 
  }
  void save(graphlab::oarchive& arc) const { 
    arc << *dynamic_cast<const base*>(this); 
  }
};




void run_gibbs(const size_t niters, 
               const size_t ntopics,
               const double alpha,
               const double beta,
               const corpus& corpus,
               std::vector<topic_id_type>& topics,
               boost::unordered_map<doc_id_type, topic_vector>& n_dt,
               graphlab::delta_dht<word_id_type, topic_vector>& n_wt,
               graphlab::delta_dht<char, topic_vector>& n_t) {

  // Preallocate the vector to store the conditional
  std::vector<double> conditional(ntopics, 0);
  size_t nchanges = 0;

  for(size_t i = 0; i < topics.size(); ++i) {
    // Get the word and document for the ith token
    const word_id_type w = corpus.tokens[i].word;
    const doc_id_type  d = corpus.tokens[i].doc;
    const topic_id_type old_topic = topics[i];

    topic_vector& vec_n_dt(n_dt[d]);
    topic_vector& vec_n_wt(n_wt[w]);
    topic_vector& vec_n_t(n_t[0]);

    // resize the corresponding topic vectors if they are not currently allocated
    if(vec_n_dt.empty()) vec_n_dt.resize(ntopics, 0);
    if(vec_n_wt.empty()) vec_n_wt.resize(ntopics, 0);

    

    // Remove the word from the current counters
    if(old_topic != NULL_TOPIC) {
      --vec_n_dt[old_topic]; 
      --vec_n_wt[old_topic];
      --vec_n_t[old_topic];
    }

    // Construct the conditional
    double normalizer = 0;
    for(size_t t = 0; t < ntopics; ++t) {
      conditional[t] = (alpha + vec_n_dt[t]) * (beta + vec_n_wt[t]) /
        (beta * corpus.nwords + vec_n_t[t]);        
      normalizer += conditional[t];
    }
    ASSERT_GT(normalizer, 0);

    // Draw a new value
    topic_id_type new_topic = 0;
    // normalize and then sample
    for(size_t t = 0; t < ntopics; ++t) conditional[t] /= normalizer;
    new_topic = graphlab::random::multinomial(conditional);      
    ASSERT_LT(new_topic, ntopics);

    // Update the topic assignment and counters
    topics[i] = new_topic;
    if(new_topic != old_topic) nchanges++;
    ++vec_n_dt[new_topic]; 
    ++vec_n_wt[new_topic];
    ++vec_n_t[new_topic];

    if((i+1) % 100000 == 0) {
      std::cout << "(" << n_wt.procid() << ", "
                << "[" << n_wt.cache_hits() << ", " << n_wt.cache_misses() << "] "
                << "[" << n_t.cache_hits() << ", " << n_t.cache_misses() << "])"
                << std::endl;
    }
  } // end of loop over tokens
  n_wt.flush();
  n_t.flush();
  // n_wt.synchronize();
  // n_t.synchronize();
} // end of run gibbs






int main(int argc, char** argv) {
  std::cout << "Running DHT based LDA" << std::endl;
    
  ///! Initialize control plain using mpi
  graphlab::mpi_tools::init(argc, argv);
  graphlab::dc_init_param rpc_parameters;
  graphlab::init_param_from_mpi(rpc_parameters);
  graphlab::distributed_control dc(rpc_parameters);
 
  // Configure command line options
  std::string dictionary_fname("dictionary.txt");
  std::string counts_fname("counts.tsv");
  size_t ntopics = 50;
  size_t niters = 2;
  double alpha(1.0/double(ntopics));
  double beta(0.1);
  graphlab::command_line_options
    clopts("Apply the LDA model to estimate topic "
           "distributions for each document.", true);
  clopts.attach_option("dictionary",
                       &dictionary_fname, dictionary_fname,
                       "Dictionary file");
  clopts.attach_option("counts", 
                       &counts_fname, counts_fname, 
                       "Counts file");
  clopts.attach_option("ntopics", 
                       &ntopics, ntopics, "Number of topics");
  clopts.attach_option("niters",
                       &niters, niters, "Number of iterations");
  clopts.attach_option("alpha",
                       &alpha, alpha, "Alpha prior");
  clopts.attach_option("beta", 
                       &beta, beta, "Beta prior");
  // Parse the command line input
  if(!clopts.parse(argc, argv)) {
    std::cout << "Error in parsing input." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Loading Corpus" << std::endl;
  // Load only the local documents
  corpus corpus(dictionary_fname, counts_fname, 
                dc.procid(), dc.numprocs());
  corpus.shuffle_tokens();
  std::cout << "(" << corpus.ndocs << ", "
            << corpus.nwords << ", "
            << corpus.ntokens << ")" 
            << std::endl;
  
  // Initialize the topic assignments
  std::vector<topic_id_type> topic_asgs(corpus.ntokens, NULL_TOPIC);
  // Initialize the per-document topic counts
  boost::unordered_map<doc_id_type, topic_vector> n_dt(corpus.ndocs);
  

  std::cout << "Construct DHTs" << std::endl;
  // Initialize the shared word and topic counts
  graphlab::delta_dht<word_id_type, topic_vector> n_wt(dc, corpus.nwords + 10);
  graphlab::delta_dht<char, topic_vector > n_t(dc, 2);
  n_t.fresh_predicate().max_uses = 10000;
  n_wt.fresh_predicate().max_uses = 1000;
  // Initialize the topic vector
  n_t[0].resize(ntopics,0);
  // for(size_t t = 0; t < ntopics; ++t) n_t[t] = 0;


  // Run the gibbs sampler
  for(size_t iteration = 0; iteration < niters; ++iteration) {
    std::cout << "Beginning Iteraton " << iteration << std::endl;
    run_gibbs(niters, ntopics, alpha, beta, corpus, 
              topic_asgs, n_dt, n_wt, n_t);
  }
  
  // Wait for all communication to finish
  dc.full_barrier();







  graphlab::mpi_tools::finalize();
  return EXIT_SUCCESS;
}


