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


#ifndef LDA_HPP
#define LDA_HPP


#include <graphlab.hpp>
#include "corpus.hpp"




#include <graphlab/macros_def.hpp>


extern size_t ntopics;
extern size_t nwords;
extern size_t ndocs;
extern double alpha;
extern double beta;
extern size_t global_lag;





enum vertex_type {DOCUMENT, WORD};


class vertex_data {
public:
  vertex_type type;
  std::vector<count_type> n_t;
  vertex_data(vertex_type type = DOCUMENT) : type(type) { }
  size_t lag() { return global_lag; }
  void apply_diff(const vertex_data& changed, 
                  const vertex_data& old) {
    ASSERT_EQ(n_t.size(), changed.n_t.size());
    ASSERT_EQ(n_t.size(), old.n_t.size());
    ASSERT_EQ(type, changed.type);
    ASSERT_EQ(type, old.type);
    for(size_t t = 0; t < n_t.size(); ++t) 
      n_t[t] += (changed.n_t[t] - old.n_t[t]);
  } // end of apply diff 
};

  
struct edge_data {
  count_type count;
  std::vector<count_type> n_t; 
  edge_data(count_type count = 0) : count(count) { }  
};

  
typedef graphlab::graph<vertex_data, edge_data> graph_type;



class lda_update : 
  public graphlab::iupdate_functor<graph_type, lda_update> {
public:
  typedef graphlab::iupdate_functor<graph_type, lda_update> base;
  static bool use_factorized; 
public: 
  size_t iters_remaining;
private: 
  std::vector<count_type> delta_n_t;
  

public:

  lda_update(const size_t iters = 100) : iters_remaining(iters) { } 
    
  void operator+=(const lda_update& other) { }
  bool is_factorizable() const { return use_factorized; }
  edge_set gather_edges() { return graphlab::OUT_EDGES; }
  consistency_model gather_consistency() const { 
    return graphlab::FULL_CONSISTENCY;
  }
  edge_set scatter_edges() { return graphlab::NO_EDGES; }
 
  void operator()(icontext_type& context) {
    ASSERT_GT(iters_remaining, 0);

    // Update global and local n_t counts
    if(!context.is_local("n_t")) {
      context.add_local("n_t", count_type(0), ntopics);
      context.add_local("old_n_t", count_type(0), ntopics);
      context.add_local("age", size_t(0));
    }
    std::vector<count_type>& n_t = context.get_local_vec<count_type>("n_t");
    if(context.get_local<size_t>("age") > 5) {
      std::vector<count_type>& old_n_t = 
        context.get_local_vec<count_type>("old_n_t");    
      for(size_t t = 0; t < ntopics; ++t) {
        context.increment_global("n_t", n_t[t] - old_n_t[t], t);
        n_t[t] = old_n_t[t] = context.get_global<count_type>("n_t", t);
      }
      context.get_local<size_t>("age") = 0;
    } else context.get_local<size_t>("age")++;
    
    // Get local data structures
    vertex_data& doc = context.vertex_data(); ASSERT_EQ(doc.type, DOCUMENT);
    // Initialize counts if necessary
    if(doc.n_t.size() != ntopics) doc.n_t.resize(ntopics);


    // Loop over the words in the document (encoded by out edges)
    const edge_list_type out_edges = context.out_edges();
    std::vector<double> prob(ntopics); 
    std::vector<count_type> new_n_t(ntopics, 0);
    double normalizer = 0; 
    foreach(edge_type edge, out_edges) {
      // Get the data ---------------------------------------------------------
      const vertex_id_type word_vid = edge.target();
      const word_id_type word_id = word_vid; ASSERT_LT(word_id, nwords);      
      vertex_data& word = context.vertex_data(word_vid);
      ASSERT_EQ(word.type, WORD);
      if(word.n_t.size() != ntopics) word.n_t.resize(ntopics);
      edge_data& edata = context.edge_data(edge);
      if(edata.n_t.size() != ntopics) edata.n_t.resize(ntopics);          

      // Compute the probability table ----------------------------------------      
      for(size_t t = 0; t < ntopics; ++t) {
        // number of tokens from document d with topic t
        const double cav_n_dt = std::max(doc.n_t[t] - edata.n_t[t], 0);
        // number of token with word w and topic t
        const double cav_n_wt = std::max(word.n_t[t] - edata.n_t[t], 0);
        // number of tokens with topic t
        const double cav_n_t = std::max(n_t[t] - edata.n_t[t], 0);
        // compute the final probability
        prob[t] = (alpha + cav_n_dt) * (beta + cav_n_wt) / 
          (beta*nwords + cav_n_t);
        normalizer += prob[t];
      }
      ASSERT_GT(normalizer, 0);
      // Normalize the probability
      for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;
      
      // Draw new topic assignments -------------------------------------------
      for(size_t t = 0; t < ntopics; ++t) new_n_t[t] = 0;
      for(count_type i = 0; i < edata.count; ++i) 
        ++new_n_t[graphlab::random::multinomial(prob)];
      
      // Update all the tables ------------------------------------------------
      for(size_t t = 0; t < ntopics; ++t) {
        const count_type delta = new_n_t[t] - edata.n_t[t];
        doc.n_t[t] += delta; word.n_t[t] += delta; n_t[t] += delta;
        edata.n_t.swap(new_n_t);
      }
    } // end of for loop

    // Reschedule self if necessary
    if(--iters_remaining > 0) context.schedule(context.vertex_id(), *this);    
  } // end of operator()


  /**
   * This is called before the gather and is used to allocate local
   * data structures.  In this case it is the internal counter
   */
  void init_gather(iglobal_context_type& context) { 
    delta_n_t.resize(ntopics, 0); 
    // Update global and local n_t counts
    if(!context.is_local("n_t")) {
      context.add_local("n_t", count_type(0), ntopics);
      context.add_local("old_n_t", count_type(0), ntopics);
      context.add_local("age", size_t(0));
    }
    std::vector<count_type>& n_t = context.get_local_vec<count_type>("n_t");
    if(context.get_local<size_t>("age") > 5) {
      std::vector<count_type>& old_n_t = 
        context.get_local_vec<count_type>("old_n_t");    
      for(size_t t = 0; t < ntopics; ++t) {
        context.increment_global("n_t", n_t[t] - old_n_t[t], t);
        n_t[t] = old_n_t[t] = context.get_global<count_type>("n_t", t);
      }
      context.get_local<size_t>("age") = 0;
    } else context.get_local<size_t>("age")++;

  }


  /**
   * Gather resamples each edge and accumulates the result in new_nt
   */ 
  void gather(icontext_type& context, edge_type edge) {  
    // Get the data ---------------------------------------------------------
    // Get local data structures
    std::vector<count_type>& n_t = context.get_local_vec<count_type>("n_t");
    const vertex_data& doc = context.vertex_data();
    ASSERT_EQ(doc.type, DOCUMENT);

    const vertex_id_type word_vid = edge.target();
    const word_id_type word_id = word_vid; ASSERT_LT(word_id, nwords);      
    vertex_data& word      = context.vertex_data(word_vid);
    ASSERT_EQ(word.type, WORD);
    edge_data& edata = context.edge_data(edge);
    if(edata.n_t.size() != ntopics) edata.n_t.resize(ntopics);          

    std::vector<double> prob(ntopics); 
    std::vector<count_type> new_n_t(ntopics, 0);
    double normalizer = 0; 

    // Compute the probability table ----------------------------------------
    for(size_t t = 0; t < ntopics; ++t) {      
      const count_type doc_n_t = 
        std::max(((doc.n_t.size()==ntopics)?doc.n_t[t]:0)+delta_n_t[t],0);
      // number of tokens from document d with topic t
      const double cav_n_dt = std::max(doc_n_t - edata.n_t[t], 0);
      // number of token with word w and topic t
      const double cav_n_wt = std::max(word.n_t[t] - edata.n_t[t], 0);
      // number of tokens with topic t
      const double cav_n_t =  std::max(n_t[t] - edata.n_t[t], 0);
      // compute the final probability
      prob[t] = (alpha + cav_n_dt) * (beta + cav_n_wt) / 
        (beta*nwords + cav_n_t);
      normalizer += prob[t];
    }
    ASSERT_GT(normalizer, 0);
    // Normalize the probability
    for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;      
    // Draw new topic assignments -------------------------------------------
    for(size_t t = 0; t < ntopics; ++t) new_n_t[t] = 0;
    for(count_type i = 0; i < edata.count; ++i) 
      ++new_n_t[graphlab::random::multinomial(prob)];
    
    // Update all the tables ------------------------------------------------
    for(size_t t = 0; t < ntopics; ++t) {
      const count_type delta = new_n_t[t] - edata.n_t[t];
      delta_n_t[t] += delta; word.n_t[t] += delta; n_t[t] += delta;
      edata.n_t.swap(new_n_t);
    }
  } // end of gather

  /**
   * Merge the topic counters
   */
  void merge(const lda_update& other) {
    ASSERT_EQ(delta_n_t.size(), ntopics);
    for(size_t t = 0; t < other.delta_n_t.size(); ++t) 
      delta_n_t[t] += other.delta_n_t[t];
  } // end of merge


  /**
   * Update the document count 
   */
  void apply(icontext_type& context) {
    ASSERT_EQ(delta_n_t.size(), ntopics);
    vertex_data& vdata = context.vertex_data();
    for(size_t t = 0; t < ntopics; ++t) vdata.n_t[t] += delta_n_t[t];
  } // end of apply

 
}; // end of lda_update



#include <graphlab/macros_undef.hpp>


#endif
