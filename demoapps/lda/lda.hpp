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

#include <graphlab/shared_data/sharedsum.hpp>
#include <graphlab.hpp>
#include "corpus.hpp"




#include <graphlab/macros_def.hpp>


extern size_t ntopics;
extern size_t nwords;
extern size_t ndocs;
extern double alpha;
extern double beta;
extern size_t global_lag;
extern std::vector< graphlab::glshared<count_type> > global_n_t;
extern std::vector< graphlab::sharedsum<count_type> > shared_n_t;




enum vertex_type {DOCUMENT, WORD};

#ifdef DIFFABLE
class vertex_data : public graphlab::idiffable<vertex_data>  
#else
class vertex_data 
#endif
{

  vertex_type m_type;
  count_type m_iterations;
  std::vector<count_type> m_nt;
public:
  vertex_data(vertex_type type = DOCUMENT) : m_type(type), m_iterations(0) { }
 
  size_t lag() { return global_lag; }

  void set_type(vertex_type new_type) { m_type = new_type; }
  const vertex_type& type() const { return m_type; }

  bool is_finished() const { return m_iterations <= 0; }
  void finished_iteration() { m_iterations--; }


  void init(size_t ntopics, count_type iterations) {
    if(m_nt.size() != ntopics) {
      m_nt.clear(); m_nt.resize(ntopics, 0);
    }
    m_iterations = iterations;
  } //end of init


  void apply_diff(const vertex_data& changed, 
                  const vertex_data& old) {
    ASSERT_EQ(m_nt.size(), changed.m_nt.size());
    ASSERT_EQ(m_nt.size(), old.m_nt.size());
    ASSERT_EQ(m_type, changed.m_type);
    ASSERT_EQ(m_type, old.m_type);
    //    std::cout << "Running diff" << std::endl;
    for(size_t t = 0; t < m_nt.size(); ++t) 
      m_nt[t] += (changed.m_nt[t] - old.m_nt[t]);
    m_iterations += (changed.m_iterations - old.m_iterations);
  } // end of apply diff 


  count_type get(const topic_id_type topic) const { 
    ASSERT_LT(topic, m_nt.size());
    return m_nt[topic];
  } // end of get

  void set(const topic_id_type topic, count_type count) { 
    ASSERT_LT(topic, m_nt.size());
    m_nt[topic] = count;
  } // end of set

  void add(const topic_id_type topic, count_type count) { 
    ASSERT_LT(topic, m_nt.size());
    m_nt[topic] += count;
  } // end of set

  void subtract(const topic_id_type topic, count_type count) { 
    ASSERT_LT(topic, m_nt.size());
    m_nt[topic] -= count;
  } // end of set
};

  
struct edge_data {
  count_type count;
  std::vector<count_type> m_nt; 
public:
  edge_data(count_type count = 0) : count(count) { }
  
  void set_count(const count_type new_count) { 
    count = new_count; 
  }
  
  count_type get_count() const { return count; }

  void init(size_t ntopics) {
    if(m_nt.size() != ntopics) {
      m_nt.clear(); m_nt.resize(ntopics, 0); 
    }
  }

  count_type get(const topic_id_type topic) const { 
    ASSERT_LT(topic, m_nt.size());
    return m_nt[topic];
  }
  void set(const topic_id_type topic, count_type count) {
    ASSERT_LT(topic, m_nt.size());
    m_nt[topic] = count;
  }
};

  
typedef graphlab::graph<vertex_data, edge_data> graph_type;



class lda_update : 
  public graphlab::iupdate_functor<graph_type, lda_update> {
  size_t iters_remaining;
public:
  typedef graphlab::iupdate_functor<graph_type, lda_update> base;
  typedef base::icontext_type icontext_type;

  typedef base::edge_id_type edge_id_type;
  typedef base::edge_list_type edge_list_type;
  typedef base::vertex_id_type vertex_id_type;

  static bool use_factorized;

private:
  std::vector<count_type> old_global_n_t;
  std::vector<count_type> local_n_t;
  

public:

  lda_update(const size_t iters = 100) : iters_remaining(iters) { } 
    
  void operator+=(const lda_update& other) { }
  bool is_factorizable() const { return use_factorized; }
  base::edge_set gather_edges() { return base::OUT_EDGES; }
  bool writable_gather() const { return true; }
  base::edge_set scatter_edges() { return base::NO_EDGES; }
 
  void operator()(icontext_type& context) {
    ASSERT_GT(iters_remaining, 0);
    // Make a local copy of the global topic counts
    std::vector<count_type> local_n_t(ntopics);
#ifdef SHAREDSUM
    for(size_t t = 0; t < ntopics; ++t) 
      local_n_t[t] = shared_n_t[t].val();
#else
    for(size_t t = 0; t < ntopics; ++t) 
      local_n_t[t] = global_n_t[t].get_val();
#endif
    const std::vector<count_type> old_global_n_t = local_n_t;

    // Get local data structures
    vertex_data& doc       = context.vertex_data();
    // only gather on a document
    ASSERT_EQ(doc.type(), DOCUMENT);

    // Loop over the words in the document (encoded by out edges)
    const edge_list_type out_edges = context.out_edge_ids();
    std::vector<double> prob(ntopics); 
    std::vector<count_type> n_dwt(ntopics);
    double normalizer = 0; 
    foreach(edge_id_type eid, out_edges) {
      // Get the data ---------------------------------------------------------
      const vertex_id_type word_vid = context.target(eid);
      const word_id_type word_id   = word_vid;
      vertex_data& word      = context.neighbor_vertex_data(word_vid);
      edge_data& edata       = context.edge_data(eid);
      edata.init(ntopics); // ensure that the edge data is initialized
      ASSERT_LT(word_id, nwords);      
      ASSERT_EQ(word.type(), WORD);

      // Compute the probability table ----------------------------------------
      for(size_t t = 0; t < ntopics; ++t) {
        // number of tokens from document d with topic t
        const double n_dt = 
          std::max(doc.get(t) - edata.get(t), 0);
        // number of token with word w and topic t
        const double n_wt = 
          std::max(word.get(t) - edata.get(t), 0);
        // number of tokens with topic t
        const double n_t = 
          std::max(local_n_t[t] - edata.get(t), 0);
        // compute the final probability
        prob[t] = (alpha + n_dt) * (beta + n_wt) / 
          (beta * nwords + n_t);
        normalizer += prob[t];
      }
      ASSERT_GT(normalizer, 0);
      // Normalize the probability
      for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;
      
      // Draw new topic assignments -------------------------------------------
      n_dwt.clear(); n_dwt.resize(ntopics, 0);
      for(count_type i = 0; i < edata.get_count(); ++i) {
        const size_t topic_id = graphlab::random::multinomial(prob); 
        ASSERT_LT(topic_id, ntopics);
        n_dwt[topic_id]++;
      }
      
      // Update all the tables ------------------------------------------------
      for(size_t t = 0; t < ntopics; ++t) {
        doc.add(t, n_dwt[t] - edata.get(t));
        word.add(t, n_dwt[t] - edata.get(t));
        local_n_t[t] += (n_dwt[t] - edata.get(t));
        edata.set(t, n_dwt[t]);
      }
    } // end of for loop
#ifdef SHAREDSUM
    // update the global variables
    for(size_t t = 0; t < ntopics; ++t) {
      shared_n_t[t] += (local_n_t[t] - old_global_n_t[t]);
    }
#else
    // update the global variables
    for(size_t t = 0; t < ntopics; ++t) {
      global_n_t[t] += (local_n_t[t] - old_global_n_t[t]);
    }
#endif

    // Reschedule self if necessary
    if(--iters_remaining > 0) 
      context.schedule(context.vertex_id(), *this);
    
  } // end of operator()


  void gather(icontext_type& context, edge_id_type eid) {

    if(local_n_t.empty()) {
      local_n_t.resize(ntopics);
      old_global_n_t.resize(ntopics);
      for(size_t t = 0; t < ntopics; ++t) 
        old_global_n_t[t] = local_n_t[t] = global_n_t[t].get_val();
    }
    

    // Get the data ---------------------------------------------------------
    // Get local data structures
    vertex_data& doc       = context.vertex_data();
    // only gather on a document
    ASSERT_EQ(doc.type(), DOCUMENT);

    const vertex_id_type word_vid = context.target(eid);
    const word_id_type word_id   = word_vid;
    vertex_data& word      = context.neighbor_vertex_data(word_vid);
    edge_data& edata       = context.edge_data(eid);
    edata.init(ntopics); // ensure that the edge data is initialized
    ASSERT_LT(word_id, nwords);      
    ASSERT_EQ(word.type(), WORD);

    std::vector<double> prob(ntopics); 
    std::vector<count_type> n_dwt(ntopics);
    double normalizer = 0; 


    // Compute the probability table ----------------------------------------
    for(size_t t = 0; t < ntopics; ++t) {
      // number of tokens from document d with topic t
      const double n_dt = 
        std::max(doc.get(t) - edata.get(t), 0);
      // number of token with word w and topic t
      const double n_wt = 
        std::max(word.get(t) - edata.get(t), 0);
      // number of tokens with topic t
      const double n_t = 
        std::max(local_n_t[t] - edata.get(t), 0);
      // compute the final probability
      prob[t] = (alpha + n_dt) * (beta + n_wt) / 
        (beta * nwords + n_t);
      normalizer += prob[t];
    }
    ASSERT_GT(normalizer, 0);
    // Normalize the probability
    for(size_t t = 0; t < ntopics; ++t) prob[t] /= normalizer;
      
    // Draw new topic assignments -------------------------------------------
    n_dwt.clear(); n_dwt.resize(ntopics, 0);
    for(count_type i = 0; i < edata.get_count(); ++i) {
      const size_t topic_id = graphlab::random::multinomial(prob); 
      ASSERT_LT(topic_id, ntopics);
      n_dwt[topic_id]++;
    }

    // Update all the tables ------------------------------------------------
    for(size_t t = 0; t < ntopics; ++t) {
      doc.add(t, n_dwt[t] - edata.get(t));
      word.add(t, n_dwt[t] - edata.get(t));
      local_n_t[t] += (n_dwt[t] - edata.get(t));
      edata.set(t, n_dwt[t]);
    }

    // Reschedule self if necessary
    if(iters_remaining > 0) 
      context.schedule(context.vertex_id(), lda_update(iters_remaining-1));

  } // end of gather


  
  void apply(icontext_type& context) {
    ASSERT_EQ(global_n_t.size(), ntopics);
    if(!local_n_t.empty()) {
      // update the global variables
      for(size_t t = 0; t < ntopics; ++t) {
        global_n_t[t] += (local_n_t[t] - old_global_n_t[t]);
      }    
    }
  }

 
}; // end of lda_update



#include <graphlab/macros_undef.hpp>


#endif
