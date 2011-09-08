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




enum vertex_type {DOCUMENT, WORD};

#ifdef DIFFABLE
class vertex_data : public graphlab::idiffable<vertex_data>  {
#else
class vertex_data {
#endif

  vertex_type m_type;
  count_type m_iterations;
  std::vector<count_type> m_nt;
public:
  vertex_data(vertex_type type = DOCUMENT) : m_type(type), m_iterations(0) { }
 

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
class lda_update;
typedef graphlab::types<graph_type, lda_update> gl;


extern size_t ntopics;
extern size_t nwords;
extern size_t ndocs;
extern double alpha;
extern double beta;
extern std::vector< graphlab::glshared<count_type> > global_n_t;









#endif
