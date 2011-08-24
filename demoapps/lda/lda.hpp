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

class vertex_data : public graphlab::idivisible<vertex_data>  {
  vertex_type m_type;
  std::vector<count_type> m_nt;
  std::vector<count_type> m_delta;
public:
  vertex_data(vertex_type type = DOCUMENT) : type(type) { }

  vertex_data split() const { 
    vertex_data vdata = *this;
    vdata.m_delta.
  }

  void set_type(vertex_type new_type) { type = newtype; }
  void set_ntopics(size_t ntopics) {
    m_nt.clear(); m_nt.resize(ntopics, 0);
  }

  vertex_type& type() const { return type; }

  count_type get(const topic_id_type topic) const { 
    ASSERT_LT(topic < m_nt.size());
    if(m_delta.empty()) {
      return m_nt[topic];
    } else {
      ASSERT_LT(topic, m_delta.size());
      return m_nt[topic] + m_delta[topic];
    }
  } // end of get

  void set(const topic_id_type topic, count_type count) { 
    ASSERT_LT(topic < m_nt.size());
    if(m_delta.empty()) {
      m_nt[topic] = count;
    } else {
      ASSERT_LT(topic < m_delta.size());    
      m_delta[topic] += (count - m_nt[topic]);
    }
  } // end of set

  void add(const topic_id_type topic, count_type count) { 
    ASSERT_LT(topic < m_nt.size());
    if(m_delta.empty()) {
      m_nt[topic] += count;
    } else {
      ASSERT_LT(topic < m_delta.size());    
      m_delta[topic] += count;
    }
  } // end of set
};


struct edge_data { 
  count_type count;
  std::vector<count_type> n_t;
  edge_data() : count(0) { }
};


typedef graphlab::graph<vertex_data, edge_data> graph_type;
class lda_update;
typedef graphlab::types<graph_type, lda_update> gl;







#endif
