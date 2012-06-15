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

#ifndef CVB0_LDA_HPP
#define CVB0_LDA_HPP
 
#include <vector>
#include <algorithm>
#include <graphlab.hpp>

#include <graphlab/macros_def.hpp>
typedef float count_type;
typedef uint16_t topic_id_type;
#define NULL_TOPIC (topic_id_type(-1))

typedef std::vector< count_type > factor_type;
typedef std::vector< topic_id_type > assignment_type;

extern double ALPHA;
extern double BETA;
extern size_t NTOPICS;
extern size_t NWORDS;
extern size_t TOPK;
extern size_t INTERVAL;
extern factor_type GLOBAL_TOPIC_COUNT;
extern std::vector<std::string> DICTIONARY;
extern size_t MAX_COUNT;



inline factor_type& operator+=(factor_type& lvalue, const factor_type& rvalue) {
  if(!rvalue.empty()) {
    if(lvalue.size() != NTOPICS) lvalue = rvalue;
    else {
      for(size_t t = 0; t < lvalue.size(); ++t) lvalue[t] += rvalue[t];
    }
  }
  return lvalue;
} // end of operator +=



/**
 * The vertex data type
 */
struct vertex_data {
  factor_type factor;
  size_t nupdates;
  float nchanges;
  vertex_data() : nupdates(0), nchanges(0) { }
  void save(graphlab::oarchive& arc) const { 
    arc << factor << nupdates << nchanges; 
  }
  void load(graphlab::iarchive& arc) { 
    arc >> factor >> nupdates >> nchanges; 
  } 
}; // end of vertex_data


/**
 * The edge data type
 */
struct edge_data {
  uint32_t count; factor_type belief;
  edge_data(uint32_t count = 0) : count(count), belief(NTOPICS) { }
  void save(graphlab::oarchive& arc) const { 
    arc << count << belief;
  }
  void load(graphlab::iarchive& arc) { 
    arc >> count >> belief;
  }
}; // end of edge data



/**
 * \brief The graph type;
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;




bool graph_loader(graph_type& graph, const std::string& fname, 
                  const std::string& line);

inline void initialize_vertex_data(graph_type::vertex_type& vertex) {
  vertex.data().factor.resize(NTOPICS);
}


bool load_and_initialize_graph(graphlab::distributed_control& dc,
                               graph_type& graph,
                               const std::string& matrix_dir);




/** populate the global dictionary */
bool load_dictionary(const std::string& fname);


inline bool is_word(const graph_type::vertex_type& vertex) {
  return vertex.num_in_edges() > 0 ? 1 : 0;
}

inline bool is_doc(const graph_type::vertex_type& vertex) {
  return vertex.num_out_edges() > 0 ? 1 : 0;
}

inline graph_type::vertex_type 
get_other_vertex(graph_type::edge_type& edge, 
                 const graph_type::vertex_type& vertex) {
  return vertex.id() == edge.source().id()? edge.target() : edge.source();
}


template<typename IContext>
class topk_aggregator {
  typedef IContext icontext_type;
  typedef std::pair<float, graphlab::vertex_id_type> cw_pair_type;
private:
  std::vector< std::set<cw_pair_type> > top_words;
  float nchanges, nupdates;
public:
  topk_aggregator(size_t nchanges = 0, size_t nupdates = 0) : 
    nchanges(nchanges), nupdates(nupdates) { }

  void save(graphlab::oarchive& arc) const { 
    arc << top_words << nchanges << nupdates; 
  }
  void load(graphlab::iarchive& arc) { 
    arc >> top_words >> nchanges >> nupdates;
  }

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
        // std::cout << DICTIONARY[pair.second] 
        //           << "(" << pair.first << ")" << ", "; 
        std::cout << DICTIONARY[pair.second] << ",  ";
      }
      std::cout << std::endl;
    }
    std::cout << "\nNumber of token changes: " << total.nchanges << std::endl;
    std::cout << "\nNumber of updates:       " << total.nupdates << std::endl;
  } // end of finalize
}; // end of topk_aggregator struct




template<typename IContext>
struct global_counts_aggregator {
  typedef graph_type::vertex_type vertex_type;
  static factor_type map(IContext& context, const vertex_type& vertex) {
    return vertex.data().factor;
  } // end of map function

  static void finalize(IContext& context, const factor_type& total) {
    for(size_t t = 0; t < total.size(); ++t)
      GLOBAL_TOPIC_COUNT[t] =
        std::max(count_type(total[t]/2), count_type(0));
  } // end of finalize
}; // end of global_counts_aggregator struct




template<typename IContext>
struct selective_signal {
 static graphlab::empty 
 docs(IContext& context, graph_type::vertex_type& vertex) {
   if(is_doc(vertex)) context.signal(vertex);
   return graphlab::empty();
 } // end of signal_docs
 static graphlab::empty 
 words(IContext& context, graph_type::vertex_type& vertex) {
    if(is_word(vertex)) context.signal(vertex);
    return graphlab::empty();
  } // end of signal_words
}; // end of selective_signal








#include <graphlab/macros_undef.hpp>
#endif
