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

typedef uint32_t count_type;
typedef uint16_t topic_id_type;



/**
 * The vertex data type
 */
struct vertex_data {
  static size_t NTOPICS;
  std::vector< graphlab::atomic<count_type> > topic_count;
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc); 
}; // end of vertex_data


/**
 * The edge data type
 */
struct edge_data {
  count_type ntokens;
  std::vector<topic_id_type> topic_assignment;
  edge_data(count_type ntokens = 0);
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);
}; // end of edge data


/**
 * \brief The graph type;
 */
typedef graphlab::distributed_graph<vertex_data, edge_data> graph_type;



std::pair<edge_data, vertex_data> edge_pair_type;
typedef std::map<graphlab::vertex_id_type, edge_pair_type> neighborhood_map_type


/**
 * \brief The gather type used to accumulate information about the
 * words in a document
 */
struct gather_type {   
  neighborhood_map_type neighborhood_map;
  std::vector<count_type> topic_count;
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);
  gather_type& operator+=(const gather_type& other);
}; // end of gather type


class cgs_lda_vertex_program :
  public graphlab::ivertex_program<graph_type, gather_type> {
private:
  gather_type total;
public:
  edge_dir_type gather_edges(icontext_type& context, 
                             const vertex_type& vertex) const;
  gather_type gather(icontext_type& context, const vertex_type& vertex, 
                     edge_type& edge) const;
  void apply(icontext_type& context, vertex_type& vertex,
             const gather_type& sum);
  edge_dir_type scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const;
  void scatter(icontext_type& context, const vertex_type& vertex, 
               edge_type& edge) const;
  void save(graphlab::oarchive& arc) const;
  void load(graphlab::iarchive& arc);
private:
  bool is_doc(const vertex_type& vertex) const;
  bool is_word(const vertex_type& vertex) const;
  /** Since the edges are undirected we use this helper function to
      get the vertex on the other side of the edge */
  vertex_type get_other_vertex(edge_type& edge, const vertex_type& vertex) const;
}; // end of cgs_lda_vertex_program

