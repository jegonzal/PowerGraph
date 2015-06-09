/*  
 * Copyright (c) 2013 Shanghai Jiao Tong University. 
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
 *      http://ipads.se.sjtu.edu.cn/projects/powerlyra.html
 *
 *
 * 2014.04  implement bipartite-aware random partitioning
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_BIPARTITE_RANDOM_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BIPARTITE_RANDOM_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges using randoming hash function on favorite.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_bipartite_random_ingress : 
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData edge_data_type;


    typedef distributed_ingress_base<VertexData, EdgeData> base_type;
    
    bool favorite_source;
  public:
    distributed_bipartite_random_ingress(distributed_control& dc, graph_type& graph, const std::string& favorite) :
    base_type(dc, graph) {
      favorite_source = (favorite == "source") ? true : false;
    } // end of constructor

    ~distributed_bipartite_random_ingress() { }

    /** Add an edge to the ingress object using random of "favorite" assignment. */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      typedef typename base_type::edge_buffer_record edge_buffer_record;
      vertex_id_type favorite = favorite_source ? source : target;
      const procid_t owning_proc = graph_hash::hash_vertex(favorite) % base_type::rpc.numprocs();
      const edge_buffer_record record(source, target, edata);
      base_type::edge_exchange.send(owning_proc, record);
    } // end of add edge
  }; // end of distributed_bipartite_random_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
