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

#ifndef GRAPHLAB_DISTRIBUTED_BIPARTITE_HASH_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BIPARTITE_HASH_INGRESS_HPP

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
   * \brief Ingress object assigning edges using randoming hash function.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_bipartite_hash_ingress : 
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData edge_data_type;


    typedef distributed_ingress_base<VertexData, EdgeData> base_type;
    bool source_is_special;
  public:
    distributed_bipartite_hash_ingress(distributed_control& dc, graph_type& graph,
                                           const std::string& specialvertex) :
    base_type(dc, graph) {
      if(specialvertex=="source")
        source_is_special=true;
      else
        source_is_special=false;
    } // end of constructor

    ~distributed_bipartite_hash_ingress() { 
      
    }

    /** Add an edge to the ingress object using source or target assignment. */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      typedef typename base_type::edge_buffer_record edge_buffer_record;
      if(source_is_special){
        const procid_t owning_proc = graph_hash::hash_vertex(source) % base_type::rpc.numprocs();
        const edge_buffer_record record(source, target, edata);
        base_type::edge_exchange.send(owning_proc, record);
      }
      else{
        const procid_t owning_proc = graph_hash::hash_vertex(target) % base_type::rpc.numprocs();
        const edge_buffer_record record(source, target, edata);
        base_type::edge_exchange.send(owning_proc, record);
      }
    } // end of add edge
  }; // end of distributed_bipartite_hash_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
