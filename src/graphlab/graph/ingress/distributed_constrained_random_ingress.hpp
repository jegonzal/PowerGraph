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

#ifndef GRAPHLAB_DISTRIBUTED_CONSTRAINED_RANDOM_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_CONSTRAINED_RANDOM_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/ingress/sharding_constraint.hpp>
#include <graphlab/graph/ingress/ingress_edge_decision.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges using randoming hash function.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_constrained_random_ingress : 
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData edge_data_type;


    typedef distributed_ingress_base<VertexData, EdgeData> base_type;

    sharding_constraint* constraint;
    boost::hash<vertex_id_type> hashvid;

  public:
    distributed_constrained_random_ingress(distributed_control& dc, graph_type& graph,
                                           const std::string& method) :
    base_type(dc, graph) {
      constraint = new sharding_constraint(dc.numprocs(), method);
    } // end of constructor

    ~distributed_constrained_random_ingress() { 
      delete constraint;
    }

    /** Add an edge to the ingress object using random assignment. */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      typedef typename base_type::edge_buffer_record edge_buffer_record;

      const std::vector<procid_t>& candidates = constraint->get_joint_neighbors(graph_hash::hash_vertex(source) % base_type::rpc.numprocs(),
                                                                                graph_hash::hash_vertex(target) % base_type::rpc.numprocs());

      const procid_t owning_proc = 
          base_type::edge_decision.edge_to_proc_random(source, target, candidates);


      const edge_buffer_record record(source, target, edata);
      base_type::edge_exchange.send(owning_proc, record);
    } // end of add edge
  }; // end of distributed_constrained_random_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
