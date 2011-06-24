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


#ifndef GRAPHLAB_GENERAL_SCOPE_HPP
#define GRAPHLAB_GENERAL_SCOPE_HPP

#include <boost/bind.hpp>

#include <graphlab/scope/iscope.hpp>
#include <graphlab/scope/iscope_factory.hpp>



namespace graphlab {
  template<typename Graph> class general_scope_factory;

  /**
   * This defines a general scope type 
   */
  template<typename Graph>
  class general_scope :
    public iscope<Graph> {
  public:
    typedef iscope<Graph> base;
    typedef typename Graph::vertex_data_type vertex_data_type;
    typedef typename Graph::edge_data_type edge_data_type;
    typedef typename Graph::edge_list_type   edge_list_type;

    using base::_vertex;
    using base::_graph_ptr;

    scope_range::scope_range_enum stype;
    iscope_factory<Graph>* factory;
  public:
    general_scope() :
      base(NULL,NULL) { }

    general_scope(Graph* graph_ptr, vertex_id_t vertex,
                  iscope_factory<Graph>* factory,
                  scope_range::scope_range_enum s = scope_range::USE_DEFAULT) :
      base(graph_ptr, vertex), stype(s), factory(factory)  {
    }

    scope_range::scope_range_enum scope_type() const {
      return stype;
    }

    ~general_scope() { }

    void commit() {}
    

    void init(Graph* graph, vertex_id_t vertex) {
      base::_graph_ptr = graph;
      base::_vertex = vertex;
    }

    vertex_data_type& vertex_data() {
      return (_graph_ptr->vertex_data(_vertex));
    }


    const vertex_data_type& vertex_data() const {
      return const_vertex_data();
    }

    const vertex_data_type& const_vertex_data() const {
      return (_graph_ptr->vertex_data(_vertex));
    }

    /// Direct calls to access edge data
    const edge_data_type& edge_data(edge_id_t eid) const {
      return const_edge_data(eid);
    }

    const edge_data_type& const_edge_data(edge_id_t eid) const {
      return (_graph_ptr->edge_data(eid));
    }
    
    edge_data_type& edge_data(edge_id_t eid) {
      return (_graph_ptr->edge_data(eid));
    }


    const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {
      return const_neighbor_vertex_data(vertex);
    }


    const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {
      return _graph_ptr->vertex_data(vertex);
    }

    // warning. Guarantee free!
    vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
      return _graph_ptr->vertex_data(vertex);
    }
    
    bool experimental_scope_upgrade(scope_range::scope_range_enum newrange) { 
      assert(factory != NULL);
      factory->release_scope(this);
      ASSERT_TRUE(factory->get_scope(thread::thread_id(),base::_vertex,newrange) == this);
      stype = newrange;
      return true;
    }

    friend class general_scope_factory<Graph>;
  }; // end of ext_locked_scope

} // end of graphlab namespace


#endif

