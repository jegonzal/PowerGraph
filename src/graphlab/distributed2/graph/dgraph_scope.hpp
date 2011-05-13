/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GRAPHLAB_GRAPH_SCOPE_HPP
#define GRAPHLAB_GRAPH_SCOPE_HPP
#include <graphlab/scope/iscope.hpp>

namespace graphlab {
template <typename Graph>
class dgraph_scope : public iscope<Graph> {
 public:
  typedef iscope<Graph> base;
  typedef typename Graph::vertex_data_type vertex_data_type;
  typedef typename Graph::edge_data_type   edge_data_type;
  typedef typename Graph::edge_list_type   edge_list_type;
  
  using base::_vertex;
  using base::_graph_ptr;
  
  scope_range::scope_range_enum stype;
  
  bool own_modified;    // if current vertex is modified
  bool owned_outedges_modified; // if owned outedges are modified
  bool remote_outedges_modified;  // if outgoing edges to remote machines are modified
  bool inedges_modified;  // if incoming edges into current vertex is modifed
  bool owned_nbr_vertices_modified;  // if owned neighboring vertices are modified
  bool remote_nbr_vertices_modified;  // if ghost vertices are modified
  
 public:
  dgraph_scope(): base(NULL, 0) { 
    reset_tracking();
  }
  
  dgraph_scope(Graph* graph_ptr, vertex_id_t vertex,
                scope_range::scope_range_enum s = scope_range::USE_DEFAULT) :
      base(graph_ptr, vertex), stype(s) { 
    reset_tracking();    
  }


  scope_range::scope_range_enum scope_type() const {
    return stype;
  }

  ~dgraph_scope() { }

  bool ghost_sync_required() const {
    return remote_outedges_modified ||
        remote_nbr_vertices_modified || 
        (inedges_modified && _graph_ptr->globalvid_to_replicas(_vertex).size() > 1) ;
  }
  
  void reset_tracking() {
    own_modified = false; 
    owned_outedges_modified = false; 
    remote_outedges_modified = false; 
    inedges_modified = false; 
    owned_nbr_vertices_modified = false; 
    remote_nbr_vertices_modified = false; 
  }

  void commit_impl(bool async, bool untracked ) {
    if (owned_nbr_vertices_modified ||
        owned_outedges_modified) {
        push_owned(async, untracked);
    }
    else if (own_modified) {
      push_current_vertex(async, untracked);  
    }
    
    if (ghost_sync_required()) {
      commit_ghosts(async);
    }

    reset_tracking();
  }

  void commit() {
    commit_impl(false, false);
  }
  void commit_async() {
    commit_impl(true, false);
  }
  
  void commit_async_untracked() {
    commit_impl(true, true);
  }
  
  void commit_ghosts(bool async) { 
    _graph_ptr->synchronize_scope(_vertex, async);
  }
  
  /**
  Push all the owned data in the scope
  */
  void push_owned(bool async, bool untracked) {
    _graph_ptr->push_owned_scope_to_replicas(_vertex,
                                             true, // modified only 
                                             true, // clear modified
                                             async,
                                             untracked);
  }

  /**
  Push the current vertex 
  */
  void push_current_vertex(bool async, bool untracked) {
    _graph_ptr->push_owned_vertex_to_replicas(_vertex,
                                               async,
                                               untracked);
  }
  
  
  void init(Graph* graph, vertex_id_t vertex) {
    base::_graph_ptr = graph;
    base::_vertex = vertex;
    reset_tracking();
  }


  vertex_data_type& vertex_data() {
    own_modified = true;
    _graph_ptr->vertex_is_modified(_vertex);
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
    vertex_id_t etarget = _graph_ptr->target(eid);
    if (etarget == _vertex)  {
      inedges_modified = true;
    }
    else if (_graph_ptr->is_owned(etarget) == false){
      remote_outedges_modified = true;
    }
    else {
      owned_outedges_modified = true;
    }
    _graph_ptr->edge_is_modified(eid);
    return (_graph_ptr->edge_data(eid));
  }


  const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {
    return const_neighbor_vertex_data(vertex);
  }


  const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {
    return _graph_ptr->vertex_data(vertex);
  }


  vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
    if (_graph_ptr->is_owned(vertex) == false) {
      remote_nbr_vertices_modified = true;
    }
    else {
      owned_nbr_vertices_modified = true;
    }
    _graph_ptr->vertex_is_modified(vertex);
    return _graph_ptr->vertex_data(vertex);
  }
  
  bool experimental_scope_upgrade(scope_range::scope_range_enum newrange) { 
    return false;
  }

};


} // end namespace graphlab
#endif

