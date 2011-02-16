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

 public:
  dgraph_scope(): base(NULL, NULL) { }
  
  dgraph_scope(Graph* graph_ptr, vertex_id_t vertex,
                scope_range::scope_range_enum s = scope_range::USE_DEFAULT) :
      base(graph_ptr, vertex), stype(s) { }

  scope_range::scope_range_enum scope_type() const {
    return stype;
  }

  ~dgraph_scope() { }

  void commit_ghosts() { 
    graph_ptr->synchronize_scope(_vertex);
  }
  
  void commit_ghosts_async() { 
    graph_ptr->synchronize_scope(_vertex, true);
  }

  void push_owned() {
  }
  
  void push_owned_async() { 
  }
  
  void init(Graph* graph, vertex_id_t vertex) {
    base::_graph_ptr = graph;
    base::_vertex = vertex;
  }


  vertex_data_type& vertex_data() {
    _graph_ptr->increment_vertex_version(_vertex);
    _graph_ptr->vertex_modified(_vertex);
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
    _graph_ptr->increment_edge_version(eid);
    _graph_ptr->vertex_modified(eid);
    return (_graph_ptr->edge_data(eid));
  }


  const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {
    return const_neighbor_vertex_data(vertex);
  }


  const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {
    return _graph_ptr->vertex_data(vertex);
  }


  vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
    _graph_ptr->increment_vertex_version(vertex);
    _graph_ptr->vertex_modified(vertex);
    return _graph_ptr->vertex_data(vertex);
  }
  
  bool experimental_scope_upgrade(scope_range::scope_range_enum newrange) { 
    return false;
  }

};


} // end namespace graphlab
#endif
