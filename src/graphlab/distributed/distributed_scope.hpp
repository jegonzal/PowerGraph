#ifndef DISTRIBUTED_SCOPE_HPP
#define DISTRIBUTED_SCOPE_HPP

#include <set>
#include <vector>
#include <cassert>

#include <graphlab/scope/iscope.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {



  template<typename Graph, typename DLockManager>
  class distributed_scope: public iscope<Graph> {
  public:
    typedef iscope<Graph> base;
    typedef typename Graph::vertex_data_type vertex_data_type;
    typedef typename Graph::edge_data_type edge_data_type;

    using base::_vertex;
    using base::_graph_ptr;
    DLockManager &lock_manager;
    scope_range::scope_range_enum stype;
    size_t deferred_lock_id;
  public:    
    distributed_scope(DLockManager &lock_manager) :
      base(NULL,-1), 
      lock_manager(lock_manager), 
      stype(scope_range::USE_DEFAULT),
      deferred_lock_id(-1) { }
     
    ~distributed_scope() { }
    
    void commit() {
      if (deferred_lock_id != size_t(-1)) {
        lock_manager.release(deferred_lock_id);
      }
      deferred_lock_id = -1;
    }
    /**
      Called only when the deferred lock has been acquired
    */
    void init(Graph* graph, vertex_id_t vertex, size_t deferred_lock = -1) {
      base::_graph_ptr = graph;
      base::_vertex = vertex;
      deferred_lock_id = deferred_lock;
    }

    vertex_data_type& vertex_data() {
      return lock_manager.get_vertex(base::_vertex);
    }


    const vertex_data_type& vertex_data() const {
      return lock_manager.get_const_vertex(base::_vertex);
    }

    const vertex_data_type& const_vertex_data() const {
      return lock_manager.get_const_vertex(base::_vertex);
    }

    /// Direct calls to access edge data
    const edge_data_type& edge_data(edge_id_t eid) const {
       // TODO: remove this... relies on graph partitioning assumptions.
//       if (base::_graph_ptr->has_constant_edges())  return base::_graph_ptr->edge_data(eid);
      return lock_manager.get_const_edge(eid);
    }

    /// Direct calls to access edge data
    const edge_data_type& const_edge_data(edge_id_t eid) const {
       // TODO: remove this... relies on graph partitioning assumptions.
       //if (base::_graph_ptr->has_constant_edges())  return base::_graph_ptr->edge_data(eid);
      return lock_manager.get_const_edge(eid);
    }
    
    edge_data_type& edge_data(edge_id_t eid) {
      // TODO: remove this... relies on graph partitioning assumptions.
      //if (base::_graph_ptr->has_constant_edges())  return base::_graph_ptr->edge_data(eid);
      return lock_manager.get_edge(eid);
    }


    const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {    
      return lock_manager.get_const_vertex(vertex);
    }

    const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {    
      return lock_manager.get_const_vertex(vertex);
    }

    // warning. Guarantee free!
    vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
      return lock_manager.get_vertex(vertex);
    }
  }; // end of distributed_scope
  
} // end of namespace
#include <graphlab/macros_undef.hpp>


#endif