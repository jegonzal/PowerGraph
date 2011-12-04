/*  
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


#ifndef GRAPHLAB_DISTRIBUTED_DGRAPH_CONTEXT_HPP
#define GRAPHLAB_DISTRIBUTED_DGRAPH_CONTEXT_HPP
#include <graphlab/context/iglobal_context.hpp>
#include <graphlab/context/icontext.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {
  
template <typename Engine>
class dgraph_context: public icontext<typename Engine::graph_type, 
                                     typename Engine::update_functor_type> {
 public:

  typedef Engine engine_type;

  typedef typename engine_type::graph_type graph_type;
  typedef typename engine_type::update_functor_type update_functor_type;
  typedef typename engine_type::shared_data_type shared_data_type;

  typedef typename graph_type::vertex_id_type      vertex_id_type;
  typedef typename graph_type::vertex_color_type   vertex_color_type;
  typedef typename graph_type::edge_id_type        edge_id_type;
  typedef typename graph_type::edge_type edge_type;
  typedef typename graph_type::vertex_data_type    vertex_data_type;
  typedef typename graph_type::edge_data_type      edge_data_type;
  typedef typename graph_type::edge_list_type      edge_list_type;
  typedef typename graph_type::vertex_list_type    vertex_list_type;
  

  /** a pointer to the engine */
  engine_type* engine_ptr;
  /** A pointer to the underlying graph datastructure */
  graph_type* graph_ptr;   
  /** A pointer to the scheduler */
  shared_data_type* shared_data_ptr;
  /** the cpuid of this context */
  size_t _threadid;

  
  /** The vertex that this graph represents*/
  vertex_id_type _vertex;

  /** The consistency model that this context ensures */
  consistency_model::model_enum _consistency;


  bool own_modified;    // if current vertex is modified
  bool owned_outedges_modified; // if owned outedges are modified
  bool remote_outedges_modified;  // if outgoing edges to remote machines are modified
  bool inedges_modified;  // if incoming edges into current vertex is modifed
  bool owned_nbr_vertices_modified;  // if owned neighboring vertices are modified
  bool remote_nbr_vertices_modified;  // if ghost vertices are modified

  
 public:

  dgraph_context(engine_type* engine_ptr = NULL,
                 graph_type* graph_ptr = NULL,
                 shared_data_type* shared_data_ptr = NULL):
      engine_ptr(engine_ptr), graph_ptr(graph_ptr), 
      shared_data_ptr(shared_data_ptr),
      _consistency(consistency_model::EDGE_CONSISTENCY) { 
    reset_tracking();    
  }


  //! Get the consistency model under which this context was acquired
  consistency_model::model_enum consistency() const { 
    return _consistency;
  }
    

  ~dgraph_context() { }

  bool ghost_sync_required() const {
    return remote_outedges_modified ||
        remote_nbr_vertices_modified || 
        (inedges_modified && graph_ptr->globalvid_to_replicas(_vertex).size() > 1) ;
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
    } else if (own_modified) {
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
    graph_ptr->synchronize_scope(_vertex, async);
  }
  
  /**
  Push all the owned data in the scope
  */
  void push_owned(bool async, bool untracked) {
    graph_ptr->push_owned_scope_to_replicas(_vertex,
                                             true, // modified only 
                                             true, // clear modified
                                             async,
                                             untracked);
  }

  /**
  Push the current vertex 
  */
  void push_current_vertex(bool async, bool untracked) {
    graph_ptr->push_owned_vertex_to_replicas(_vertex,
                                               async,
                                               untracked);
  }
  
  
  void init(vertex_id_type vertex, size_t threadid) {
    _vertex = vertex;
    _threadid = threadid;
    reset_tracking();
  }


  vertex_data_type& vertex_data() {
    own_modified = true;
    graph_ptr->vertex_is_modified(_vertex);
    return (graph_ptr->vertex_data(_vertex));
  }


  const vertex_data_type& vertex_data() const {
    return const_vertex_data();
  }

  const vertex_data_type& const_vertex_data() const {
    return (graph_ptr->vertex_data(_vertex));
  }

  /// Direct calls to access edge data
  const edge_data_type& edge_data(const edge_type& eid) const {
    return const_edge_data(eid);
  }

  const edge_data_type& const_edge_data(const edge_type& eid) const {
    return (graph_ptr->edge_data(eid.eid()));
  }
  
  edge_data_type& edge_data(const edge_type& eid) {
    vertex_id_type etarget = eid.target();
    if (etarget == _vertex) {
      inedges_modified = true;
    } else if (graph_ptr->is_owned(etarget) == false){
      remote_outedges_modified = true;
    } else {
      owned_outedges_modified = true;
    }
    graph_ptr->edge_is_modified(eid.eid());
    return (graph_ptr->edge_data(eid.eid()));
  }


  edge_data_type& edge_data(vertex_id_type source, vertex_id_type target) {
    if (target == _vertex) {
      inedges_modified = true;
    } else if (graph_ptr->is_owned(target) == false){
      remote_outedges_modified = true;
    } else {
      owned_outedges_modified = true;
    }
    graph_ptr->edge_is_modified(graph_ptr->find(source, target).eid());

    return graph_ptr->edge_data(source, target);
  }

  const edge_data_type& const_edge_data(vertex_id_type source, 
                                        vertex_id_type target) const {
    return graph_ptr->edge_data(source, target);
  }
  
  const edge_data_type& edge_data(vertex_id_type source, 
                                  vertex_id_type target) const {
    return graph_ptr->edge_data(source, target);
  }

  const vertex_data_type& vertex_data(vertex_id_type vertex) const {
    return const_vertex_data(vertex);
  }

  const vertex_data_type& const_vertex_data(vertex_id_type vertex) const {
    return graph_ptr->vertex_data(vertex);
  }

  vertex_data_type& vertex_data(vertex_id_type vertex) {
    if (vertex == _vertex) return vertex_data();
    if (graph_ptr->is_owned(vertex) == false) {
      remote_nbr_vertices_modified = true;
    } else {
      owned_nbr_vertices_modified = true;
    }
    graph_ptr->vertex_is_modified(vertex);
    return graph_ptr->vertex_data(vertex);
  }

  vertex_color_type color() const { 
    return graph_ptr->get_color(_vertex); 
  }

  vertex_color_type color(vertex_id_type vid) const { 
    return graph_ptr->get_color(vid); 
  }

  vertex_id_type vertex_id() const { 
    return _vertex; 
  }
  
  edge_id_type edge(vertex_id_type source,
                    vertex_id_type target) const {
    return graph_ptr->edge_id(source, target);
  }

  bool edge_exists(vertex_id_type source,
                    vertex_id_type target) const {
    return graph_ptr->find(source, target).first;
  }

  edge_type reverse_edge(const edge_type& eid) const {      
    return graph_ptr->reverse_edge(eid);
  }


  virtual edge_type find(vertex_id_type source,
                          vertex_id_type target) const  {
    return graph_ptr->find(source, target);
  }


  edge_list_type in_edges() const {
    return graph_ptr->in_edges(_vertex);
  }

  edge_list_type in_edges(vertex_id_type v) const {
    return graph_ptr->in_edges(v);
  }


  size_t num_in_edges(vertex_id_type v) const {
    return graph_ptr->num_in_neighbors(v);
  }


  size_t num_out_edges(vertex_id_type v) const {
    return graph_ptr->num_out_neighbors(v);
  }


  size_t num_in_edges() const {
    return graph_ptr->num_in_neighbors(_vertex);
  }


  size_t num_out_edges() const {
    return graph_ptr->num_out_neighbors(_vertex);
  }

  edge_list_type out_edges() const {
    return graph_ptr->out_edges(_vertex);
  }
  
  edge_list_type out_edges(vertex_id_type v) const {
    return graph_ptr->out_edges(v);
  }

  std::vector<vertex_id_type> in_vertices() const {
    return graph_ptr->in_vertices(_vertex);
  }

  vertex_list_type in_vertices_list() const {
    return graph_ptr->in_vertices_list(_vertex);
  }

  std::vector<vertex_id_type> out_vertices() const {
    return graph_ptr->out_vertices(_vertex);
  }

  vertex_list_type out_vertices_list() const {
    return graph_ptr->out_vertices_list(_vertex);
  }

  //! Get the source vertex of the edge id argument
  vertex_id_type source(edge_id_type edge_id) const {
    return graph_ptr->source(edge_id);
  }

  //! get the target vertex of the edge id argument
  vertex_id_type target(edge_id_type edge_id) const {
    return graph_ptr->target(edge_id);
  }
  
  void acquire_lock(const std::string& key, size_t index = 0) { 
    shared_data_ptr->acquire_lock(key, index); 
  }

  void release_lock(const std::string& key, size_t index = 0) { 
    shared_data_ptr->release_lock(key, index); 
  }

  void commit_change(const std::string& key, size_t index = 0) { 
    shared_data_ptr->commit_change(key, index);
  }

  void get_global(const std::string& key,                                       
                  graphlab::any_vector*& ret_vec_ptr,
                  bool& ret_is_const) {
    shared_data_ptr->get_global(key, ret_vec_ptr, ret_is_const);
  } // end of get_any_pair    



  void schedule(const vertex_id_type& vertex, 
                const update_functor_type& update_fun) {
    engine_ptr->schedule_from_context(_threadid, vertex, update_fun);
  }

  void schedule_in_neighbors(const vertex_id_type& vertex, 
                              const update_functor_type& update_fun) {
    const edge_list_type edges = graph_ptr->in_edges(vertex);
    foreach(const edge_type& edge, edges) {
      schedule(edge.source(), update_fun);
    }
  }

  void schedule_out_neighbors(const vertex_id_type& vertex, 
                              const update_functor_type& update_fun) {
    const edge_list_type edges = graph_ptr->out_edges(vertex);
    foreach(const edge_type& edge, edges) {
      schedule(edge.target(), update_fun);
    }
  }

  void schedule_neighbors(const vertex_id_type& vertex, 
                          const update_functor_type& update_fun) {
    schedule_in_neighbors(vertex, update_fun);
    schedule_out_neighbors(vertex, update_fun);
  }

  /**
    * Get the number of vertices in the graph.
    */
  size_t num_vertices() const {
    return graph_ptr->num_vertices();
  }

  /**
    * Get the number of edges in the graph
    */
  size_t num_edges() const {
    return graph_ptr->num_edges();
  }

  /**
    * Get an estimate of the number of update functions executed up
    * to this point.
    */
  size_t num_updates() const {
    return engine_ptr->last_update_count();
  }

  void terminate() {
    engine_ptr->stop();
  }
 
  float elapsed_time() const { return engine_ptr->elapsed_time(); }
 
  
};


} // end namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif

