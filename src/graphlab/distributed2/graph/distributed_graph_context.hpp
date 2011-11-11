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


#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_CONTEXT_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_CONTEXT_HPP
#include <graphlab/context/icontext.hpp>

namespace graphlab {

template<typename Engine>
class dgraph_context: public icontext<typename Engine::graph_type, 
                                     typename Engine::update_functor_type> {

 public:
  typedef icontext<typename Engine::graph_type, 
                   typename Engine::update_functor_type> base;
                   
  typedef Engine engine_type;

  typedef typename engine_type::graph_type graph_type;
  typedef typename engine_type::update_functor_type update_functor_type;
  typedef typename engine_type::ischeduler_type ischeduler_type;

  typedef typename graph_type::vertex_id_type      vertex_id_type;
  typedef typename graph_type::vertex_color_type   vertex_color_type;
  typedef typename graph_type::edge_id_type        edge_id_type;
  typedef typename graph_type::vertex_data_type    vertex_data_type;
  typedef typename graph_type::edge_data_type      edge_data_type;
  typedef typename graph_type::edge_list_type      edge_list_type;
  
  /** A pointer to the underlying graph datastructure */
  graph_type* graph_ptr;   

  /** a pointer to the engine */
  engine_type* eng_ptr;

  size_t cpuid;
  
   /** The vertex that this graph represents*/
  vertex_id_type vid;

  /** The consistency model that this context ensures */
  consistency_model::model_enum consistency;

  bool own_modified;    // if current vertex is modified
  bool owned_outedges_modified; // if owned outedges are modified
  bool remote_outedges_modified;  // if outgoing edges to remote machines are modified
  bool inedges_modified;  // if incoming edges into current vertex is modifed
  bool owned_nbr_vertices_modified;  // if owned neighboring vertices are modified
  bool remote_nbr_vertices_modified;  // if ghost vertices are modified
  
 public:
  dgraph_context(graph_type* graph_ptr = NULL,
                 engine_type* eng_ptr = NULL,
                 size_t cpuid = -1):graph_ptr(graph_ptr), 
                                    eng_ptr(eng_ptr), cpuid(cpuid),
                                    vid(-1), consistency(consistency_model::EDGE_CONSISTENCY) {
    reset_tracking();    
  }

  ~dgraph_context() { }

  bool ghost_sync_required() const {
    return remote_outedges_modified ||
        remote_nbr_vertices_modified || 
        (inedges_modified && graph_ptr->globalvid_to_replicas(vid).size() > 1) ;
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
    graph_ptr->synchronize_scope(vid, async);
  }
  
  /**
  Push all the owned data in the scope
  */
  void push_owned(bool async, bool untracked) {
    graph_ptr->push_owned_scope_to_replicas(vid,
                                             true, // modified only 
                                             true, // clear modified
                                             async,
                                             untracked);
  }

  /**
  Push the current vertex 
  */
  void push_current_vertex(bool async, bool untracked) {
    graph_ptr->push_owned_vertex_to_replicas(vid,
                                               async,
                                               untracked);
  }
  
  
  void init(const vertex_id_type vertex, 
            consistency_model::model_enum consistency) {             
    vid = vertex;
    consistency = consistency;
  } 
  
  
  size_t num_vertices() const { 
    return graph_ptr->num_vertices(); 
  }
  size_t num_edges() const { 
    return graph_ptr->num_edges(); 
  }
  void terminate() { 
    eng_ptr->stop(); 
  }



  vertex_data_type& vertex_data() {
    own_modified = true;
    graph_ptr->vertex_is_modified(vid);
    return graph_ptr->vertex_data(vid);
  }


  const vertex_data_type& vertex_data() const {
    return const_vertex_data();
  }

  const vertex_data_type& const_vertex_data() const {
    return graph_ptr->vertex_data(vid);
  }

  /// Direct calls to access edge data
  const edge_data_type& edge_data(edge_id_t eid) const {
    return const_edge_data(eid);
  }

  const edge_data_type& const_edge_data(edge_id_t eid) const {
    return (graph_ptr->edge_data(eid));
  }
  
  edge_data_type& edge_data(edge_id_t eid) {
    vertex_id_t etarget = graph_ptr->target(eid);
    if (etarget == vid)  {
      inedges_modified = true;
    }
    else if (graph_ptr->is_owned(etarget) == false){
      remote_outedges_modified = true;
    }
    else {
      owned_outedges_modified = true;
    }
    graph_ptr->edge_is_modified(eid);
    return (graph_ptr->edge_data(eid));
  }


  const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {
    return const_neighbor_vertex_data(vertex);
  }


  const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {
    return graph_ptr->vertex_data(vertex);
  }


  vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
    if (graph_ptr->is_owned(vertex) == false) {
      remote_nbr_vertices_modified = true;
    }
    else {
      owned_nbr_vertices_modified = true;
    }
    graph_ptr->vertex_is_modified(vertex);
    return graph_ptr->vertex_data(vertex);
  }
  

  vertex_color_type color() const { return graph_ptr->get_color(vid); }
  
  vertex_id_type vertex_id() const { return vid; }

  edge_id_type edge(vertex_id_type source,
                    vertex_id_type target) const {
    return graph_ptr->edge_id(source, target);
  }

  bool edge_exists(vertex_id_type source,
                   vertex_id_type target) const {
    return graph_ptr->find(source, target).first;
  }

  edge_id_type reverse_edge(edge_id_type eid) const {      
    return graph_ptr->rev_edge_id(eid);
  }

  edge_list_type in_edge_ids() const {
    return graph_ptr->in_edge_ids(vid);
  }

  edge_list_type in_edge_ids(vertex_id_type v) const {
    return graph_ptr->in_edge_ids(v);
  }

  edge_list_type out_edge_ids() const {
    return graph_ptr->out_edge_ids(vid);
  }

  edge_list_type out_edge_ids(vertex_id_type v) const {
    return graph_ptr->out_edge_ids(v);
  }

  //! Get the source vertex of the edge id argument
  vertex_id_type source(edge_id_type edge_id) const {
    return graph_ptr->source(edge_id);
  }

  //! get the target vertex of the edge id argument
  vertex_id_type target(edge_id_type edge_id) const {
    return graph_ptr->target(edge_id);
  }
  

  //! Get the consistency model under which this context was acquired
  consistency_model::model_enum consistency() const { 
    return consistency;
  }
  

  void schedule(const vertex_id_type& vertex, 
                const update_functor_type& update_fun) {
    engine_ptr->schedule(cpuid, vertex, update_fun);
  }

  void schedule_in_neighbors(const vertex_id_type& vertex, 
                             const update_functor_type& update_fun) {
    engine_ptr->schedule_in_neighbors(vertex, update_fun);
  }

  void schedule_out_neighbors(const vertex_id_type& vertex, 
                             const update_functor_type& update_fun) {
    engine_ptr->schedule_out_neighbors(vertex, update_fun);
  }

};


} // end namespace graphlab
#endif


