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

/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#ifndef GRAPHLAB_CONTEXT_HPP
#define GRAPHLAB_CONTEXT_HPP

//#include <boost/bind.hpp>

#include <graphlab/context/iglobal_context.hpp>
#include <graphlab/context/icontext.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab {



  /**
   * This defines a general scope type 
   */
  template<typename Engine>
  class context : 
    public icontext<typename Engine::graph_type, 
                    typename Engine::update_functor_type> {
  public:
   

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
    typedef typename graph_type::edge_wrapper_type   edge_wrapper_type;
    
  private:    
    /** a pointer to the engine */
    engine_type* engine_ptr;
    /** A pointer to the underlying graph datastructure */
    graph_type* graph_ptr;   
    /** A pointer to the scheduler */
    ischeduler_type* scheduler_ptr;
    /** the cpuid of this context */
    size_t cpuid;

    
    /** The vertex that this graph represents*/
    vertex_id_type vid;

    /** The consistency model that this context ensures */
    consistency_model::model_enum _consistency;
 
  public:

    // Cache related members --------------------------------------------------
    // vertex data cache
    struct cache_entry {
      vertex_data_type old;
      vertex_data_type current;
      uint16_t reads, writes;
      cache_entry() : reads(0), writes(0) { }
    };
    typedef std::map<vertex_id_type, cache_entry> cache_map_type;  
    cache_map_type cache;



  public:
    
    /** 
     * \brief construct an icontext from a graph This is called by the
     * engine when creating an icontext to be passed into an update
     * function.
     */
    context(engine_type* engine_ptr = NULL,
            graph_type* graph_ptr = NULL,
            ischeduler_type* scheduler_ptr = NULL,
            size_t cpuid = -1) :
      engine_ptr(engine_ptr), graph_ptr(graph_ptr), 
      scheduler_ptr(scheduler_ptr), cpuid(cpuid),
      vid(-1), _consistency(consistency_model::EDGE_CONSISTENCY) { }
    

    void init(const vertex_id_type vertex, 
              consistency_model::model_enum consistency) {             
      vid = vertex;
      _consistency = consistency;
    } // end of init_vertex

    size_t num_vertices() const { return graph_ptr->num_vertices(); }
    size_t num_edges() const { return graph_ptr->num_edges(); }
    size_t num_updates() const { return engine_ptr->last_update_count(); } 
    void terminate() { engine_ptr->stop(); }

   

    vertex_data_type& vertex_data(const vertex_id_type vid) {
      typedef typename cache_map_type::iterator iterator_type;
      iterator_type iter = cache.find(vid);
      if(iter != cache.end()) {
        return iter->second.current;
      } else {
        return graph_ptr->vertex_data(vid);
      }
    }

    const vertex_data_type& vertex_data(const vertex_id_type vid) const {
      typedef typename cache_map_type::const_iterator iterator_type;
      iterator_type iter = cache.find(vid);
      if(iter != cache.end()) {
        return iter->second.current;
      } else {
        return graph_ptr->vertex_data(vid);
      }
    }



    vertex_data_type& vertex_data() {
      return vertex_data(vid);
    }

    const vertex_data_type& vertex_data() const {
      return vertex_data(vid);
    }

    const vertex_data_type& const_vertex_data() const {
      return vertex_data(vid);
    }

    const vertex_data_type& const_vertex_data(vertex_id_type vid) const {
      return vertex_data(vid);
    }

    // edge_data_type& edge_data(edge_id_type eid) {
    //   return graph_ptr->edge_data(eid);  
    // }

    // const edge_data_type& edge_data(edge_id_type eid) const {
    //   return graph_ptr->edge_data(eid);  
    // }

    // const edge_data_type& const_edge_data(edge_id_type eid) const {
    //   return graph_ptr->edge_data(eid);  
    // }

    edge_data_type& edge_data(vertex_id_type source,
        vertex_id_type target) {
      return graph_ptr->edge_data(source, target);
    }

    const edge_data_type& edge_data(vertex_id_type source,
        vertex_id_type target) const {
      return graph_ptr->edge_data(source, target);
    }


    void commit() { }

    vertex_color_type color() const { return graph_ptr->get_color(vid); }
    
    vertex_id_type vertex_id() const { return vid; }

    
    bool edge_exists(vertex_id_type source,
                     vertex_id_type target) const {
      return graph_ptr->find(source, target).first;
    }

    /** Interfaces about edge_ids are deprecated. */
    // edge_id_type edge(vertex_id_type source,
    //                   vertex_id_type target) const {
    //   return graph_ptr->edge_id(source, target);
    // }


    // edge_id_type reverse_edge(edge_id_type eid) const {      
    //   return graph_ptr->rev_edge_id(eid);
    // }

    // edge_list_type in_edge_ids() const {
    //   return graph_ptr->in_edge_ids(vid);
    // }

    // edge_list_type in_edge_ids(vertex_id_type v) const {
    //   return graph_ptr->in_edge_ids(v);
    // }

    // edge_list_type out_edge_ids() const {
    //   return graph_ptr->out_edge_ids(vid);
    // }

    // edge_list_type out_edge_ids(vertex_id_type v) const {
    //   return graph_ptr->out_edge_ids(v);
    // }

    // std::vector<vertex_id_type> in_vertices() const {
    //   return graph_ptr->in_vertices(vid);
    // }

    // vertex_list_type in_vertices_list() const {
    //   return graph_ptr->in_vertices_list(vid);
    // }

    // std::vector<vertex_id_type> out_vertices() const {
    //   return graph_ptr->out_vertices(vid);
    // }

    // vertex_list_type out_vertices_list() const {
    //   return graph_ptr->out_vertices_list(vid);
    // }

    // //! Get the source vertex of the edge id argument
    // vertex_id_type source(edge_id_type edge_id) const {
    //   return graph_ptr->source(edge_id);
    // }

    // //! get the target vertex of the edge id argument
    // vertex_id_type target(edge_id_type edge_id) const {
    //   return graph_ptr->target(edge_id);
    // }
    
    edge_list_type get_in_edges (vertex_id_type v) {
      return graph_ptr->get_in_edges(v);
    }

    edge_list_type get_out_edges (vertex_id_type v) {
      return graph_ptr->get_out_edges(v);
    }

    edge_list_type get_in_edges () {
      return graph_ptr->get_in_edges(vid);
    }

    edge_list_type get_out_edges () {
      return graph_ptr->get_out_edges(vid);
    }


    //! Get the consistency model under which this context was acquired
    consistency_model::model_enum consistency() const { 
      return _consistency;
    }
    

    void schedule(const vertex_id_type& vertex, 
                  const update_functor_type& update_fun) {
      scheduler_ptr->schedule(cpuid, vertex, update_fun);
    }

    void schedule_in_neighbors(const vertex_id_type& vertex, 
                               const update_functor_type& update_fun) {
      const edge_list_type edges = graph_ptr->get_in_edges(vertex);
      foreach(const edge_wrapper_type& ewrapper, edges) 
        schedule(ewrapper.src, update_fun);
    }

    void schedule_out_neighbors(const vertex_id_type& vertex, 
                               const update_functor_type& update_fun) {
      const edge_list_type edges = graph_ptr->get_out_edges(vertex);
      foreach(const edge_wrapper_type& ewrapper, edges) 
        schedule(ewrapper.target, update_fun);
    }

    void schedule_neighbors(const vertex_id_type& vertex, 
                            const update_functor_type& update_fun) {
      schedule_in_neighbors(vertex, update_fun);
      schedule_out_neighbors(vertex, update_fun);
    }



  protected:

    void acquire_lock(const std::string& key, size_t index = 0) { 
      engine_ptr->acquire_global_lock(key, index); 
    }

    void release_lock(const std::string& key, size_t index = 0) { 
      engine_ptr->release_global_lock(key, index); 
    }

    void commit_change(const std::string& key, size_t index = 0) { }

    void get_global(const std::string& key,                                       
                    graphlab::any_vector*& ret_vec_ptr,
                    bool& ret_is_const) {
      engine_ptr->get_global(key, ret_vec_ptr, ret_is_const);
    } // end of get_any_pair    


  }; // end of context

} // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif

