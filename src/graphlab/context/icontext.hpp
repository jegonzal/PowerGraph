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


/** \file 
 *
 * This file describes the icontext interface as well as the the
 * context_range_enum.
 *
 */
#ifndef GRAPHLAB_ICONTEXT_HPP
#define GRAPHLAB_ICONTEXT_HPP

#include <set>
#include <vector>
#include <cassert>

#include <graphlab/contex/iglobal_context.hpp>
#include <graphlab/graph/graph.hpp>

#include <graphlab/context/consistency_model.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {


  inline consistency_model::lock_type_enum 
  central_vertex_lock_type(consistency_model::model_enum srange) {
    switch (srange) {
    case consistency_model::NULL_CONSISTENCY:
      return consistency_model::NO_LOCK;
    case consistency_model::VERTEX_READ_CONSISTENCY:
    case consistency_model::READ_CONSISTENCY:
      return consistency_model::READ_LOCK;
    case consistency_model::VERTEX_CONSISTENCY:
    case consistency_model::EDGE_CONSISTENCY:
    case consistency_model::FULL_CONSISTENCY:
      return consistency_model::WRITE_LOCK;
    case consistency_model::USE_DEFAULT:
      logstream(LOG_FATAL) 
        << "USE_DEFAULT not supported for lock requests!" << std::endl;
    default:
      logstream(LOG_FATAL) << "UNREACHABLE STATE!" << std::endl;
      // unreachable
      return consistency_model::NO_LOCK;
    }
  } // end of central_vertex_lock_type



  inline consistency_model::lock_type_enum 
  adjacent_vertex_lock_type(consistency_model::model_enum srange) {
    switch (srange) {
    case consistency_model::NULL_CONSISTENCY:
    case consistency_model::VERTEX_READ_CONSISTENCY:
    case consistency_model::VERTEX_CONSISTENCY:
      return consistency_model::NO_LOCK;
    case consistency_model::READ_CONSISTENCY:
    case consistency_model::EDGE_CONSISTENCY:
      return consistency_model::READ_LOCK;
    case consistency_model::FULL_CONSISTENCY:
      return consistency_model::WRITE_LOCK;
    case consistency_model::USE_DEFAULT:
      logstream(LOG_FATAL) 
        << "USE_DEFAULT not supported for lock requests!" << std::endl;
    default:
      logstream(LOG_FATAL) << "UNREACHABLE STATE!" << std::endl;
      // unreachable
      return consistency_model::NO_LOCK;
    }
  } // end of adjacent_vertex_lock_type


  inline bool context_is_subset_of(consistency_model::model_enum A,
                                   consistency_model::model_enum B) {
    /*
      if (A==consistency_model::READ_CONSISTENCY && B ==
      consistency_model::VERTEX_CONSISTENCY) return false; else return
      (A < B);
    */  
    return (!(A==consistency_model::READ_CONSISTENCY 
              && B == consistency_model::VERTEX_CONSISTENCY)) 
      && (A < B);
  } // end of context_is_subset_of








  /**
   * \brief represents the data associated with a vertex its adjacent
   * edges and neighbors.
   *
   * The update function is passed an instance of icontext and is uses
   * that instance to obtain information about the graph structure and
   * to read and modify the graph data.
   */
  template<typename Graph, typename UpdateFunctor>
  class icontext : iglobal_context {
  public:
    //! The type of graph that the icontext operates on
    typedef Graph           graph_type;
    typedef UpdateFunctor   update_functor_type;
    typedef typename graph_type::vertex_id_type    vertex_id_type;
    typedef typename graph_type::edge_id_type      edge_id_type;
    typedef typename graph_type::vertex_color_type vertex_color_type;
    //! The edge data type associated with the graph
    typedef typename graph_type::edge_list_type    edge_list_type;
    //! The vertex data type associated with the graph
    typedef typename graph_type::vertex_data_type  vertex_data_type;
    //! The edge data type associated with the graph
    typedef typename graph_type::edge_data_type    edge_data_type;

  public:    
    
    /** icontext destructor */
    virtual ~icontext() { }
    
        
    /**
     * \brief Returns the vertex id of the base vertex in this context.
     *
     * This method is used by the update function to get the base
     * vertex of the context.  The base vertex is the vertex that the
     * update function is being applied to.
     */
    virtual vertex_id_type vertex_id() const = 0;
        
 
    /**
     * \brief Get the vertex color of the base vertx in this context
     */
    virtual vertex_color_type color() const = 0;


    /** 
     * \brief edge lookup from source target pair to edge id.
     *
     * This is used to get structural information from the graph.  If
     * the edge is not present this method will fail. 
     */
    virtual edge_id_type edge(vertex_id_type source,
                              vertex_id_type target) const = 0;

    /**
     * \brief test whether an edge is present
     *
     * This method tests whether the edge exists.  If the edge exists
     * this method returns true.  
     */
    virtual bool edge_exists(vertex_id_type source,
                             vertex_id_type target) const  = 0;
    
    /**
     * \brief Get the reverse edge.
     *
     * Get the reverse edge id.  If no such edge exists this method
     * will fail. 
     */
    virtual edge_id_type reverse_edge(edge_id_type eid) const  = 0;

    /** 
     * \brief get all in edges to the base vertex of this context. 
     * 
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    virtual edge_list_type in_edge_ids() const = 0;

    /** 
     * \brief get all in edge ids to the vertex argument
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    virtual edge_list_type in_edge_ids(vertex_id_type v) const = 0;

    /** 
     * \brief get all out edge ids to the base vertex of this context
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    virtual edge_list_type out_edge_ids() const = 0;

    /** 
     * \brief get all out ede ids to the vertex argument.
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    virtual edge_list_type out_edge_ids(vertex_id_type v) const = 0;

    //! Get the source vertex of the edge id argument
    virtual vertex_id_type source(edge_id_type edge_id) const = 0;

    //! get the target vertex of the edge id argument
    virtual vertex_id_type target(edge_id_type edge_id) const = 0; 
    
    //! Get the consistency model under which this context was acquired
    virtual consistency_model::model_enum consistency() const = 0; 

    /**
     * \brief Get a mutable reference to the data associated with the
     * base vertex
     *
     * Certain optimizations may be made in future version of graphlab
     * that could benefit from immutable vertex_data requests.
     * Therefore if the vertex data does not need to be mutable use
     * the const reference version of vertex_data.
     */
    virtual vertex_data_type& vertex_data() = 0;

    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     * \deprecated use const_vertex_data
     * This should be called if the data does not need to be modified.
     */
    virtual const vertex_data_type& vertex_data() const = 0;
    
    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const vertex_data_type& const_vertex_data() const = 0; 
    
    /**
     * \brief Get a mutable reference to the data associated with the
     * edge.
     *
     * This should only be invoked on edges that are adjacent to the
     * base vertex.  If the edge_data is not going to be modified the
     * const version of this function should be used to permit further
     * optimization.
     */
    virtual edge_data_type& edge_data(edge_id_type eid) = 0;

    /**
     * \brief Get an immutable reference to the data associated with
     * the edge.
     * \deprecated use const_edge_data
     * This should only be invoked on edges that are adjacent to the
     * base vertex. 
     */
    virtual const edge_data_type& edge_data(edge_id_type eid) const = 0;    
    
    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const edge_data_type& const_edge_data(edge_id_type eid) const = 0; 
    
    /**
     * \brief get a mutable reference to the data associated with a
     * neighboring vertex.
     *
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     * If the neighboring vertex data is not going to be modified the
     * const version of this function should be called to permit
     * further optimization by the graphlab engine.
     */
    virtual vertex_data_type& neighbor_vertex_data(vertex_id_type vertex) = 0; 

    /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     * \deprecated Use const_neighbor_vertex_data
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& neighbor_vertex_data(vertex_id_type vertex) const = 0;
        
    /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     *
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& 
    const_neighbor_vertex_data(vertex_id_type vertex) const = 0;


    /**
     * Adds a task to execute the update function on the vertex with
     * the given priority.
     */
    virtual void schedule(const vertex_id_type& vertex, 
                          const update_functor_type& update_fun) = 0;    


    /**
     * Schedule an update on all the neighbors of a particular vertex
     */
    virtual void schedule_in_neighbors(const vertex_id_type& vertex, 
                                       const update_functor_type& update_fun) = 0;

    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    virtual void schedule_out_neighbors(const vertex_id_type& vertex, 
                                        const update_functor_type& update_fun) = 0;
                                                  
 
    // /**
    //    Experimental context upgrade scheme. Returns true if context upgrade is 
    //    successful. If this ever returns false, you are hosed. Should work
    //    with general_context. Note that after context_upgrade is called, any graph
    //    data within the context may change due to a race between releasing and 
    //    reacquiring the upgraded context.
    // */
    // virtual bool 
    // experimental_context_upgrade(consistency_model::model_enum newrange) { 
    //   return false;
    // }
    

  }; // end of icontexty
  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

