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
 * This file describes the iscope interface as well as the the
 * scope_range_enum.
 *
 */
#ifndef GRAPHLAB_SCOPE_HPP
#define GRAPHLAB_SCOPE_HPP

#include <set>
#include <vector>
#include <cassert>

#include <graphlab/graph/graph.hpp>

#include <graphlab/scope/consistency_model.hpp>

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


  inline bool scope_is_subset_of(consistency_model::model_enum A,
                                 consistency_model::model_enum B) {
    /*
      if (A==consistency_model::READ_CONSISTENCY && B ==
      consistency_model::VERTEX_CONSISTENCY) return false; else return
      (A < B);
    */  
    return (!(A==consistency_model::READ_CONSISTENCY 
              && B == consistency_model::VERTEX_CONSISTENCY)) 
      && (A < B);
  } // end of scope_is_subset_of




  /**
   * \brief represents the data associated with a vertex its adjacent
   * edges and neighbors.
   *
   * The update function is passed an instance of iscope and is uses
   * that instance to obtain information about the graph structure and
   * to read and modify the graph data.
   */
  template<typename Graph>
  class iscope {
  public:
    //! The type of graph that the iscope operates on
    typedef Graph graph_type;
    typedef typename graph_type::vertex_id_type    vertex_id_type;
    typedef typename graph_type::edge_id_type      edge_id_type;
    typedef typename graph_type::vertex_color_type vertex_color_type;
    //! The edge data type associated with the graph
    typedef typename graph_type::edge_list_type   edge_list_type;


    //! The vertex data type associated with the graph
    typedef typename graph_type::vertex_data_type vertex_data_type;

    //! The edge data type associated with the graph
    typedef typename graph_type::edge_data_type   edge_data_type;




  public:    

    /** 
     * \brief construct an iscope from a graph This is called by the
     * engine when creating an iscope to be passed into an update
     * function.
     */
    iscope(Graph* graph_ptr = NULL, vertex_id_type vertex = -1,
           consistency_model::model_enum consistency = 
           consistency_model::EDGE_CONSISTENCY) : 
      _graph_ptr(graph_ptr), _vertex(vertex), 
      _consistency(consistency) { }
    
    /** iscope destructor */
    virtual ~iscope() { }
    
    /**
     * \brief commits all changes.
     * This is called by the engine after the update function returns.
     */
    virtual void commit() { }

    /**
     * Get the number of vertices in the graph
     */
    size_t num_vertices() const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->num_vertices();  
    }

    vertex_color_type color() const {
      return _graph_ptr->get_color(_vertex);
    }
    
    /**
     * \brief Returns the vertex id of the base vertex in this scope.
     *
     * This method is used by the update function to get the base
     * vertex of the scope.  The base vertex is the vertex that the
     * update function is being applied to.
     */
    vertex_id_type vertex() const { return _vertex; }

    /** 
     * \brief edge lookup from source target pair to edge id.
     *
     * This is used to get structural information from the graph.  If
     * the edge is not present this method will fail. 
     */
    edge_id_type edge(vertex_id_type source,
                      vertex_id_type target) const {
      assert(_graph_ptr != NULL);
      // No cheating
      // assert(source == _vertex || target == _vertex);
      return _graph_ptr->edge_id(source, target);
    }

    /**
     * \brief test whether an edge is present
     *
     * This method tests whether the edge exists.  If the edge exists
     * this method returns true.  
     */
    bool edge_exists(vertex_id_type source,
                     vertex_id_type target) const {
      assert(_graph_ptr != NULL);
      // No cheating
      // assert(source == _vertex || target == _vertex);
      return _graph_ptr->find(source, target).first;
    }
    

    /**
     * \brief Get the reverse edge.
     *
     * Get the reverse edge id.  If no such edge exists this method
     * will fail. 
     */
    edge_id_type reverse_edge(edge_id_type eid) const {      
      assert(_graph_ptr != NULL);      
      //       // No cheating
      //       assert(_graph_ptr->source(eid) == _vertex ||
      //              _graph_ptr->target(eid) == _vertex);
      return _graph_ptr->rev_edge_id(eid);
    }

    /** 
     * \brief get all in edges to the base vertex of this scope. 
     * 
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type in_edge_ids() const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->in_edge_ids(_vertex);
    }

    /** 
     * \brief get all in edge ids to the vertex argument
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type in_edge_ids(vertex_id_type v) const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->in_edge_ids(v);
    }


    /** 
     * \brief get all out edge ids to the base vertex of this scope
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type out_edge_ids() const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->out_edge_ids(_vertex);
    }

    /** 
     * \brief get all out ede ids to the vertex argument.
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    edge_list_type out_edge_ids(vertex_id_type v) const {
      assert(_graph_ptr != NULL);
      return _graph_ptr->out_edge_ids(v);
    }

    //! Get the source vertex of the edge id argument
    vertex_id_type source(edge_id_type edge_id) const {
      //      assert(_graph_ptr != NULL);
      return _graph_ptr->source(edge_id);
    }

    //! get the target vertex of the edge id argument
    vertex_id_type target(edge_id_type edge_id) const {
      //       assert(_graph_ptr != NULL);
      return _graph_ptr->target(edge_id);
    }
    

    //! Get the consistency model under which this scope was acquired
    consistency_model::model_enum consistency() const { 
      return _consistency;
    }

    /**
     * \brief Get a mutable reference to the data associated with the
     * base vertex
     *
     * Certain optimizations may be made in future version of graphlab
     * that could benefit from immutable vertex_data requests.
     * Therefore if the vertex data does not need to be mutable use
     * the const reference version of vertex_data.
     */
    virtual vertex_data_type& vertex_data() {
      return _graph_ptr->vertex_data(_vertex);
    }

    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     * \deprecated use const_vertex_data
     * This should be called if the data does not need to be modified.
     */
    virtual const vertex_data_type& vertex_data() const {
      return const_vertex_data();
    }
    
    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const vertex_data_type& const_vertex_data() const {
      return _graph_ptr->vertex_data(_vertex);
    }
    
    
    /**
     * \brief Get a mutable reference to the data associated with the
     * edge.
     *
     * This should only be invoked on edges that are adjacent to the
     * base vertex.  If the edge_data is not going to be modified the
     * const version of this function should be used to permit further
     * optimization.
     */
    virtual edge_data_type& edge_data(edge_id_type eid) {
      return _graph_ptr->edge_data(eid);
    }

    /**
     * \brief Get an immutable reference to the data associated with
     * the edge.
     * \deprecated use const_edge_data
     * This should only be invoked on edges that are adjacent to the
     * base vertex. 
     */
    virtual const edge_data_type& edge_data(edge_id_type eid) const {
      return const_edge_data(eid);
    }
    
    
    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const edge_data_type& const_edge_data(edge_id_type eid) const {
      return _graph_ptr->edge_data(eid);
    }
    
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
    virtual vertex_data_type& neighbor_vertex_data(vertex_id_type vertex) {
      return _graph_ptr->vertex_data(vertex);
    }

    /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     * \deprecated Use const_neighbor_vertex_data
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& 
    neighbor_vertex_data(vertex_id_type vertex) const {
      return const_neighbor_vertex_data(vertex);
    }
        

    /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     *
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& 
    const_neighbor_vertex_data(vertex_id_type vertex) const {
      return _graph_ptr->vertex_data(vertex);
    }


    /**
       Experimental scope upgrade scheme. Returns true if scope upgrade is 
       successful. If this ever returns false, you are hosed. Should work
       with general_scope. Note that after scope_upgrade is called, any graph
       data within the scope may change due to a race between releasing and 
       reacquiring the upgraded scope.
    */
    virtual bool 
    experimental_scope_upgrade(consistency_model::model_enum newrange) { 
      return false;
    }

  protected:
    /** A pointer to the underlying graph datastructure */
    Graph* _graph_ptr;

    /** The vertex that this graph represents*/
    vertex_id_type _vertex;

    /** The consistency model that this scope ensures */
    consistency_model::model_enum _consistency;
    

  }; // end of iscope
  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

