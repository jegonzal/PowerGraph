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

#include <graphlab/context/iglobal_context.hpp>
#include <graphlab/context/consistency_model.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {



  /**
   * \brief represents the data associated with a vertex its adjacent
   * edges and neighbors.
   *
   * The update function is passed an instance of icontext and is uses
   * that instance to obtain information about the graph structure and
   * to read and modify the graph data.
   */
  template<typename Graph, typename UpdateFunctor>
  class icontext : public iglobal_context {
  public:
    //! The type of graph that the icontext operates on
    typedef Graph           graph_type;
    typedef UpdateFunctor   update_functor_type;
    typedef typename graph_type::vertex_id_type    vertex_id_type;
    typedef typename graph_type::vertex_color_type vertex_color_type;
    //! The edge data type associated with the graph
    typedef typename graph_type::edge_list_type    edge_list_type;
    typedef typename graph_type::edge_type edge_type;

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

    
    //! Reverse an edge
    virtual edge_type reverse_edge(const edge_type& edge) const = 0;

    /**
     * \brief test whether an edge is present
     *
     * This method tests whether the edge exists.  If the edge exists
     * this method returns true.  
     */
    virtual edge_type find(vertex_id_type source,
                           vertex_id_type target) const  = 0;



    /**
     * Get the in edges associated with the center vertex
     */
    virtual edge_list_type in_edges() const = 0;

    /**
     * Get the number of in edges associated with the center vertex 
     */
    virtual size_t num_in_edges() const = 0;

    /**
     * Get the in edges associated with the vertex v
     */
    virtual edge_list_type in_edges(vertex_id_type v) const = 0;

    /**
     * Get the number of in edges associated with the vertex v
     */
    virtual size_t num_in_edges(vertex_id_type v) const = 0;


    /**
     * Get the out edges associated with the center vertex
     */
    virtual edge_list_type out_edges() const = 0;

    /**
     * Get the number of out edges associated with the center vertex
     */
    virtual size_t num_out_edges() const = 0;

    /**
     * Get the out edges associated with the vertex v
     */
    virtual edge_list_type out_edges(vertex_id_type v) const = 0;

    /**
     * Get the number out edges associated with the vertex v
     */
    virtual size_t num_out_edges(vertex_id_type v) const = 0;


    
   
    //! Get the consistency model under which this context was acquired
    virtual consistency_model consistency() const = 0; 

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
    virtual edge_data_type& edge_data(const edge_type& edge) = 0;
    virtual edge_data_type& 
    edge_data(vertex_id_type source, vertex_id_type target) = 0;

    /**
     * \brief Get an immutable reference to the data associated with
     * the vase vertex.
     *
     * This should be called if the data does not need to be modified.
     */    
    virtual const edge_data_type& 
    const_edge_data(const edge_type& eid) const = 0; 
    virtual const edge_data_type& 
    const_edge_data(vertex_id_type source, vertex_id_type target) const = 0;
    
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
    virtual vertex_data_type& vertex_data(vertex_id_type vertex) = 0; 

    /**
     * \brief get an immutable reference to the data associated with a
     * neighboring vertex.
     *
     * This function should only be invoked on neighboring
     * vertices. Unfortunately, due to the Log(d) lookup required to
     * enforce the adjacency constraint we do not check at this time.
     */
    virtual const vertex_data_type& 
    const_vertex_data(vertex_id_type vertex) const = 0;


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
                                                  

    /**
     * Schedule an update on all the out neighbors of a particular vertex
     */
    virtual void schedule_neighbors(const vertex_id_type& vertex, 
                                    const update_functor_type& update_fun) = 0;
                                                  

    

  }; // end of icontexty
  
} // end of namespace
#include <graphlab/macros_undef.hpp>

#endif

