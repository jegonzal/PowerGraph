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

#ifndef GRAPHLAB_GPU_SCOPE_HPP
#define GRAPHLAB_GPU_SCOPE_HPP

#include <graphlab/graph/graph.hpp>



namespace graphlab {
  /**
   * This defines a general scope type 
   */
  template<typename GPUProgram>
  class gpu_scope {
  public:
    typedef typename GPUProgram::vertex_data_type vertex_data_type;
    typedef typename GPUProgram::edge_data_type edge_data_type;

    
  public:

    __device__ gpu_scope() { }
    

  
    /**
     * \brief Returns the vertex id of the base vertex in this scope.
     *
     * This method is used by the update function to get the base
     * vertex of the scope.  The base vertex is the vertex that the
     * update function is being applied to.
     */
    __device__ vertex_id_t vertex() { }


    /** Get the number of vertices */
    __device__ size_t num_vertices() { }


    /** 
     * \brief edge lookup from source target pair to edge id.
     *
     * This is used to get structural information from the graph.  If
     * the edge is not present this method will fail. 
     */
    __device__ edge_id_t edge(vertex_id_t source,
                              vertex_id_t target) const {  }

    /**
     * \brief test whether an edge is present
     *
     * This method tests whether the edge exists.  If the edge exists
     * this method returns true.  
     */
    __device__ bool edge_exists(vertex_id_t source,
                                vertex_id_t target) const {  }
    

    /**
     * \brief Get the reverse edge.
     *
     * Get the reverse edge id.  If no such edge exists this method
     * will fail. 
     */
    __device__ edge_id_t reverse_edge(edge_id_t eid) const {  }

    /** 
     * \brief get all in edges to the base vertex of this scope. 
     * 
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    __device__ edge_list in_edge_ids() const {  }

    /** 
     * \brief get all in edge ids to the vertex argument
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    __device_ edge_list in_edge_ids(vertex_id_t v) const { }


    /** 
     * \brief get all out edge ids to the base vertex of this scope
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    __device__ edge_list out_edge_ids() const {  }

    /** 
     * \brief get all out ede ids to the vertex argument.
     *
     * This method returns an immutable vector of edge ids sorted in
     * order of <source id, dest id> pairs.
     */
    __device__ edge_list out_edge_ids(vertex_id_t v) const { }

    //! Get the source vertex of the edge id argument
    __device__ vertex_id_t source(edge_id_t edge_id) const { }

    //! get the target vertex of the edge id argument
    __device__ vertex_id_t target(edge_id_t edge_id) const { }





    /** Get a reference to the vertex data */
    __device__ vertex_data_type& vertex_data() { }
    __device__ const vertex_data_type& const_vertex_data() const { }

    /// Direct calls to access edge data
    __device__
    edge_data_type& edge_data(edge_id_t eid) { }
    __device__
    const edge_data_type& const_edge_data(edge_id_t eid) const { }
    

    __device__
    const vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) const {
    }

    __device__
    const vertex_data_type& const_neighbor_vertex_data(vertex_id_t vertex) const {
    }

    __device__
    vertex_data_type& neighbor_vertex_data(vertex_id_t vertex) {
    }
    

  }; // end of ext_locked_scope

} // end of graphlab namespace


#endif
