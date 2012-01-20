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
 * \file
 *
 * This file describes the graph interface that all graphlab graphs
 * must satisfy.  
 *
 */

#ifndef GRAPHLAB_IGRAPH_HPP
#define GRAPHLAB_IGRAPH_HPP

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>


#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>
#include <functional>


#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional.hpp>
#include <boost/iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>



#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/graph/graph_basic_types.hpp>


#include <graphlab/util/random.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab { 


  // CLASS GRAPH ==============================================================>
  /**
     \brief The GraphLab primary Graph container templatized over the
     vertex and edge types.
    
     Every vertex and edge in the graph is assigned a unique integer
     ID.  The type of the vertex id is
     <code>graphlab::vertex_id_type</code> and the type of the edge id
     is <code>graphlab::edge_id_type</code>. Both
     <code>vertex_id_type</code> and <code>edge_id_type</code> are
     currently defined as <code>uint32_t</code>.  While this limits
     the graphs to 4 billion vertices it also helps reduce the storage
     overhead. We encourage users to use the
     <code>vertex_id_type</code> and <code>edge_id_type</code> types
     as they may change in larger distributed systems.

     <h2> Graph Creation </h2>
  
     Vertices and edges are added using the graph::add_vertex()
     and graph::add_edge() member functions:
  
  
     \code
     vertex_id_type graph::add_vertex(const VertexData& vdata = VertexData()) 
     edge_id_type graph::add_edge(vertex_id_type source, vertex_id_type target, 
     const EdgeData& edata = EdgeData()) 
     \endcode
  
     The functions return the id's of the added vertex and edge
     respectively.  An edge can only be added if both the source and
     target vertex id's are already in the graph. Duplicate edges are not 
     supported and may result in undefined behavior.

     The graph behaves like an STL container by storing a local copy of
     any vertex or edge data.  This data can be accessed using
     the useful graph routines described below.

     The Graph datastructure is stored internally as a sorted adjacency
     list.  Where each vertex contains two vectors listing all of its
     in-edges and all of its out-edges. The in-edges are sorted by the
     id of the source vertex, while the out-edges are sorted by the id
     of the destination vertex. This allows for <i>O(log(n))</i> time
     look up of an edge id given its source and destination vertex.
   
     However, this invariant is very expensive to maintain during graph
     consruction.  Therefore, the Graph datastructure allows the
     invariant to be violated during graph::add_vertex() or
     graph::add_edge(). A final call to the member function
     graph::finalize() is needed after graph construction to restore
     the invariant.  The engine routines will defensively call
     graph::finalize() it is not first called by the user.
  */
  
  template<typename Graph>
  class igraph {
  public:
    /// The type of a vertex is a simple size_t
    typedef typename Graph::vertex_id_type vertex_id_type;
    /// The type of an edge id 
    typedef typename Graph::edge_id_type edge_id_type;
    /// Type for vertex colors 
    typedef typename Graph::vertex_color_type vertex_color_type;        
    /** The type of the vertex data stored in the graph */
    typedef typename Graph::vertex_data_type vertex_data_type;
    /** The type of the edge data stored in the graph */
    typedef typename Graph::edge_data_type  edge_data_type;

    /** This class represents an edge with source() and target()*/
    class edge_type : Graph::edge_type {
      typename Graph::edge_type* edge_ptr;
    public:
      inline vertex_id_type source() const { return edge_ptr->source(); }
      inline vertex_id_type target() const { return edge_ptr->target(); }
      inline bool empty() const { return edge_ptr->empty(); }
    };

    /** This class represents a lazy list of edge_type. */
    class edge_list_type : Graph::edge_list_type {
      typedef typename Graph::edge_list_type base_type;
      typedef typename base_type::iterator iterator;
      typedef typename base_type::const_iterator const_iterator;
      typedef typename base_type::value_type value_type;
      base_type elist_ptr;
    public:
      size_t size() const { return elist_ptr->size(); }
      edge_type operator[](size_t i) const {  return (*elist_ptr)[i]; }
      const_iterator begin() const { return elist_ptr->begin(); }
      const_iterator end() const { return elist_ptr->end(); }
      bool empty() const { return elist_ptr->empty(); }  
    }; // end of class edge_list_type

  private:
    Graph* base_graph_ptr;

  public:

    // METHODS =================================================================>
    /**
     * \brief Resets the graph state.
     */
    void clear() { base_graph_ptr->clear(); };

    /**
     * Finalize a graph by sorting its edges to maximize the
     * efficiency of graphlab.  
     * This function takes O(|V|log(degree)) time and will 
     * fail if there are any duplicate edges.
     * This is also automatically invoked by the engine at
     * start.
     */
    void finalize() { base_graph_ptr->finalize(); }
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return base_graph_ptr->num_vertices(); }

    /** \brief Get the number of vertices local to this machine */
    size_t local_vertices() const { return base_graph_ptr->local_vertices(); } 

    /** \brief Get the number of edges */
    size_t num_edges() const { return base_graph_ptr->num_edges(); }

    /** \brief Finds an edge.  
     *
     * The value of the first element of the pair will be true if an
     *  edge from src to target is found and false otherwise. If the
     *  edge is found, the edge ID is returned in the second element
     *  of the pair. 
     */
    edge_type find(const vertex_id_type source, 
                   const vertex_id_type target) const { 
      base_graph_ptr->find(source,target); }


    //! Return the edge in the opposite direction
    edge_type reverse_edge(const edge_type& edge) const { 
      base_graph_ptr->reverse_edge(edge); }

    /** 
     * \brief Creates a vertex containing the vertex data and returns
     * the id of the new vertex id. Vertex ids are assigned in
     * increasing order with the first vertex having id 0.
     */
    vertex_id_type add_vertex(const vertex_data_type& vdata = 
                              vertex_data_type() ) {
      base_graph_ptr->add_vertex(vdata); 
    }

    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) { base_graph_ptr->resize(num_vertices); }
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const edge_data_type& edata = edge_data_type()) {
      base_graph_ptr->add_edge(source, target, edata);
    }


    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_type vertex) {
      return base_graph_ptr->color(vertex);
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_type vid) const {
      return base_graph_ptr->color(vid);
    }
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    vertex_data_type& vertex_data(vertex_id_type vid) {
      return base_graph_ptr->vertex_data(vid);
    }
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const vertex_data_type& vertex_data(vertex_id_type vid) const {
      return base_graph_ptr->vertex_data(vid);
    }

    /** \brief Returns a reference to the data stored on the edge source->target. */
    edge_data_type& edge_data(vertex_id_type source, vertex_id_type target) {
      return base_graph_ptr->edge_data(source, target);
    }
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const edge_data_type& edge_data(vertex_id_type source, 
                                    vertex_id_type target) const {
      return base_graph_ptr->edge_data(source, target);
    }

    /** \brief Returns a reference to the data stored on the edge e */
    edge_data_type& edge_data(edge_type edge) {
      return base_graph_ptr->edge_data(edge);
    }
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const edge_data_type& edge_data(edge_type edge) const {
      return base_graph_ptr->edge_data(edge);
    }

    //! Get all the edge which edge.target() == vid
    edge_list_type in_edges(const vertex_id_type vid) const {
      return base_graph_ptr->in_edges(vid);
    }

    //! Get the number of edges which edge.target() == vid
    size_t num_in_edges(const vertex_id_type vid) const {
      return base_graph_ptr->num_in_edges(vid);
    }

    //! Get all the edges which edge.source() == vid
    edge_list_type out_edges(const vertex_id_type vid) const {
      return base_graph_ptr->out_edges(vid);
    }

    //! Get the number of edges which edge.source() == vid
    size_t num_out_edges(const vertex_id_type vid) const {
      return base_graph_ptr->num_out_edges(vid);
    }

   
    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {  base_graph_ptr->load(arc); }

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const { base_graph_ptr->save(arc); }

    /** \brief Load the graph from a file */
    void load(const std::string& filename) { 
      base_graph_ptr->load(filename); }

    /** \brief save the graph to the file given by the filename */    
    void save(const std::string& filename) const {
      base_graph_ptr->save(filename); }
    
  }; // End of graph







} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

