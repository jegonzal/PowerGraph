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


/**
 * \file
 *
 * This file contains the template for the graphlab graph
 * data-structure.
 *
 */

#ifndef GRAPHLAB_GRAPH_HPP
#define GRAPHLAB_GRAPH_HPP

#include <omp.h>
#include <cmath>

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
  
  template<typename VertexData, typename EdgeData>
  class graph {
  public:


    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;
    
    /// The type of an edge id 
    typedef graphlab::edge_id_type edge_id_type;
    
    /// Type for vertex colors 
    typedef graphlab::vertex_color_type vertex_color_type;
        
    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
    
    /** This class represents an edge with: src, target and pointer to
        edge_data. */
    class edge_type {
      const graph* graph_ptr;
      const edge_id_type _id;
    public:
      edge_type(const graph* graph_ptr = NULL, 
                const edge_id_type id = -1) : 
        graph_ptr(graph_ptr), _id(id) { }     
      const vertex_id_type& source() const {
        ASSERT_FALSE(empty()); 
        ASSERT_LT(_id, graph_ptr->edges.size());
        return graph_ptr->edges[_id].source(); 
      }
      const vertex_id_type& target() const { 
        ASSERT_FALSE(empty());
        ASSERT_LT(_id, graph_ptr->edges.size());
        return graph_ptr->edges[_id].target(); 
      }
      bool empty() const { return graph_ptr == NULL; }
      friend class graph;
    }; // end of class edge_type.
    friend class edge_type;


    /** This class represents a lazy list of edge_type. */
    class edge_list_type {
    public:
      /**
       * The edge functor takes an edge id and returns an edge object
       */ 
      struct edge_functor : 
        public std::unary_function<edge_id_type, edge_type> {
        const graph* graph_ptr;
        edge_functor(const graph* graph_ptr) : graph_ptr(graph_ptr) { }
        edge_type operator()(const edge_id_type eid) const {
          return edge_type(graph_ptr, eid);
        }
      }; // end of edge_functor
      typedef boost::transform_iterator<edge_functor, const edge_id_type*> iterator;
      typedef iterator const_iterator;
      typedef edge_type value_type;
    private:
      iterator begin_iter, end_iter;
    public:
      edge_list_type(const graph* graph_ptr = NULL, 
                     const edge_id_type* begin_iter = NULL, 
                     const edge_id_type* end_iter = NULL) : 
        begin_iter(begin_iter, edge_functor(graph_ptr)),
        end_iter(end_iter, edge_functor(graph_ptr)) { }
      size_t size() const { return end_iter - begin_iter; }
      edge_type operator[](size_t i) const {
        ASSERT_LT(i, size()); return *(begin_iter + i);
      }
      const_iterator begin() const { return begin_iter; }
      const_iterator end() const { return end_iter; }
      bool empty() const { return size() == 0; }  
    }; // end of class edge_list_type


  private:

    /** Internal edge class  */   
    class edge_info {
      vertex_id_type _source, _target;
      EdgeData _data;
    public:
      edge_info(vertex_id_type source = -1, vertex_id_type target = -1) :
        _source(source), _target(target) { }
      edge_info(vertex_id_type source, vertex_id_type target, EdgeData data) : 
        _source(source), _target(target), _data(data) { }
      bool operator<(const edge_info& other) const {
        return (_source < other._source) || 
          (_source == other._source && _target < other._target); 
      }
      inline const vertex_id_type& source() const { return _source; }
      inline const vertex_id_type& target() const { return _target; }   
      inline EdgeData& data() { return _data; }
      inline const EdgeData& data() const { return _data; }
      void load(iarchive& arc) {
        arc >> _source >> _target >> _data;
      }
      void save(oarchive& arc) const {
        arc << _source << _target << _data;
      }
    }; // end of edge_data
      
      


    
 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data */
    std::vector<VertexData> vertices;

    /** The edge data is a vector of edges where each edge stores its
        source, destination, and data. */
    std::vector<edge_info> edges;
    
    /** A map from src_vertex -> dest_vertex -> edge index */   
    std::vector< std::vector<edge_id_type> >  in_edge_ids;
    
    /** A map from src_vertex -> dest_vertex -> edge index */   
    std::vector< std::vector<edge_id_type> >  out_edge_ids;
    
    /** The vertex colors specified by the user. **/
    std::vector< vertex_color_type > vcolors;  
    
    /** Mark whether the graph is finalized.  Graph finalization is a
        costly procedure but it can also dramatically improve
        performance. */
    bool finalized;
    

     

   
  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    graph() : finalized(true) {  }

    /**
     * Create a graph with nverts vertices.
     */
    graph(size_t nverts) : 
      vertices(nverts), in_edge_ids(nverts), out_edge_ids(nverts), vcolors(nverts),
      finalized(true)  { }


    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
      vertices.clear();
      edges.clear();
      in_edge_ids.clear();
      out_edge_ids.clear();
      vcolors.clear();
      finalized = true;
    }

    /* Free the memory reserved by data members */
    void clear_memory() {
      clear();
      std::vector<VertexData>().swap(vertices);
      std::vector<edge_info>().swap(edges);
      std::vector< std::vector<edge_id_type> >().swap(in_edge_ids);
      std::vector< std::vector<edge_id_type> >().swap(out_edge_ids);
      std::vector<vertex_color_type>().swap(vcolors);
    }



    /**
     * Finalize a graph by sorting its edges to maximize the
     * efficiency of graphlab.  
     * This function takes O(|V|log(degree)) time and will 
     * fail if there are any duplicate edges.
     * This is also automatically invoked by the engine at
     * start.
     */
    void finalize() {   
      // check to see if the graph is already finalized
      if(finalized) return;
      edge_id_less_functor edge_id_less(*this);      
      // Sort all in edges set
#pragma omp parallel for
      for(ssize_t i = 0; i < ssize_t(in_edge_ids.size()); ++i) {
        std::vector<edge_id_type>& eset(in_edge_ids[i]);
        // Sort the edge vector
        std::sort(eset.begin(), eset.end(), edge_id_less);
        // Test for duplicate edges
        if (eset.size() > 1) {
          for(size_t j = 0; j < eset.size()-1; ++j) {
            // Duplicate edge test
            if(!edge_id_less(eset[j], eset[j+1])) {
              logstream(LOG_FATAL)
                << "Duplicate edge "
                << "(" << edges[eset[j]].source() << ", " 
                << edges[eset[j]].target() << ") "
                << "found!  GraphLab does not support graphs "
                << "with duplicate edges." << std::endl;
            }
          }
        }  
      } // end of for loop
      // Sort all out edges sets
#pragma omp parallel for
      for(ssize_t i = 0; i < ssize_t(out_edge_ids.size()); ++i) {
        std::vector<edge_id_type>& eset(out_edge_ids[i]);
        std::sort(eset.begin(), eset.end(), edge_id_less);
        // Test for dupliate edges
        if (eset.size() > 1) {
          for(size_t j = 0; j < eset.size()-1; ++j) {
            // Duplicate edge test
            if(!edge_id_less(eset[j], eset[j+1])) {
              logstream(LOG_FATAL)
                << "Duplicate edge "
                << "(" << edges[eset[j]].source() << ", " 
                << edges[eset[j]].target() << ") "
                << "found!  GraphLab does not support graphs "
                << "with duplicate edges." << std::endl;
            }
          }
        }
      }
      finalized = true;
    } // End of finalize
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return vertices.size(); } 

    /** \brief Get the number of vertices local to this machine */
    size_t local_vertices() const { return vertices.size(); } 

    /** \brief Get the number of edges */
    size_t num_edges() const { return edges.size(); }


    /** \brief Get the number of in edges of a particular vertex */
    size_t num_in_neighbors(vertex_id_type v) const {
      ASSERT_LT(v, vertices.size());
      return in_edge_ids[v].size();
    } // end of num vertices
    
    /** \brief Get the number of out edges of a particular vertex */
    size_t num_out_neighbors(vertex_id_type v) const  {
      ASSERT_LT(v, vertices.size());
      return out_edge_ids[v].size();
    } // end of num vertices

    /** \brief Finds an edge.  
     *
     * The value of the first element of the pair will be true if an
     *  edge from src to target is found and false otherwise. If the
     *  edge is found, the edge ID is returned in the second element
     *  of the pair. 
     */

    edge_type find(const vertex_id_type source, 
                   const vertex_id_type target) const {
      ASSERT_LT(source, in_edge_ids.size());
      ASSERT_LT(target, out_edge_ids.size());
      // Check the base case that the souce or target have no edges
      if (in_edge_ids[target].size() == 0 || out_edge_ids[source].size() == 0) {
        return edge_type();
      } else if(finalized) { // O( log degree ) search ========================>
        // if it is finalized then do the search using a binary search
        // If their are fewer in edges into the target search the in
        // edges
        if(in_edge_ids[target].size() < out_edge_ids[source].size()) {
          // search the source vertices for the edge
          size_t index = binary_search(in_edge_ids[target], source, target);
          if(index < in_edge_ids[target].size())
            return edge_type(this, in_edge_ids[target][index]);
          else return edge_type();
        } else { // If their are fewer edges out of the source binary
                 // search there
          // search the source vertices for the edge
          size_t index = binary_search(out_edge_ids[source], source, target);
          if(index < out_edge_ids[source].size())
            return edge_type(this, out_edge_ids[source][index]);
          else return edge_type();
        }
      } else { // O( degree ) search ==========================================>
        // if there are few in edges at the target search there
        if(in_edge_ids[target].size() < out_edge_ids[source].size()) {
          // linear search the in_edges at the target 
          foreach(const edge_id_type& e, in_edge_ids[target]) {
            if(edges[e].source() == source && edges[e].target() == target) 
              return edge_type(this, e);
          }
          return edge_type();
        } else { // fewer out edges at the source
          // linear search the out_edges at the source
          foreach(const edge_id_type& e, out_edge_ids[source]) {
            if(edges[e].source() == source && edges[e].target() == target) 
              return edge_type(this, e);
          }
          return edge_type();
        }
      } // End of else 
    } // end of find
    

    // edge_type find(const vertex_id_type source, 
    //                const vertex_id_type target) const {
    //   ASSERT_LT(source, in_edge_ids.size());
    //   ASSERT_LT(target, out_edge_ids.size());
    //   // Check the base case that the souce or target have no edges
    //   if (in_edge_ids[target].size() == 0 || out_edge_ids[source].size() == 0)
    //     return edge_type();
    //   // Find the edge
    //   edge_id_less_functor less_functor(*this);
    //   typedef std::vector<edge_id_type>::const_iterator iterator_type;
    //   // if out edges is a smaller set
    //   if(out_edge_ids[source].size() < in_edge_ids[target].size()) {
    //     const iterator_type iter = finalized? 
    //       std::lower_bound(out_edge_ids[source].begin(), out_edge_ids[source].end(), 
    //                        less_functor) :
    //       std::find(out_edge_ids[source].begin(), out_edge_ids[source].end(), 
    //                 less_functor);
    //     if(iter == out_edge_ids[source].end()) return edge_type();
    //     else return edge_type(this, *iter);
    //   } else {
    //     const iterator_type iter = finalized? 
    //       std::lower_bound(in_edge_ids[target].begin(), in_edge_ids[target].end(), 
    //                        less_functor) :
    //       std::find(in_edge_ids[target].begin(), in_edge_ids[target].end(), 
    //                 less_functor);
    //     if(iter == in_edge_ids[target].end()) return edge_type();
    //     else return edge_type(this, *iter);
    //   }
    // } // end of find

    edge_type reverse_edge(const edge_type& edge) const {
      ASSERT_EQ(edge.graph_ptr, this);
      return find(edge.target(), edge.source());
    }

    // vertex_id_type source(const edge_type& edge) const {
    //   ASSERT_EQ(edge.graph_ptr, this);
    //   return edge.source();
    // }

    // vertex_id_type target(const edge_type& edge) const {
    //   ASSERT_EQ(edge.graph_ptr, this);
    //   return edge.target();
    // }


    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
      vertices.push_back(vdata);
      // Resize edge maps
      out_edge_ids.push_back(std::vector<edge_id_type>()); // resize(vertices.size());
      in_edge_ids.push_back(std::vector<edge_id_type>()); // resize(vertices.size());
      vcolors.push_back(vertex_color_type()); // resize(vertices.size());
      return (vertex_id_type)vertices.size() - 1;
    } // End of add vertex;


    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) {
      ASSERT_GE(num_vertices, vertices.size());
      vertices.resize(num_vertices);
      // Resize edge maps
      out_edge_ids.resize(vertices.size());
      in_edge_ids.resize(vertices.size());
      vcolors.resize(vertices.size());
    } // End of resize
    
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData()) {
      if ( source >= vertices.size() 
           || target >= vertices.size() ) {
        logstream(LOG_FATAL) 
          << "Attempting add_edge (" << source
          << " -> " << target
          << ") when there are only " << vertices.size() 
          << " vertices" << std::endl;
        ASSERT_MSG(source < vertices.size(), "Invalid source vertex!");
        ASSERT_MSG(target < vertices.size(), "Invalid target vertex!");
      }
      if(source == target) {
        logstream(LOG_FATAL) 
          << "Attempting to add self edge (" 
          << source << " -> " << target <<  ").  "
          << "This operation is not permitted in GraphLab!" << std::endl;
        ASSERT_MSG(source != target, "Attempting to add self edge!");
      }

      // Add the edge to the set of edge data (this copies the edata)
      edges.push_back( edge_info( source, target, edata ) );

      // Add the edge id to in and out edge maps
      edge_id_type edge_id = (edge_id_type)edges.size() - 1;
      in_edge_ids[target].push_back(edge_id);
      out_edge_ids[source].push_back(edge_id);

      // Determine if the graph is still finalized A graph is
      // finalized if it was finalized and the newly added edge_id is
      // in the correct location in the in and out edge lists (which
      // is true if either the lists only contain a single element or
      // the last two elements are in the correct order).
      edge_id_less_functor edge_id_less(*this);
      finalized = finalized &&
        ((in_edge_ids[target].size() < 2) ||
         edge_id_less(*(in_edge_ids[target].end()-2),
                      *(in_edge_ids[target].end()-1))) &&
        ((out_edge_ids[source].size() < 2) ||
         edge_id_less(*(out_edge_ids[source].end()-2),
                      *(out_edge_ids[source].end()-1)));
    } // End of add edge
        
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(vertex_id_type v) {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_type v) const {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target) {
      ASSERT_LT(source, vertices.size());
      ASSERT_LT(target, vertices.size());
      const edge_type edge = find(source, target);
      // We must find the edge!
      if(edge.empty()) {
        logstream(LOG_FATAL) 
          << "Edge " << source << "-->" << target << " not found!" << std::endl;
      }
      // the edge id should be valid!
      ASSERT_LT(edge._id, edges.size());
      return edges[edge._id].data();
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, 
                              vertex_id_type target) const {
      ASSERT_LT(source, vertices.size());
      ASSERT_LT(target, vertices.size());
      const edge_type edge = find(source, target);
      // We must find the edge!
      if(edge.empty()) {
        logstream(LOG_FATAL) 
          << "Edge " << source << "-->" << target << " not found!" << std::endl;
      }
      // the edge id should be valid!
      ASSERT_LT(edge._id, edges.size());
      return edges[edge._id].data();
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_type edge) {
      ASSERT_EQ(edge.graph_ptr, this);
      ASSERT_LT(edge._id, edges.size());
      return edges[edge._id].data();
    }
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(edge_type edge) const {
      ASSERT_EQ(edge.graph_ptr, this);
      ASSERT_LT(edge._id, edges.size());
      return edges[edge._id].data();
    }

    size_t estimate_sizeof() const {
      const size_t eid_size = sizeof(edge_id_type);
      const size_t vlist_size = sizeof(vertices) + 
        vertices.capacity() *sizeof(VertexData);
      const size_t vcolor_size = sizeof(vcolors) + 
        vcolors.capacity() * sizeof(vertex_color_type);
      const size_t elist_size = sizeof(edges) + 
        edges.capacity() * sizeof(edge_info);
      const size_t inout_shell_size = sizeof(in_edge_ids) + 
        in_edge_ids.capacity() * sizeof(std::vector<edge_id_type>) + 
        sizeof(out_edge_ids) + 
        out_edge_ids.capacity() * sizeof(std::vector<edge_id_type>);
      size_t inout_content_size = 0;
      foreach(std::vector<edge_id_type> eid_list, in_edge_ids) {
        inout_content_size += sizeof(eid_size) * eid_list.capacity();
      }
      inout_content_size  *= 2;
      return vlist_size + vcolor_size + elist_size + inout_shell_size + 
        inout_content_size;
    } // end of estimate sizeof
    
    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_type vertex) const {
      ASSERT_LT(vertex, vertices.size());
      return vcolors[vertex];
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_type vertex) {
      ASSERT_LT(vertex, vertices.size());
      return vcolors[vertex];
    }

    vertex_color_type get_color(vertex_id_type vid) const { 
      return color(vid); 
    }
    
    void set_color(vertex_id_type vid, vertex_color_type col) {
      color(vid) = col;
    }
    
    /** \brief This function constructs a heuristic coloring for the 
        graph and returns the number of colors */
    size_t compute_coloring() {
      // Reset the colors
      for(vertex_id_type v = 0; v < num_vertices(); ++v) color(v) = 0;
      // construct a permuation of the vertices to use in the greedy
      // coloring. \todo Should probably sort by degree instead when
      // constructing greedy coloring.
      std::vector<std::pair<vertex_id_type, vertex_id_type> > 
	permutation(num_vertices());

      for(vertex_id_type v = 0; v < num_vertices(); ++v) 
        permutation[v] = std::make_pair(-num_in_neighbors(v), v);
      //      std::random_shuffle(permutation.begin(), permutation.end());
      std::sort(permutation.begin(), permutation.end());
      // Recolor
      size_t max_color = 0;
      std::set<vertex_color_type> neighbor_colors;
      for(size_t i = 0; i < permutation.size(); ++i) {
        neighbor_colors.clear();
        const vertex_id_type& vid = permutation[i].second;
        // Get the neighbor colors
        foreach(edge_type edge, in_edges(vid)){
          const vertex_id_type& neighbor_vid = edge.source();
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          neighbor_colors.insert(neighbor_color);
        }
        foreach(edge_type edge, out_edges(vid)){
          const vertex_id_type& neighbor_vid = edge.target();
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          neighbor_colors.insert(neighbor_color);
        }

        vertex_color_type& vertex_color = color(vid);
        vertex_color = 0;
        foreach(vertex_color_type neighbor_color, neighbor_colors) {
          if(vertex_color != neighbor_color) break;
          else vertex_color++;
          // Ensure no wrap around
          ASSERT_NE(vertex_color, 0);                
        }
        max_color = std::max(max_color, size_t(vertex_color) );

      }
      // Return the NUMBER of colors
      return max_color + 1;           
    } // end of compute coloring


    /**
     * \brief Check that the colors satisfy a valid coloring of the graph.
     * return true is coloring is valid;
     */
    bool valid_coloring() const {
      for(vertex_id_type vid = 0; vid < num_vertices(); ++vid) {
        const vertex_color_type& vertex_color = color(vid);
        // Get the neighbor colors
        foreach(const edge_type& edge, in_edges(vid)){
          const vertex_id_type& neighbor_vid = edge.source();
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          if(vertex_color == neighbor_color) return false;
        }
      }
      return true;
    }
    
   
    /**
     * Get all the edge which edge.target() == v
     */
    edge_list_type in_edges(const vertex_id_type v) const {
      ASSERT_LT(v, in_edge_ids.size());
      return edge_list_type(this, &(*(in_edge_ids[v].begin())), 
                            &(*in_edge_ids[v].end()));
    }

    /**
     * Get all the edges which edge.source() == v
     */
    edge_list_type out_edges(const vertex_id_type v) const {
      ASSERT_LT(v, out_edge_ids.size());
      return edge_list_type(this, &(*out_edge_ids[v].begin()), 
                            &(*out_edge_ids[v].end()));
    }

   
    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();    
      // read the vertices and colors
      arc >> vertices
          >> edges
          >> in_edge_ids
          >> out_edge_ids
          >> vcolors
          >> finalized;
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << vertices
          << edges
          << in_edge_ids
          << out_edge_ids
          << vcolors
          << finalized;
    } // end of save
    

    /** \brief Load the graph from a file */
    void load(const std::string& filename) {
      std::ifstream fin(filename.c_str());
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();
    } // end of load


    /**
     * \brief save the graph to the file given by the filename
     */    
    void save(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save
    
  private:    

    // PRIVATE HELPERS =========================================================>
    struct edge_id_less_functor {
      const graph& graph_ref;
      edge_id_less_functor(const graph& graph_ref) : graph_ref(graph_ref) { }
      bool operator()(const edge_id_type a, const edge_id_type b) {
        return graph_ref.edges[a] < graph_ref.edges[b];
      }
    }; // end of edge_id_less_functor
    friend class edge_id_less_functor;


    /**
     * This function tries to find the edge in the vector.  If it
     * fails it returns size_t(-1)
     * TODO: switch to stl binary search
     */
    size_t binary_search(const std::vector<edge_id_type>& vec,
                         vertex_id_type source, vertex_id_type target) const {
      // Ensure that the graph is finalized before calling this function
      //      finalize();
      ASSERT_TRUE(finalized);
      // Compare to the middle of the list
      size_t first = 0;
      size_t last = vec.size() - 1;
      while(first <= last) {
        size_t mid = (first+last)/2;
        ASSERT_LT(mid, vec.size());
        vertex_id_type mid_source = edges[vec[mid]].source();
        vertex_id_type mid_target = edges[vec[mid]].target();
        // Edge found
        if(mid_source == source && mid_target == target) {
          return mid;
        } else if (first == last) {
          // Can't search further so we fail
          return -1;
        } 
        // otherwise search further
        if(std::make_pair(source, target) <
           std::make_pair(mid_source, mid_target) ) {
          ASSERT_GT(mid, 0);
          // Search left
          last = mid - 1;
        } else {
          // Search right
          first = mid + 1;
        }
      }
      // We failed to find
      return -1;
    } // end of binary search 
   
    
  }; // End of graph

















  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const graph<VertexData, EdgeData>& graph) {
    typedef typename graphlab::graph<VertexData, EdgeData>::vertex_id_type 
      vertex_id_type;
    typedef typename graphlab::graph<VertexData, EdgeData>::edge_id_type 
      edge_id_type;
    typedef typename graphlab::graph<VertexData, EdgeData>::edge_type
      edge_type;
    for(vertex_id_type vid = 0; vid < graph.num_vertices(); ++vid) {
      foreach(edge_type ewrapper, graph.out_edges(vid))
        out << vid << ", " << ewrapper.target << '\n';      
    }
    return out;
  }

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

