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
 * The original graph.hpp contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Inteface Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */


#ifndef GRAPHLAB_GRAPH3_HPP
#define GRAPHLAB_GRAPH3_HPP

#include <omp.h>
#include <cmath>
#include <stdio.h>
#include <string>
#include <list>
#include <vector>

#include <fstream>

#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/graph/graph_storage.hpp>
#include <graphlab/macros_def.hpp>



uint mmap_from_file(std::string filename, uint *& array);
uint array_from_file(std::string filename, uint *& array);

namespace graphlab { 
  struct edge_type_impl{
   uint _source; 
   uint _target;
   edge_type_impl(uint source, uint target) : _source(source), _target(target) {}
   edge_type_impl() : _source(-1), _target(-1) { }
   uint source() const { return _source; }
   uint target() const { return _target; }
  };

enum iterator_type {INEDGE, OUTEDGE}; 

  typedef edge_type_impl edge_type;
  //typedef edge_type_impl edge_id_type;

    // Internal iterator on edge_types.
    class edge_iterator : 
      public std::iterator<std::forward_iterator_tag, edge_type> {
    public:
      typedef edge_type reference;
    public:
      // Cosntructors
      edge_iterator () : offset(-1), empty(true) { }
     
      edge_iterator (vertex_id_type _center, size_t _offset, 
                     iterator_type _itype, const uint* _gstore_ptr) : 
        center(_center), offset(_offset), itype(_itype), empty(false), 
        gstore_ptr(_gstore_ptr) { }
      
      edge_iterator (const edge_iterator& it) :
        center(it.center), offset(it.offset), itype(it.itype), 
        empty(it.empty), gstore_ptr(it.gstore_ptr) { }
  
      inline edge_type operator*() const  {
        ASSERT_TRUE(!empty);
        return make_value();
      }

      typedef boost::detail::
      operator_arrow_result<edge_type, edge_type, edge_type*> arrow_type;
      inline arrow_type::type operator->() const {
        return arrow_type::make(make_value());
      }


      inline bool operator==(const edge_iterator& it) const {
        return (empty && it.empty) || 
          (empty == it.empty && itype == it.itype && center == it.center && 
           offset == it.offset);
      }

      inline bool operator!=(const edge_iterator& it) const { 
        return !(*this == it);
      }

      inline edge_iterator& operator++() {
        ASSERT_TRUE(!empty);
        ++offset;
        return *this;
      }

      inline edge_iterator operator++(int) {
        ASSERT_TRUE(!empty);
        const edge_iterator copy(*this);
        operator++();
        return copy;
      }


      inline int operator-(const edge_iterator& it) const {
        ASSERT_TRUE(!empty && itype == it.itype && center == it.center);
        return offset - it.offset;
      }

      inline edge_iterator operator+(size_t i) const {
        edge_iterator retval(center, offset+i, itype, gstore_ptr);
        return retval;
      }

    private:
      // Generate the ret value of the iterator.
      inline edge_type make_value() const {
        edge_type ret;
        if (itype == INEDGE) {
          edge_type rvalue(gstore_ptr[offset], center);
          ret = rvalue;
        } else if (itype == OUTEDGE) {
          edge_type rvalue(center, gstore_ptr[offset]);
          ret = rvalue;
        } else {
          logstream(LOG_FATAL) << "Edge iterator type is invalid." 
                               << std::endl;
        }
        return ret;
      }
    private:
      vertex_id_type center;
      uint offset;
      iterator_type itype;
      bool empty;
      const uint* gstore_ptr;
    }; // end of class edge_iterator.



  struct edge_list{
    uint * start_ptr;
    uint * end_ptr;
    uint _size;
    uint source;
    typedef edge_iterator iterator;
    typedef edge_iterator const_iterator;
    typedef edge_type value_type;
    iterator_type itype;

    edge_list(): start_ptr(NULL), end_ptr(NULL), _size(0), source(-1) { }
    edge_list(uint * _start_ptr, uint * _end_ptr, uint size, uint _source){
      start_ptr = _start_ptr; end_ptr = _end_ptr;
      _size = size; source = _source;
    }
    uint size() const { return _size; }
    edge_type operator[](uint i) const{
      ASSERT_LT(i, _size);
      return edge_type(source, *(start_ptr+i));
    }
    edge_iterator begin() const { return edge_iterator(source, 0, itype, start_ptr); }
    edge_iterator end() const { return edge_iterator(source, 0, itype, end_ptr); }
     bool empty() const { return size() == 0; }


  };




  template<typename VertexData, typename EdgeData>
  class graph3 {

   // typedef graph_storage<VertexData, EdgeData> gstore_type;

  public:

    /// The type of a vertex is a simple size_t
    //typedef graphlab::vertex_id_type vertex_id_type;
    typedef uint vertex_id_type;    

    /// The type of an edge id 
    //typedef size_t edge_id_type;
    
    /// Type for vertex colors 
    typedef char vertex_color_type;

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
 
    /** Represents an edge with source() and target()*/
    /** The type of the edge list */
    //typedef typename gstore_type::edge_list edge_list;
    typedef edge_type_impl edge_type;
    typedef edge_type_impl edge_id_type;
 
    /** Interface for iupdate functor.*/
    //typedef typename gstore_type::edge_list edge_list_type;
    typedef edge_list edge_list_type;

    /** Edge list type for temporary insertion. */
    //typedef typename gstore_type::edge_info edge_info;


    size_t num_nodes, _num_edges;
    uint * node_in_degrees;
    uint * node_out_degrees;
    uint * node_in_edges;
    uint * node_out_edges;
    std::vector<VertexData> node_vdata_array;
    char _color; //not implement yet
    EdgeData _edge;

  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    graph3() : num_nodes(0),_num_edges(0),finalized(false), node_in_edges(NULL), node_in_degrees(NULL), node_out_edges(NULL), node_out_degrees(NULL) {  }

    /**
     * Create a graph with nverts vertices.
     */
    graph3(size_t nverts) : finalized(false) { }

    graph3(const graph<VertexData, EdgeData>& g) { (*this) = g; }

    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
      /*finalized = false;
      gstore.clear();
      vertices.clear();
      edges_tmp.clear();
      vcolors.clear();
      ++changeid;*/
    }

    void clear_reserve() {
      /*clear();
      edges_tmp.clear();
      std::vector<VertexData>().swap(vertices);
      std::vector<vertex_color_type>().swap(vcolors);
      gstore.clear_reserve();*/
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
      finalized = true;
    } // End of finalize
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const {
      return num_nodes;
    } // end of num vertices

    /** \brief Get the number of vertices local to this machine */
    size_t local_vertices() const {
      return num_nodes;
    } 

    /** \brief Get the number of edges */
    size_t num_edges() const {
       return undirected? _num_edges/2 : _num_edges;
    } 

    /** \brief Finds an edge.
        The value of the first element of the pair will be true if an 
        edge from src to target is found and false otherwise. If the 
        edge is found, the edge ID is returned in the second element of the pair. */
    edge_type find(const vertex_id_type source,
                   const vertex_id_type target) const {
       return edge_type(-1,-1); //todo
    } // end of find

    edge_type reverse_edge(const edge_type& edge) const {
      //return gstore.find(edge.target(), edge.source());
      return edge_type(-1,-1);
    }


    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
      if (finalized)
      {
        logstream(LOG_FATAL)
          << "Attempting add vertex"
          << "to a finalized graph." << std::endl;
        ASSERT_MSG(false, "Add vertex to a finalized graph.");
      }
    } // End of add vertex;


    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) {
      ASSERT_GE(num_vertices, num_nodes());
      //TODO
    } // End of resize
    
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_type add_edge(vertex_id_type source, vertex_id_type target, 
                          const EdgeData& edata = EdgeData()) {
      return edge_id_type(source, target); //TODO
    } // End of add edge
        
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      return node_vdata_array[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      return node_vdata_array[v];
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target) {
     ASSERT_TRUE(finalized);
     return _edge;
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, vertex_id_type target) const {
     ASSERT_TRUE(finalized);
     return _edge;
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_type edge) { 
      ASSERT_TRUE(finalized);
      return _edge;
    }
    const EdgeData& edge_data(edge_type edge) const {
      //return 
    }

    size_t num_in_edges(const vertex_id_type v) const {
      return node_in_degrees[v];
    }

    size_t num_out_edges(const vertex_id_type v) const {
      return node_out_degrees[v];
    }

    edge_list_type in_edges(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      if (undirected)
         return out_edges(v);
      return edge_list_type(&node_in_edges[node_in_degrees[v]], &node_in_edges[node_in_degrees[v+1]], node_in_degrees[v],v);  
    }

    edge_list_type out_edges(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      return edge_list_type(&node_out_edges[node_out_degrees[v]], &node_out_edges[node_out_degrees[v+1]], node_out_degrees[v+1]-node_out_degrees[v],v);
    }

    const edge_list_type in_edges(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      if (undirected)
        return out_edges(v);
      return edge_list_type(&node_in_edges[node_in_degrees[v]], &node_in_edges[node_in_degrees[v+1]], node_in_degrees[v],v);  
     }

    const edge_list_type out_edges(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      return edge_list_type(&node_out_edges[node_out_degrees[v]], &node_out_edges[node_out_degrees[v+1]], node_out_degrees[v],v);
     }



   
    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_type vertex) const {
      //ASSERT_LT(vertex, num_nodes());
      //return vcolors[vertex];
      return _color;
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_type vertex) {
      ASSERT_LT(vertex, num_nodes);
      //return vcolors[vertex];
      return _color;
    }

    vertex_color_type get_color(vertex_id_type vid) const{
      return char(0);
    }
    
    void set_color(vertex_id_type vid, vertex_color_type col) {
    }
    
    /** \brief This function constructs a heuristic coloring for the 
        graph and returns the number of colors */
    size_t compute_coloring() {
      return -1;
    } // end of compute coloring


    /**
     * \brief Check that the colors satisfy a valid coloring of the graph.
     * return true is coloring is valid;
     */
    bool valid_coloring() const {
      return true;
    }
    
    /** \brief count the number of times the graph was cleared and rebuilt */
    size_t get_changeid() const {
      return 0;
    }

    size_t estimate_sizeof() const {
      /*size_t vlist_size = sizeof(vertices) + sizeof(VertexData) * vertices.capacity();
      size_t vcolor_size = sizeof(vcolors) + sizeof(vertex_color_type) * vcolors.capacity();
      size_t elist_size = edges_tmp.estimate_sizeof(); 
      size_t store_size = gstore.estimate_sizeof();*/

//      printf("graph3: tmplist size: %u, gstoreage size: %u \n", elist_size, store_size);
       return num_nodes*sizeof(uint)+num_edges*sizeof(int);
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
    } // end of save
   

    /** \brief Load the graph from a file */
    void load(const std::string& filename, bool nodes) {
      if (nodes){
         int rc =array_from_file(filename + ".nodes", node_out_degrees);
	 num_nodes = (rc/4)-1;
	 node_vdata_array.resize(num_nodes);
 	 logstream(LOG_INFO) << "Read " << num_nodes << " nodes" << std::endl;
      }
      else {
         int rc = array_from_file(filename + ".edges", node_out_edges);
         _num_edges = (rc/4)-1;
 	 logstream(LOG_INFO) << "Read " << (undirected? _num_edges/2 : _num_edges) << " edges" << std::endl;
      }
    } // end of load

    /**
     * \brief save the adjacency structure to a text file.
     *
     * Save the adjacency structure as a text file in:
     *    src_Id, dest_Id \n
     *    src_Id, dest_Id \n
     * format.
     */
    void save_adjacency(const std::string& filename) const {
    }



    /**
     * builds a topological_sort of the graph returning it in topsort. 
     * 
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    bool topological_sort(std::vector<vertex_id_type>& topsort) const {
      return true;
    } // end of topological sort

    void set_undirected(){ undirected = true; }

    
  private:    
    /** Internal edge class  */   

 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data */
    //std::vector<VertexData> vertices;


    /** Mark whether the graph is finalized.  Graph finalization is a
        costly procedure but it can also dramatically improve
        performance. */
    bool finalized;
    bool undirected;
    
  }; // End of class graph3

  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const graph3<VertexData, EdgeData>& graph) {
    typedef typename graphlab::graph3<VertexData, EdgeData>::vertex_id_type 
      vertex_id_type;
    typedef typename graphlab::graph3<VertexData, EdgeData>::edge_id_type 
      edge_id_type;
    for(vertex_id_type vid = 0; vid < graph.num_vertices(); ++vid) {
      foreach(edge_id_type eid, graph.out_edge_ids(vid))
        out << vid << ", " << graph.target(eid) << '\n';      
    }
    return out;
  }
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif

