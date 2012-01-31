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
 * Implementation by Danny Bickson, CMU
 *
 */


#ifndef GRAPHLAB_GRAPH3_HPP
#define GRAPHLAB_GRAPH3_HPP
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <omp.h>
#include <cmath>
#include <stdio.h>
#include <string>
#include <vector>

#include <fstream>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/graph/graph_storage.hpp>
#include <graphlab/macros_def.hpp>



template<typename T>
uint array_from_file(std::string filename, T *& array){
          struct stat sb;
          int fd = open (filename.c_str(), O_RDONLY);
          if (fd == -1) {
                  perror ("open");
                  logstream(LOG_FATAL) << "Failed to open input file: " << filename << std::endl;
          }

          if (fstat (fd, &sb) == -1) {
                  perror ("fstat");
                  logstream(LOG_FATAL) << "Failed to get size of  input file: " << filename << std::endl;
          }

          if (!S_ISREG (sb.st_mode)) {
                  logstream(LOG_FATAL) << "Input file: " << filename 
              << " is not a regular file and can not be mapped " << std::endl;
          }
	  close(fd);
 
	  int toread = sb.st_size/sizeof(T); 
          array = new T[toread];
          int total = 0;
	  FILE * f = fopen(filename.c_str(), "r");
          if (f == NULL){
	     perror("fopen");
             logstream(LOG_FATAL) << "Failed to open input file: " << filename << std::endl;
          }
         
          while(total < toread){
	     int rc = fread(array+total, sizeof(T), toread-total,f);
	     if (rc < 0 ){
	       perror("fread");
               logstream(LOG_FATAL) << "Failed to read from input file: " << filename << std::endl;
	     }
	     total += rc; 
          }
          return sb.st_size;
}

enum iterator_type {INEDGE, OUTEDGE}; 


uint mmap_from_file(std::string filename, uint *& array);

namespace graphlab { 
  struct edge_type_impl{
   uint _source; 
   uint _target;
   uint _offset;
   edge_type_impl(uint source, uint target, uint offset) : _source(source), _target(target), _offset(offset) {}
   edge_type_impl() : _source(-1), _target(-1), _offset(-1) { }
   uint source() const { return _source; }
   uint target() const { return _target; }
   uint offset() const { return _offset; }
  };


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
          edge_type rvalue(gstore_ptr[offset], center, offset);
          ret = rvalue;
        } else if (itype == OUTEDGE) {
          edge_type rvalue(center, gstore_ptr[offset], offset);
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
    uint abs_offset;
    typedef edge_iterator iterator;
    typedef edge_iterator const_iterator;
    typedef edge_type value_type;
    iterator_type itype;

    edge_list(): start_ptr(NULL), end_ptr(NULL), _size(0), source(-1), abs_offset(0), itype(OUTEDGE) { }
    edge_list(uint * _start_ptr, uint * _end_ptr, uint size, uint _source, uint _abs_offset, iterator_type _itype){
      start_ptr = _start_ptr; end_ptr = _end_ptr; abs_offset = _abs_offset;
      _size = size; source = _source; itype = _itype;
    }

    uint size() const { return _size; }

    edge_type operator[](uint i) const{
      ASSERT_LT(i, _size);
      if (itype == OUTEDGE)
        return edge_type(source, *(start_ptr+i), abs_offset+i);
      else return edge_type(*(start_ptr+i), source, abs_offset+i);
    }
    edge_iterator begin() const { return edge_iterator(source, 0, itype, start_ptr ); }
    edge_iterator end() const { return edge_iterator(source, 0, itype, end_ptr ); }
     bool empty() const { return size() == 0; }


  };


  /**
 * CSR/CSC implementation of graph.
 * Assumptions: number of nodes and number of edges are below MAX_UINT32
 */

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
    std::vector<VertexData> *node_vdata_array;
    EdgeData * edge_weights;
    char _color; //not implement yet
    EdgeData _edge;

  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    graph3(){
      num_nodes = _num_edges = 0;
      node_in_edges = node_out_edges = node_in_degrees = node_out_degrees = NULL;
      edge_weights = NULL;
      node_vdata_array = NULL;
      _color = 0; //not implement yet
      undirected = false;
    }

    /**
     * Create a graph with nverts vertices.
     */
    //graph3(size_t nverts) { }

    //graph3(const graph<VertexData, EdgeData>& g) { (*this) = g; }

    // METHODS =================================================================>
    //
    uint * get_node_out_edges(){ return node_out_edges; }
    uint * get_node_in_edges(){ return node_in_edges; }
    /**
     * \brief Resets the graph state.
     */
    void clear() {
       if (node_in_degrees != NULL){
	 delete [] node_in_degrees; node_in_degrees = NULL;
       }
       if (node_out_degrees != NULL){
	 delete [] node_out_degrees; node_out_degrees = NULL;
       }
       if (node_in_edges != NULL){
         delete [] node_in_edges; node_in_edges = NULL; 
       }
       if (node_out_edges != NULL){
         delete [] node_out_edges; node_out_edges = NULL;
       }
       if (edge_weights != NULL){
         delete [] edge_weights; edge_weights = NULL;
       }
}

    void clear_reserve() {
      clear();
    }

    void set_node_vdata_array(const std::vector<VertexData> * _node_vdata_array){
      assert(_node_vdata_array);
      node_vdata_array = (std::vector<VertexData>*)_node_vdata_array;
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

    edge_type find(const vertex_id_type source,
                   const vertex_id_type _target) const {
        /*if (node_out_degrees[source] < node_out_degrees[source+1]){
          std::vector<uint> v(node_out_edges +node_out_degrees[source],node_out_edges+node_out_degrees[source+1]);
          if (binary_search (v.begin(), v.end(), _target))
             return edge_type(source, _target, 0);
        }*/

       for (uint i=node_out_degrees[source]; i< node_out_degrees[source+1]; i++)
          if (node_out_edges[i] == _target)
	     return edge_type(source, _target, i);
	  else if (node_out_edges[i] > _target) //incoming edges asssumed to be sorted
              return edge_type(-1,-1,-1);
       
        return edge_type(-1,-1,-1);
   
    } // end of find

    edge_type reverse_edge(const edge_type& edge) const {
        for (uint i=node_in_degrees[edge.source()]; i< node_in_degrees[edge.source()+1]; i++)
          if (node_in_edges[i] == edge.source())
	     return edge_type(edge.target(), edge.source(), i);

      return edge_type(-1,-1,-1);
    }


    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
       assert(false); //not implemented yet
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
     // ASSERT_GE(num_vertices, num_vertices);
       assert(false); //not implemented yet
      //TODO
    } // End of resize
    
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                          const EdgeData& edata = EdgeData()) {
       assert(false); //not implemented yet
      //return edge_id_type(source, target); //TODO
    } // End of add edge
        
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      return node_vdata_array->at(v);
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      return node_vdata_array->at(v);
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target) {
      edge_type pos = find(source, target);
      if (pos.offset() != (uint)-1) 
      return edge_weights[pos.offset()];
      else return _edge;
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, vertex_id_type target) const {
      edge_type pos = find(source, target);
      if (pos.offset() != -1) 
      return edge_weights[pos.offset()];
      else return _edge;
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_type edge) { 
      ASSERT_NE(edge.offset(), -1);
      return edge_weights[edge.offset()];
    }
    const EdgeData& edge_data(edge_type edge) const {
       ASSERT_NE(edge.offset(), -1);
       return edge_weights[edge.offset()];
    }

    size_t num_in_edges(const vertex_id_type v) const {
      return node_in_degrees[v+1]-node_in_degrees[v];
    }

    size_t num_out_edges(const vertex_id_type v) const {
      return node_out_degrees[v+1]-node_out_degrees[v];
    }

    edge_list_type in_edges(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      return edge_list_type(&node_in_edges[node_in_degrees[v]], &node_in_edges[node_in_degrees[v+1]], num_in_edges(v),v,node_in_degrees[v], INEDGE);  
    }

    edge_list_type out_edges(vertex_id_type v) {
      ASSERT_LT(v, num_nodes);
      return edge_list_type(&node_out_edges[node_out_degrees[v]], &node_out_edges[node_out_degrees[v+1]], num_out_edges(v),v, node_out_degrees[v], OUTEDGE);
    }

    const edge_list_type in_edges(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      return edge_list_type(&node_in_edges[node_in_degrees[v]], &node_in_edges[node_in_degrees[v+1]], num_in_edges(v),v,node_in_degrees[v], INEDGE);  
     }

    const edge_list_type out_edges(vertex_id_type v) const {
      ASSERT_LT(v, num_nodes);
      return edge_list_type(&node_out_edges[node_out_degrees[v]], &node_out_edges[node_out_degrees[v+1]], num_out_edges(v),v,node_out_degrees[v], OUTEDGE);
     }

    const vertex_color_type& color(vertex_id_type vertex) const {
       assert(false); //not implemented yet
      return _color;
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_type vertex) {
      ASSERT_LT(vertex, num_nodes);
       assert(false); //not implemented yet
      return _color;
    }

    vertex_color_type get_color(vertex_id_type vid) const{
       assert(false); //not implemented yet
      return char(0);
    }
    
    void set_color(vertex_id_type vid, vertex_color_type col) {
       assert(false); //not implemented yet
    }
    
    /** \brief This function constructs a heuristic coloring for the 
        graph and returns the number of colors */
    size_t compute_coloring() {
       assert(false); //not implemented yet
      return -1;
    } // end of compute coloring


    /**
     * \brief Check that the colors satisfy a valid coloring of the graph.
     * return true is coloring is valid;
     */
    bool valid_coloring() const {
       assert(false); //not implemented yet
      return true;
    }
    
    /** \brief count the number of times the graph was cleared and rebuilt */
    size_t get_changeid() const {
       assert(false); //not implemented yet
      return 0;
    }

    size_t estimate_sizeof() const {
      /*size_t vlist_size = sizeof(vertices) + sizeof(VertexData) * vertices.capacity();
      size_t vcolor_size = sizeof(vcolors) + sizeof(vertex_color_type) * vcolors.capacity();
      size_t elist_size = edges_tmp.estimate_sizeof(); 
      size_t store_size = gstore.estimate_sizeof();*/

//      printf("graph3: tmplist size: %u, gstoreage size: %u \n", elist_size, store_size);
       return 2*(num_nodes*sizeof(uint)+num_edges*sizeof(int));
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
       assert(false); //not implemented yet
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
       assert(false); //not implemented yet
    } // end of save
   

    /** \brief Load the graph from a file */
    void load(const std::string& filename, bool no_node_data, bool no_edge_data) {
         int rc =array_from_file(filename + ".nodes", node_out_degrees);
	 num_nodes = (rc/4)-1;
         if (!no_node_data){
	    if (node_vdata_array == NULL)
               node_vdata_array = new std::vector<VertexData>();
            node_vdata_array->resize(num_nodes);
         }
 	 logstream(LOG_INFO) << "Read " << num_nodes << " nodes" << std::endl;
         rc = array_from_file(filename + ".edges", node_out_edges);
         _num_edges = (rc/4)-1;
 	 logstream(LOG_INFO) << "Read " << (undirected? _num_edges/2 : _num_edges) << " edges" << std::endl;

         if (!no_edge_data){
           rc = array_from_file(filename + ".weights", edge_weights);
           assert(rc/4 == _num_edges); 
         }
    } // end of loa


    void verify_degrees(const uint * nodes, int len, int n){
       assert(nodes[0] == 0);
       for (int i=0; i< len; i++)
          assert(nodes[i] < n);
    }

    void verify_edges(const uint * edges, int len, int n){
      for (int i=0; i< len; i++)
         assert(edges[i] < (uint)n);
    }

    /** \brief Load the graph from a file */
    void load_directed(const std::string& filename, bool no_node_data, bool no_edge_data) {
      assert(!undirected);
         int rc =array_from_file(filename + ".nodes", node_out_degrees);
	 num_nodes = (rc/4)-1;
         if (!no_node_data){
            if (node_vdata_array == NULL)
              node_vdata_array = new std::vector<VertexData>();
	    node_vdata_array->resize(num_nodes);
         }
         int rc2 =array_from_file(filename + "-r.nodes", node_in_degrees);
         assert(rc == rc2);
         logstream(LOG_INFO) << filename << " Read " << num_nodes << " nodes" << std::endl;
         rc = array_from_file(filename + ".edges", node_out_edges);
         _num_edges = (rc/sizeof(uint));
         rc2 = array_from_file(filename + "-r.edges", node_in_edges);
         assert(rc == rc2);
  	 logstream(LOG_INFO) << filename << " Read " << (undirected? _num_edges/2 : _num_edges) << " edges" << std::endl;
         verify_edges(node_out_edges, _num_edges, num_nodes);
         verify_edges(node_in_edges, _num_edges, num_nodes);
         if (!no_edge_data){
           rc = array_from_file(filename + ".weights", edge_weights);
           assert(rc/sizeof(double) == size_t(_num_edges)); 
           logstream(LOG_INFO) << filename << " Read: " << _num_edges << " weights " << std::endl;
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
       assert(false); //not implemented yet
    }



    /**
     * builds a topological_sort of the graph returning it in topsort. 
     * 
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    bool topological_sort(std::vector<vertex_id_type>& topsort) const {
       assert(false); //not implemented yet
      return true;
    } // end of topological sort

    void set_undirected(){ undirected = true; }

    
  private:    
    /** Internal edge class  */   
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

