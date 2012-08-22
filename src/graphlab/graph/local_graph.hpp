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

#ifndef GRAPHLAB_LOCAL_GRAPH_HPP
#define GRAPHLAB_LOCAL_GRAPH_HPP


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

#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>
#include <boost/type_traits.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <graphlab/graph/graph_basic_types.hpp>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/graph/graph_storage.hpp>
#include <graphlab/macros_def.hpp>



namespace graphlab { 


  template<typename VertexData, typename EdgeData>
  class json_parser;

  template<typename VertexData, typename EdgeData>
  class local_graph {
    

    /** \internal
     * \brief The type of the graph structure storage of the local_graph. */
    typedef graph_storage<VertexData, EdgeData> gstore_type;
  public:
    
    /** The type of the vertex data stored in the local_graph. */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the local_graph. */
    typedef EdgeData edge_data_type;

    typedef typename gstore_type::edge_info edge_info;

    typedef graphlab::vertex_id_type vertex_id_type;
    typedef graphlab::edge_id_type edge_id_type;

    friend class json_parser<VertexData, EdgeData>;

    
    struct edge_type;
    struct vertex_type;
    struct edge_list_type;

    /** Edge object which provides access to the edge data
     * and information about it.
     */
    struct edge_type {
      local_graph& lgraph_ref;
      typename gstore_type::edge_type e;
      edge_type(local_graph& lgraph_ref, 
                typename gstore_type::edge_type e) : 
        lgraph_ref(lgraph_ref),e(e) { }

      /// \brief Returns a constant reference to the data on the edge.
      const edge_data_type& data() const {
        return lgraph_ref.gstore.edge_data(e);
      }
      /// \brief Returns a reference to the data on the edge.
      edge_data_type& data() {
        return lgraph_ref.gstore.edge_data(e);
      }
      /// \brief Returns the source vertex of the edge.
      vertex_type source() const {
        return vertex_type(lgraph_ref, e.source());
      }
      /// \brief Returns the target vertex of the edge.
      vertex_type target() const {
        return vertex_type(lgraph_ref, e.target());
      }
      /** 
       *  \brief Returns the id of the edge.*/
      edge_id_type id() const {
        return lgraph_ref.gstore.edge_id(e);
      }
    };
    
    /** Vertex object which provides access to the vertex data
     * and information about it.
     */ 
    struct vertex_type {
      local_graph& lgraph_ref;
      lvid_type vid;
      vertex_type(local_graph& lgraph_ref, lvid_type vid):lgraph_ref(lgraph_ref),vid(vid) { }
      
      /// \brief Returns a constant reference to the data on the vertex.
      const vertex_data_type& data() const {
        return lgraph_ref.vertex_data(vid);
      }
      /// \brief Returns a reference to the data on the vertex.
      vertex_data_type& data() {
        return lgraph_ref.vertex_data(vid);
      }
      /// \brief Returns the number of in edges of the vertex.
      size_t num_in_edges() const {
        return lgraph_ref.num_in_edges(vid);
      }
      /// \brief Returns the number of out edges of the vertex.
      size_t num_out_edges() const {
        return lgraph_ref.num_out_edges(vid);
      }
      /// \brief Returns the ID of the vertex. 
      lvid_type id() const {
        return vid;
      }
      /// \brief Returns a list of in edges.
      edge_list_type in_edges() {
        return edge_list_type(lgraph_ref, lgraph_ref.gstore.in_edges(vid));
      }
      /// \brief Returns a list of out edges.
      edge_list_type out_edges() {
        return edge_list_type(lgraph_ref, lgraph_ref.gstore.out_edges(vid));
      }
    };
    
    struct make_edge_type_functor {
      typedef typename gstore_type::edge_type argument_type;
      typedef edge_type result_type;
      local_graph& lgraph_ref;
      make_edge_type_functor(local_graph& lgraph_ref):lgraph_ref(lgraph_ref) { }
      result_type operator() (const argument_type et) const {
        return edge_type(lgraph_ref, et);
      }
    };
    
    /** \brief Represents an iteratable list of edge_types. */
    struct edge_list_type {
      make_edge_type_functor me_functor;
      typename gstore_type::edge_list elist;
      typedef boost::transform_iterator<make_edge_type_functor, typename gstore_type::edge_list::iterator> iterator;
      typedef iterator const_iterator;

      edge_list_type(local_graph& lgraph_ref, typename gstore_type::edge_list elist): me_functor(lgraph_ref), elist(elist) { }
      /// \brief Returns the size of the edge list.
      size_t size() const { return elist.size(); }
      /// \brief Random access to the list elements. 
      edge_type operator[](size_t i) const {return me_functor(elist[i]);}
      /// \brief Returns an iterator to the beginning of the list.
      iterator begin() const { return
          boost::make_transform_iterator(elist.begin(), me_functor); }
      /// \brief Returns an iterator to the end of the list.
      iterator end() const { return
          boost::make_transform_iterator(elist.end(), me_functor); }
      bool empty() const { return elist.empty(); }
    }; // end of class edge_list.



    
  public:

    // CONSTRUCTORS ============================================================>
    
    /** Create an empty local_graph. */
    local_graph() : finalized(false) { }

    /** Create a local_graph with nverts vertices. */
    local_graph(size_t nverts) :
      vertices(nverts),
      finalized(false) { }

    // METHODS =================================================================>

    /**
     * \brief Resets the local_graph state.
     */
    void clear() {
      finalized = false;
      gstore.clear();
      vertices.clear();
      edges_tmp.clear();
    }

    /**
     * \brief Reset the local_graph state and free up the reserved memory.
     */
    void clear_reserve() {
      clear();
      edges_tmp.clear();
      std::vector<VertexData>().swap(vertices);
      gstore.clear_reserve();
    }
    

    /**
     * \brief Finalize the local_graph data structure by
     * sorting edges to maximize the efficiency of graphlab.  
     * This function takes O(|V|log(degree)) time and will 
     * fail if there are any duplicate edges.
     * Detail implementation depends on the type of graph_storage.
     * This is also automatically invoked by the engine at start.
     */
    void finalize() {   
      if(finalized) return;
      graphlab::timer mytimer; mytimer.start();
      gstore.finalize(vertices.size(), edges_tmp);
      logstream(LOG_INFO) << "Graph finalized in " << mytimer.current_time() 
                          << " secs" << std::endl;
      finalized = true;
    } // End of finalize

    /** \brief Get the number of vertices */
    size_t num_vertices() const {
      return vertices.size();
    } // end of num vertices

    /** \brief Get the number of vertices local to this machine */
    size_t local_vertices() const {
      return vertices.size();
    } // end of num vertices

    /** \brief Get the number of edges */
    size_t num_edges() const {
      if (finalized) {
        return gstore.edge_size();
      } else {
        return edges_tmp.size();
      }
    } // end of num edges

    /** \brief Finds an edge. Returns an empty edge if not exists. */
    edge_type find(const lvid_type source,
                   const lvid_type target) const {
      return gstore.find(source, target);
    } // end of find

    /** \brief Finds the reverse of an edge. Returns an empty edge if not exists. */
    edge_type reverse_edge(const edge_type& edge) const {
      return gstore.find(edge.target(), edge.source());
    }


    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    void add_vertex(lvid_type vid, 
                    const VertexData& vdata = VertexData() ) {
        // logstream(LOG_INFO)
        //   << "Attempting add vertex to a finalized local_graph." << std::endl;
        // // ASSERT_MSG(false, "Add vertex to a finalized local_graph.");
      if(vid >= vertices.size()) {
        // Enable capacity doubling if resizing beyond capacity
        if(vid >= vertices.capacity()) {
          const size_t new_size = std::max(2 * vertices.capacity(), 
                                           size_t(vid));
          vertices.reserve(new_size);
        }
        vertices.resize(vid+1);
      }
      vertices[vid] = vdata;    
    } // End of add vertex;

    void reserve(size_t num_vertices) {
      ASSERT_GE(num_vertices, vertices.size());
      vertices.reserve(num_vertices);
    }

    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) {
      ASSERT_GE(num_vertices, vertices.size());
      vertices.resize(num_vertices);
    } // End of resize

    void reserve_edge_space(size_t n) {
      edges_tmp.reserve_edge_space(n);
    }
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared. Should not be called after finalization.
     */
    edge_id_type add_edge(lvid_type source, lvid_type target, 
                          const EdgeData& edata = EdgeData()) {
      if (finalized) {
        logstream(LOG_FATAL)
          << "Attempting add edge to a finalized local_graph." << std::endl;
        ASSERT_MSG(false, "Add edge to a finalized local_graph.");
      }

      if(source == target) {
        logstream(LOG_FATAL) 
          << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
          << "This operation is not permitted in GraphLab!" << std::endl;
        ASSERT_MSG(source != target, "Attempting to add self edge!");
      }
      if(source >= vertices.size() || target >= vertices.size()) 
        add_vertex(std::max(source, target));

      // Add the edge to the set of edge data (this copies the edata)
      edges_tmp.add_edge(source, target, edata);

      // This is not the final edge_id, so we always return 0. 
      return 0;
    } // End of add edge
    
    /**
     * \brief Add edges in block.
     */
    void add_edges(const std::vector<lvid_type>& src_arr, 
                   const std::vector<lvid_type>& dst_arr, 
                   const std::vector<EdgeData>& edata_arr) {
      ASSERT_TRUE((src_arr.size() == dst_arr.size())
                  && (src_arr.size() == edata_arr.size()));
      if (finalized) {
        logstream(LOG_FATAL)
          << "Attempting add edges to a finalized local_graph." << std::endl;
      }

      for (size_t i = 0; i < src_arr.size(); ++i) {
        lvid_type source = src_arr[i];
        lvid_type target = dst_arr[i];
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
            << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
            << "This operation is not permitted in GraphLab!" << std::endl;
          ASSERT_MSG(source != target, "Attempting to add self edge!");
        }
      }
      edges_tmp.add_block_edges(src_arr, dst_arr, edata_arr);
    } // End of add block edges


    /** \brief Returns a vertex of given ID. */
    vertex_type vertex(lvid_type vid) {
      ASSERT_LT(vid, vertices.size());
      return vertex_type(*this, vid);
    }
    
    /** \brief Returns a vertex of given ID. */
    const vertex_type vertex(lvid_type vid) const {
      ASSERT_LT(vid, vertices.size());
      return vertex_type(*this, vid);
    }
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(lvid_type v) {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v. */
    const VertexData& vertex_data(lvid_type v) const {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(lvid_type source, lvid_type target) {
      ASSERT_TRUE(finalized);
      return gstore.edge_data(source, target);
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target. */
    const EdgeData& edge_data(lvid_type source, lvid_type target) const {
      ASSERT_TRUE(finalized);
      return gstore.edge_data(source, target);
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge. */
    EdgeData& edge_data(const edge_type& edge) { 
      ASSERT_TRUE(finalized);
      return gstore.edge_data(edge.e);
    }
    /** \brief Returns a constant reference to the data stored on the edge. */
    const EdgeData& edge_data(const edge_type& edge) const {
      ASSERT_TRUE(finalized);
      return gstore.edge_data(edge.e);
    }



    /** \brief Load the local_graph from an archive */
    void load(iarchive& arc) {
      clear();    
      // read the vertices
      arc >> vertices
          >> gstore
          >> finalized;
    } // end of load

    /** \brief Save the local_graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << vertices
          << gstore
          << finalized;
    } // end of save
    
    /** swap two graphs */
    void swap(local_graph& other) {
      std::swap(vertices, other.vertices);
      std::swap(gstore, other.gstore);
      std::swap(finalized, other.finalized);
    } // end of swap


    /** \brief Load the local_graph from a file */
    void load(const std::string& filename) {
      std::ifstream fin(filename.c_str());
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();
    } // end of load



    /**
     * \brief save the local_graph to the file given by the filename
     */    
    void save(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save



    /**
     * \brief save the adjacency structure to a text file.
     *
     * Save the adjacency structure as a text file in:
     *    src_Id, dest_Id \n
     *    src_Id, dest_Id \n
     * format.
     */
    void save_adjacency(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      ASSERT_TRUE(fout.good());
      for(size_t i = 0; i < num_edges(); ++i) {
        fout << gstore.source(i) << ", " << gstore.target(i) << "\n";
        ASSERT_TRUE(fout.good());
      }          
      fout.close();
    }

 /****************************************************************************
 *                       Internal Functions                                 *
 *                     ----------------------                               *
 * These functions functions and types provide internal access to the       *
 * underlying local_graph representation. They should not be used unless you      *
 * *really* know what you are doing.                                        *
 ****************************************************************************/
    /** \internal
     *  \brief Returns the internal id of the edge*/
    edge_id_type edge_id(const edge_type& edge) const {
      return gstore.edge_id(edge.e);
    }

    /** 
     * \internal
     * \brief Returns the number of in edges of the vertex with the given id. */
    size_t num_in_edges(const lvid_type v) const {
      ASSERT_TRUE(finalized);
      return gstore.num_in_edges(v);
    }

    /** 
     * \internal
     * \brief Returns the number of in edges of the vertex with the given id. */
    size_t num_out_edges(const lvid_type v) const {
      ASSERT_TRUE(finalized);
      return gstore.num_out_edges(v);
    }

    /** 
     * \internal
     * \brief Returns a list of in edges of the vertex with the given id. */
    edge_list_type in_edges(lvid_type v) {
      return edge_list_type(*this, gstore.in_edges(v));
    }

    /** 
     * \internal
     * \brief Returns a list of out edges of the vertex with the given id. */
    edge_list_type out_edges(lvid_type v) {
      return edge_list_type(*this, gstore.out_edges(v));
    }

    /** 
     * \internal
     * \brief Returns the estimated memory footprint of the local_graph. */
    size_t estimate_sizeof() const {
      const size_t vlist_size = sizeof(vertices) + 
        sizeof(VertexData) * vertices.capacity();
      size_t elist_size = edges_tmp.estimate_sizeof();
      size_t store_size = gstore.estimate_sizeof();
      // std::cout << "local_graph: tmplist size: " << (double)elist_size/(1024*1024)
      //           << "  gstoreage size: " << (double)store_size/(1024*1024)
      //           << "  vdata list size: " << (double)vlist_size/(1024*1024)
      //           << std::endl;
      return store_size + vlist_size + elist_size;
    }
    /** \internal
     * \brief Returns the column index of CSR stored in the 
     * internal local_graph storage.
     */
    const std::vector<lvid_type>& get_out_index_storage() const {
      return gstore.get_csr_src();
    }
    /** \internal
     * \brief Returns the row index of CSC stored in the 
     * internal local_graph storage.
     */
    const std::vector<lvid_type>& get_in_index_storage() const {
      return gstore.get_csc_dst(); 
    }
    /** \internal
     * \brief Returns the row pointer of CSR stored in the 
     * internal local_graph storage.
     */
    const std::vector<lvid_type>& get_out_edge_storage() const {
      return gstore.get_csr_dst();
    }

    /** \internal
     * \brief Returns the column pointer of CSC stored in the 
     * internal local_graph storage.
     */
    const std::vector<lvid_type>& get_in_edge_storage() const {
      return gstore.get_csc_src();
    }
    /** \internal
     * \brief Returns the reference of edge data list stored in the
     * internal local_graph storage.
     */
    const std::vector<EdgeData> & get_edge_data_storage() const {
      return gstore.get_edge_data();
    }

    /** \internal
     * \brief For debug purpose, returns the largest vertex id in the edges_tmp
     */ 
    const lvid_type maxlvid() const {
      if (edges_tmp.size()) {
        lvid_type max(0);
        foreach(lvid_type i, edges_tmp.source_arr)
         max = std::max(max, i); 
        foreach(lvid_type i, edges_tmp.target_arr)
         max = std::max(max, i); 
        return max;
      } else {
        return lvid_type(-1);
      }
    }
   
  private:    
    /** Internal edge class  */   

 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data */
    std::vector<VertexData> vertices;

    /** Stores the edge data and edge relationships. */
    gstore_type gstore;

    /** The edge data is a vector of edges where each edge stores its
        source, destination, and data. Used for temporary storage. The
        data is transferred into CSR+CSC representation in
        Finalize. This will be cleared after finalized.*/
    edge_info edges_tmp;
   
    /** Mark whether the local_graph is finalized.  Graph finalization is a
        costly procedure but it can also dramatically improve
        performance. */
    bool finalized;

  }; // End of class local_graph


  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const local_graph<VertexData, EdgeData>& local_graph) {
    for(lvid_type vid = 0; vid < local_graph.num_vertices(); ++vid) {
      foreach(edge_id_type eid, local_graph.out_edge_ids(vid))
        out << vid << ", " << local_graph.target(eid) << '\n';
    }
    return out;
  }
} // end of namespace graphlab


namespace std {
  /**
   * Swap two graphs
   */
  template<typename VertexData, typename EdgeData>
  inline void swap(graphlab::local_graph<VertexData,EdgeData>& a,
                   graphlab::local_graph<VertexData,EdgeData>& b) {
    a.swap(b);
  } // end of swap

}; // end of namespace std

#include <graphlab/macros_undef.hpp>
#endif

