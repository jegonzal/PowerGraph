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
 * Modified by Jay (haijieg@cs.cmu.edu)
 */


#ifndef GRAPHLAB_GRAPH_HPP
#define GRAPHLAB_GRAPH_HPP

// #define DEBUG_GRAPH

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

#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>
#include <boost/type_traits.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/lambda/construct.hpp>
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
  class graph {
    
    typedef graph_storage<VertexData, EdgeData> gstore_type;
  public:
    
    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData edge_data_type;

    typedef typename gstore_type::edge_info edge_info;

    typedef graphlab::vertex_id_type vertex_id_type;
    typedef graphlab::edge_id_type edge_id_type;
    
    struct edge_type;
    struct vertex_type;
    struct edge_list_type;

    struct edge_type {
      graph& g;
      typename gstore_type::edge_type e;
      edge_type(graph& g, typename gstore_type::edge_type e):g(g),e(e) { }

      const edge_data_type& data() const {
        return g.gstore.edge_data(e);
      }
      edge_data_type& data() {
        return g.gstore.edge_data(e);
      }
      vertex_type source() {
        return vertex_type(g, e.source());
      }
      vertex_type target() {
        return vertex_type(g, e.target());
      }

      edge_id_type id() {
        return g.gstore.edge_id(e);
      }
    };
    
    struct vertex_type {
      graph& g;
      lvid_type vid;
      vertex_type(graph& g, lvid_type vid):g(g),vid(vid) { }
      
      const vertex_data_type& data() const {
        return g.vertex_data(vid);
      }

      vertex_data_type& data() {
        return g.vertex_data(vid);
      }

      size_t num_in_edges() const {
        return g.num_in_edges(vid);
      }

      size_t num_out_edges() const {
        return g.num_out_edges(vid);
      }

      lvid_type id() const {
        return vid;
      }

      edge_list_type in_edges() {
        return edge_list_type(g, g.gstore.in_edges(vid));
      }

      edge_list_type out_edges() {
        return edge_list_type(g, g.gstore.out_edges(vid));
      }
    };
    
    struct make_edge_type_functor {
      typedef typename gstore_type::edge_type argument_type;
      typedef edge_type result_type;
      graph& g;
      make_edge_type_functor(graph& g):g(g) { }
      result_type operator() (const argument_type et) const {
        return edge_type(g, et);
      }
    };
    
    // Represents an iteratable list of edge_types.
    struct edge_list_type {
      make_edge_type_functor me_functor;
      typename gstore_type::edge_list elist;
      typedef boost::transform_iterator<make_edge_type_functor, typename gstore_type::edge_list::iterator> iterator;
      typedef iterator const_iterator;
      // Cosntruct an edge_list with begin and end.
      edge_list_type(graph& g, typename gstore_type::edge_list elist): me_functor(g), elist(elist) { }
      size_t size() const { return elist.size(); }
      edge_type operator[](size_t i) const {return me_functor(elist[i]);}
      iterator begin() const { return
          boost::make_transform_iterator(elist.begin(), me_functor); }
      iterator end() const { return
          boost::make_transform_iterator(elist.end(), me_functor); }
      bool empty() const { return elist.empty(); }
    }; // end of class edge_list.



    
  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    graph() : finalized(false) { }

    /**
     * Create a graph with nverts vertices.
     */
    graph(size_t nverts) :
      vertices(nverts),
      finalized(false) { }

    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
      finalized = false;
      gstore.clear();
      vertices.clear();
      edges_tmp.clear();
    }

    void clear_reserve() {
      clear();
      edges_tmp.clear();
      std::vector<VertexData>().swap(vertices);
      gstore.clear_reserve();
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
      // std::cout << "considering finalize" << std::endl;
      // check to see if the graph is already finalized
      if(finalized) return;
      //      std::cout << "Finalizing" << std::endl;
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

    /** \brief Finds an edge.
        The value of the first element of the pair will be true if an 
        edge from src to target is found and false otherwise. If the 
        edge is found, the edge ID is returned in the second element of the pair. */
    edge_type find(const lvid_type source,
                   const lvid_type target) const {
      return gstore.find(source, target);
    } // end of find

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
        //   << "Attempting add vertex to a finalized graph." << std::endl;
        // // ASSERT_MSG(false, "Add vertex to a finalized graph.");
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
     * existing data will be cleared.
     */
    edge_id_type add_edge(lvid_type source, lvid_type target, 
                          const EdgeData& edata = EdgeData()) {
      if (finalized) {
        logstream(LOG_FATAL)
          << "Attempting add edge to a finalized graph." << std::endl;
        ASSERT_MSG(false, "Add edge to a finalized graph.");
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
     * \brief Add edge in block
     */
    void add_edges(const std::vector<lvid_type>& src_arr, 
                   const std::vector<lvid_type>& dst_arr, 
                   const std::vector<EdgeData>& edata_arr) {
      ASSERT_TRUE((src_arr.size() == dst_arr.size())
                  && (src_arr.size() == edata_arr.size()));
      if (finalized) {
        logstream(LOG_FATAL)
          << "Attempting add edges to a finalized graph." << std::endl;
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


    vertex_type vertex(lvid_type vid) {
      return vertex_type(*this, vid);
    }
    
    const vertex_type vertex(lvid_type vid) const {
      return vertex_type(*this, vid);
    }
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(lvid_type v) {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
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
        edge source->target */
    const EdgeData& edge_data(lvid_type source, lvid_type target) const {
      ASSERT_TRUE(finalized);
      return gstore.edge_data(source, target);
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(const edge_type& edge) { 
      ASSERT_TRUE(finalized);
      return gstore.edge_data(edge.e);
    }
    const EdgeData& edge_data(const edge_type& edge) const {
      ASSERT_TRUE(finalized);
      return gstore.edge_data(edge.e);
    }

    edge_id_type edge_id(const edge_type& edge) const {
      return gstore.edge_id(edge.e);
    }

    size_t num_in_edges(const lvid_type v) const {
      ASSERT_TRUE(finalized);
      return gstore.num_in_edges(v);
    }

    size_t num_out_edges(const lvid_type v) const {
      ASSERT_TRUE(finalized);
      return gstore.num_out_edges(v);
    }

    edge_list_type in_edges(lvid_type v) {
      return edge_list_type(*this, gstore.in_edges(v));
    }

    edge_list_type out_edges(lvid_type v) {
      return edge_list_type(*this, gstore.out_edges(v));
    }



    size_t estimate_sizeof() const {
      const size_t vlist_size = sizeof(vertices) + 
        sizeof(VertexData) * vertices.capacity();
      size_t elist_size = edges_tmp.estimate_sizeof();
      size_t store_size = gstore.estimate_sizeof();
      std::cout << "graph: tmplist size: " << (double)elist_size/(1024*1024)
                << "  gstoreage size: " << (double)store_size/(1024*1024)
                << "  vdata list size: " << (double)vlist_size/(1024*1024)
                << std::endl;
      return store_size + vlist_size + elist_size;
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();    
      // read the vertices
      arc >> vertices
          >> gstore
          >> finalized;
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << vertices
          << gstore
          << finalized;
    } // end of save
    
    /** swap two graphs */
    void swap(graph& other) {
      std::swap(vertices, other.vertices);
      std::swap(gstore, other.gstore);
      std::swap(finalized, other.finalized);
    } // end of swap


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



    const std::vector<lvid_type>& get_out_index_storage() const {
      return gstore.get_csr_src();
    }
    const std::vector<lvid_type>& get_in_index_storage() const {
      return gstore.get_csc_dst(); 
    }
    const std::vector<lvid_type>& get_out_edge_storage() const {
      return gstore.get_csr_dst();
    }
    const std::vector<lvid_type>& get_in_edge_storage() const {
      return gstore.get_csc_src();
    }

    const std::vector<EdgeData> & get_edge_data_storage() const {
      return gstore.get_edge_data();
    }

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
      for(size_t i = 0; i < num_edges; ++i) {
        fout << gstore.source(i) << ", " << gstore.target(i) << "\n";
        ASSERT_TRUE(fout.good());
      }          
      fout.close();
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
   
    /** Mark whether the graph is finalized.  Graph finalization is a
        costly procedure but it can also dramatically improve
        performance. */
    bool finalized;
  }; // End of class graph


  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const graph<VertexData, EdgeData>& graph) {
    for(lvid_type vid = 0; vid < graph.num_vertices(); ++vid) {
      foreach(edge_id_type eid, graph.out_edge_ids(vid))
        out << vid << ", " << graph.target(eid) << '\n';      
    }
    return out;
  }
} // end of namespace graphlab


namespace std {
  /**
   * Swap two graphs
   */
  template<typename VertexData, typename EdgeData>
  inline void swap(graphlab::graph<VertexData,EdgeData>& a,
                   graphlab::graph<VertexData,EdgeData>& b) {
    a.swap(b);
  } // end of swap

}; // end of namespace std

#include <graphlab/macros_undef.hpp>
#endif

