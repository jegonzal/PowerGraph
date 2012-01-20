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
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 *
 * Modified by Jay (haijieg@cs.cmu.edu)
 *
 * Describe the interface of a graph. The actual implementation is in graph_storage.hpp.
 *
 * Change interface:
 *  edge_id_type add_edge (vertex_id_type src, vertex_d_type dst, const EdgeData& edata) always return 0 as the temporary edge_id.
 *  May want to change the return type into void.
 *
 * Tested:
 * Test file: graphlab/demoapps/graphTest
 * demo, pagerank
 */


#ifndef GRAPHLAB_GRAPH2_HPP
#define GRAPHLAB_GRAPH2_HPP

#define DEBUG_GRAPH

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

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/util/random.hpp>
#include <graphlab/graph/graph_storage.hpp>
#include <graphlab/graph/graph_storage_noedge.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab { 

  template<typename VertexData, typename EdgeData>
  class graph2 {

    typedef graph_storage<VertexData, EdgeData> gstore_type;
    //typedef graph_storage_noedge<VertexData, EdgeData> gstore_type;

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
 
    /** Represents an edge with source() and target()*/
    typedef typename gstore_type::edge_type edge_type;

    /** The type of the edge list */
    typedef typename gstore_type::edge_list edge_list;
    
    /** Interface for iupdate functor.*/
    typedef typename gstore_type::edge_list edge_list_type;

    /** Edge list type for temporary insertion. */
    typedef typename gstore_type::edge_info edge_info;

  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    graph2() : finalized(false),changeid(0) { }

    /**
     * Create a graph with nverts vertices.
     */
    graph2(size_t nverts) : 
      vertices(nverts),
      vcolors(nverts),
      finalized(false),changeid(0) { }

    void set_is_directed (bool x) {gstore.set_is_directed(x);}
    bool get_is_directed () {return gstore.get_is_directed();}


    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
      finalized = false;
      gstore.clear();
      vertices.clear();
      edges_tmp.clear();
      vcolors.clear();
      ++changeid;
    }

    void clear_reserve() {
      clear();
      edges_tmp.clear();
      std::vector<VertexData>().swap(vertices);
      std::vector<vertex_color_type>().swap(vcolors);
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
      gstore.finalize(vertices.size(), edges_tmp);
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
    edge_type find(const vertex_id_type source,
                   const vertex_id_type target) const {
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
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
      if (finalized)
      {
        logstream(LOG_FATAL)
          << "Attempting add vertex"
          << "to a finalized graph." << std::endl;
        ASSERT_MSG(false, "Add vertex to a finalized graph.");
      }

      vertices.push_back(vdata);
      // Resize edge maps
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
      vcolors.resize(vertices.size());
    } // End of resize
    
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_type add_edge(vertex_id_type source, vertex_id_type target, 
                          const EdgeData& edata = EdgeData()) {
      if (finalized)
      {
        logstream(LOG_FATAL)
          << "Attempting add edge"
          << "to a finalized graph." << std::endl;
        ASSERT_MSG(false, "Add edge to a finalized graph.");
      }


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

      // Add the edge to the set of edge data (this copies the edata)
      edges_tmp.add_edge(source, target, edata);

      // This is not the final edge_id, so we always return 0. 
      return 0;
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
     ASSERT_TRUE(finalized);
     return gstore.edge_data(source, target);
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, vertex_id_type target) const {
     ASSERT_TRUE(finalized);
     return gstore.edge_data(source, target);
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_type edge) { 
      ASSERT_TRUE(finalized);
      return gstore.edge_data(edge);
    }
    const EdgeData& edge_data(edge_type edge) const {
      ASSERT_TRUE(finalized);
      return gstore.edge_data(edge);
    }

    size_t num_in_edges(const vertex_id_type v) const {
      return gstore.num_in_edges(v);
    }

    size_t num_out_edges(const vertex_id_type v) const {
      return gstore.num_out_edges(v);
    }

    edge_list_type in_edges(vertex_id_type v) {
      return gstore.in_edges(v);
    }

    edge_list_type out_edges(vertex_id_type v) {
      return gstore.out_edges(v);
    }

    const edge_list_type in_edges(vertex_id_type v) const {
      return gstore.in_edges(v);
    }

    const edge_list_type out_edges(vertex_id_type v) const {
      return gstore.out_edges(v);
    }


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

    vertex_color_type get_color(vertex_id_type vid) const{
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
        permutation[v] = std::make_pair(-num_in_edges(v), v);
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
        edge_list in_edges = in_edges(vid);
        // Get the neighbor colors
        foreach(edge_type edge, in_edges){
          const vertex_id_type& neighbor_vid = edge.source();
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          if(vertex_color == neighbor_color) return false;
        }
      }
      return true;
    }
    
    /** \brief count the number of times the graph was cleared and rebuilt */
    size_t get_changeid() const {
      return changeid;
    }

    size_t estimate_sizeof() const {
      size_t vlist_size = sizeof(vertices) + sizeof(VertexData) * vertices.capacity();
      size_t vcolor_size = sizeof(vcolors) + sizeof(vertex_color_type) * vcolors.capacity();
      size_t elist_size = edges_tmp.estimate_sizeof(); 
      size_t store_size = gstore.estimate_sizeof();

//      printf("graph2: tmplist size: %u, gstoreage size: %u \n", elist_size, store_size);
      return store_size + vlist_size + vcolor_size + elist_size;
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();    
      // read the vertices and colors
      arc >> vertices
          >> vcolors
          >> gstore
          >> finalized;
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << vertices
          << vcolors
          << gstore
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



    const std::vector<vertex_id_type>& get_out_index_storage() const {
      return gstore.get_csr_src();
    }
    const std::vector<vertex_id_type>& get_in_index_storage() const {
      return gstore.get_csc_dst(); 
    }
    const std::vector<vertex_id_type>& get_out_edge_storage() const {
      return gstore.get_csr_dst();
    }
    const std::vector<vertex_id_type>& get_in_edge_storage() const {
      return gstore.get_csc_src();
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



    /**
     * builds a topological_sort of the graph returning it in topsort. 
     * 
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    bool topological_sort(std::vector<vertex_id_type>& topsort) const {
      topsort.clear();
      topsort.reserve(num_vertices());
    
      std::vector<size_t> indeg;
      indeg.resize(num_vertices());
      std::queue<vertex_id_type> q;
      for (size_t i = 0;i < num_vertices(); ++i) {
        indeg[i] = in_edges(i).size();
        if (indeg[i] == 0) {
          q.push(i);
        }
      }
    
      while (!q.empty()) {
        vertex_id_type v = q.front();
        q.pop();
        topsort.push_back(v);
        foreach(edge_type edge, out_edges(v)) {
          vertex_id_type destv = edge.target();
          --indeg[destv];
          if (indeg[destv] == 0) {
            q.push(destv);
          }
        }
      }
    
      if (q.empty() && topsort.size() != num_vertices()) return false;
      return true;
    } // end of topological sort

    
  private:    
    /** Internal edge class  */   

 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data */
    std::vector<VertexData> vertices;

    /** The vertex colors specified by the user. **/
    std::vector< vertex_color_type > vcolors;  

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
    
    /** increments whenever the graph is cleared. Used to track the
     *  changes to the graph structure  */
    size_t changeid;
  }; // End of class graph2

  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const graph2<VertexData, EdgeData>& graph) {
    // out << "Printing vertex data:\n";
    // for(size_t i = 0; i < graph.num_vertices(); ++i) {
    //   out << "V_" << i << ": D[" << graph.vertex_data(i) << "]\n";      
    // }
    // out << "Printing edge data:\n";
    // for(size_t i = 0; i < graph.num_edges(); ++i) {
    //   out << "(V_" << graph.source(i) << "-> V_" << graph.target(i) << "): "
    //       << "D[" << graph.edge_data(i) << "]\n";      
    // }
    // return out;
    // Print adjacency List
    typedef typename graphlab::graph2<VertexData, EdgeData>::vertex_id_type 
      vertex_id_type;
    typedef typename graphlab::graph2<VertexData, EdgeData>::edge_id_type 
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

