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


#include <boost/bind.hpp>
#include <boost/unordered_set.hpp>



#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>

#include <graphlab/extern/metis/metis.hpp>

#include <graphlab/util/random.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab { 


  /**
   * \brief the partitioning methods.
   * 
   * The partition method struct contains the partition_method_enum
   * which desribes the partitioning method supported by the graphlab
   * graph partitioning framework.  Graphlab links to the metis graph
   * partitioning library and also provides its own partitioning
   * algorithms.  
   */
  struct partition_method {

    /**
      \brief the partition methods
     */
    enum partition_method_enum {
      PARTITION_RANDOM, /**< Vertices are randomly assigned to partitions. Every partition
                      will have roughly the same number of vertices. (between
                      #vertices / nparts and (#vertices / nparts + 1)) */
      PARTITION_METIS,  /**< Partitions the graph using the
                    <a href="http://glaros.dtc.umn.edu/gkhome/views/metis"> METIS</a>
                    graph partitioning package. <b>A modified version of METIS is
                    included with GraphLab, so contact us if you encounter any
                    issues.</b>*/
      PARTITION_BFS, /**< Partitions the graph using Breadth First Searches. A random
                    root vertex is selected and inserted into the first
                    partition. Additional vertices are inserted into the partition
                    by starting a BFS from the root vertex.  When the partition
                    reaches capacity: (#vertices / nparts), the procedure restarts
                    with a next partition with a new random root vertex. */
      PARTITION_EDGE_NUM, /**< Partitions the vertices such that every partition 
                        has roughly the same number of edges. */
    };

    /// Converts a partition_method_enum to a string
    static std::string enum_to_string(partition_method_enum val) {
      switch(val) {
        case PARTITION_RANDOM:
          return "random";
        case PARTITION_METIS:
          return "metis";
        case PARTITION_BFS:
          return "bfs";
        case PARTITION_EDGE_NUM:
          return "edge_num";
        default:
          return "";
      }
    }
    
    /// Converts a string to a partition_method_enum. Returns true on success
    static bool string_to_enum(std::string s, partition_method_enum &val) {
      if (s == "random") {
        val = PARTITION_RANDOM;
        return true;
      }
      else if (s == "metis") {
        val = PARTITION_METIS;
        return true;
      }
      else if (s == "bfs") {
        val = PARTITION_BFS;
        return true;
      }
      else if (s == "edge_num") {
        val = PARTITION_EDGE_NUM;
        return true;
      }
      return false;
    }
  };

  /// The type of a vertex is a simple size_t
  typedef uint32_t vertex_id_t;

  /// The type of an edge id 
  typedef uint32_t edge_id_t;

  /// Type for vertex colors 
  typedef uint8_t vertex_color_type;


  /** This class defines a set of edges */
  class edge_list {
  public:
    typedef const edge_id_t* iterator; // Should not be used
    typedef const edge_id_t* const_iterator;
    typedef edge_id_t value_type;
  private:
    const edge_id_t* begin_ptr; // Points to first element
    const edge_id_t* end_ptr; // One past end   
  public:
    /** \brief Construct an empty edge list */
    edge_list() : begin_ptr(NULL), end_ptr(NULL) { }
    /** \brief Construct an edge list from an std vector */
    edge_list(const std::vector<edge_id_t>& edges) :
      begin_ptr(&(*edges.begin())), 
      end_ptr(begin_ptr + edges.size()) { }
    /** \brief Construct an edge list from an in memory array */
    edge_list(const edge_id_t* begin_ptr, size_t len) :
      begin_ptr(begin_ptr),  end_ptr(begin_ptr + len) { }

    /** \brief Get the size of the edge list */
    size_t size() const { return (size_t)(end_ptr - begin_ptr); }

    /** \brief Get the ith edge in the edge list */
    edge_id_t operator[](size_t i) const {
      ASSERT_LT(i,  size());
      return *(begin_ptr + i);
    }

    /** \brief Returns a pointer to the start of the edge list */
    const edge_id_t* begin() const {
      return begin_ptr;
    }

    /** \brief Returns a pointer to the end of the edge list */
    const edge_id_t* end() const {
      return end_ptr;
    } 

    /** \brief Fill a vector with edge id list */
    void fill_vector(std::vector<edge_id_t>& lvalue) const {
      lvalue.clear();
      foreach(edge_id_t eid, *this) lvalue.push_back(eid);    
    }

    /** \brief test if the edge list is empty */
    bool empty() const { return size() == 0; }

  }; // End of edge list




  
  template<typename VertexData, typename EdgeData> class graph;


  // CLASS GRAPH ==============================================================>
  /**
    \brief The GraphLab primary Graph container templatized over the vertex and edge types.
    
    Every vertex and edge in the graph is assigned a unique integer
    ID.  The type of the vertex id
    is <code>graphlab::vertex_id_t</code> and the type of the edge id
    is <code>graphlab::edge_id_t</code>. Both <code>vertex_id_t</code>
    and <code>edge_id_t</code> are currently defined
    as <code>uint32_t</code>.  While this limits the graphs to 4
    billion vertices it also helps reduce the storage overhead. We
    encourage users to use the <code>vertex_id_t</code>
    and <code>edge_id_t</code> types as they may change in larger
    distributed systems.

    <h2> Graph Creation </h2>
  
    Vertices and edges are added using the graph::add_vertex()
    and graph::add_edge() member functions:
  
  
  \code
  vertex_id_t graph::add_vertex(const VertexData& vdata = VertexData()) 
  edge_id_t graph::add_edge(vertex_id_t source, vertex_id_t target, const EdgeData& edata = EdgeData()) 
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

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
      
    typedef edge_list edge_list_type;
  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    graph() : finalized(true),changeid(0) {  }

    /**
     * Create a graph with nverts vertices.
     */
    graph(size_t nverts) : 
      vertices(nverts),
      in_edges(nverts), out_edges(nverts), vcolors(nverts),
      finalized(true),changeid(0) { }

    graph(const graph<VertexData, EdgeData>& g) { (*this) = g; }

    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
      vertices.clear();
      edges.clear();
      in_edges.clear();
      out_edges.clear();
      vcolors.clear();
      finalized = true;
      ++changeid;
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
      //      std::cout << "considering finalize" << std::endl;
      // check to see if the graph is already finalized
      if(finalized) return;
      //      std::cout << "Finalizing" << std::endl;
      typedef std::vector< edge_id_t > edge_set;
      edge_id_less_functor less_functor(this);      
      // Sort all in edges sets
      //      foreach(edge_set& eset, in_edges) {
#pragma omp parallel for
      for(ssize_t i = 0; i < ssize_t(in_edges.size()); ++i) {
        edge_set& eset(in_edges[i]);
        // Sort the edge vector
        std::sort(eset.begin(),
                  eset.end(),
                  less_functor);
        // Test for duplicate edges
        if (eset.size() > 1) {
          for(size_t j = 0; j < eset.size()-1; ++j) {
            // Duplicate edge test
            if(!edge_id_less(eset[j], eset[j+1])) {
              logstream(LOG_FATAL)
                << "Duplicate edge "
                << "(" << source(eset[j]) << ", " << target(eset[j]) << ") "
                << "found!  GraphLab does not support graphs "
                << "with duplicate edges." << std::endl;
            }
          }
        }  
      } // end of for loop


      // Sort all out edges sets
      //      foreach(edge_set& eset, out_edges) {
#pragma omp parallel for
      for(ssize_t i = 0; i < ssize_t(out_edges.size()); ++i) {
        edge_set& eset(out_edges[i]);
        std::sort(eset.begin(),
                  eset.end(),
                  less_functor);
        // Test for dupliate edges
        if (eset.size() > 1) {
          for(size_t j = 0; j < eset.size()-1; ++j) {
            // Duplicate edge test
            if(!edge_id_less(eset[j], eset[j+1])) {
              logstream(LOG_FATAL)
                << "Duplicate edge "
                << "(" << source(eset[j]) << ", " << target(eset[j]) << ") "
                << "found!  GraphLab does not support graphs "
                << "with duplicate edges." << std::endl;
            }
          }
        }
      }
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
      return edges.size();
    } // end of num edges


    /** \brief Get the number of in edges of a particular vertex */
    size_t num_in_neighbors(vertex_id_t v) const {
      ASSERT_LT(v, vertices.size());
      return in_edges[v].size();
    } // end of num vertices
    
    /** \brief Get the number of out edges of a particular vertex */
    size_t num_out_neighbors(vertex_id_t v) const  {
      ASSERT_LT(v, vertices.size());
      return out_edges[v].size();
    } // end of num vertices

    /** \brief Finds an edge.
    The value of the first element of the pair will be true if an 
    edge from src to target is found and false otherwise. If the 
    edge is found, the edge ID is returned in the second element of the pair. */
    std::pair<bool, edge_id_t>
    find(vertex_id_t source, vertex_id_t target) const {
      ASSERT_LT(source, in_edges.size());
      ASSERT_LT(target, out_edges.size());
      // Check the base case that the souce or target have no edges
      if (in_edges[target].size() == 0 ||
          out_edges[source].size() == 0) {
        return std::make_pair(false,-1);
      } else if(finalized) { // O( log degree ) search ========================>
        // if it is finalized then do the search using a binary search
        // If their are fewer in edges into the target search the in
        // edges
        if(in_edges[target].size() < out_edges[source].size()) {
          // search the source vertices for the edge
          size_t index = binary_search(in_edges[target], source, target);
          if(index < in_edges[target].size())
            return std::make_pair(true, in_edges[target][index]);
          else
            return std::make_pair(false,-1);
        } else { // If their are fewer edges out of the source binary
                 // search there
          // search the source vertices for the edge
          size_t index = binary_search(out_edges[source], source, target);
          if(index < out_edges[source].size())
            return std::make_pair(true, out_edges[source][index]);
          else
            return std::make_pair(false,-1);
        }
      } else { // O( degree ) search ==========================================>
        // if there are few in edges at the target search there
        if(in_edges[target].size() < out_edges[source].size()) {
          // linear search the in_edges at the target 
          foreach(edge_id_t eid, in_edges[target]) {
            ASSERT_LT(eid, edges.size());
            if(edges[eid].source() == source 
               && edges[eid].target() == target) {
              return std::make_pair(true, eid);
            }
          }
          return std::make_pair(false, -1);
        } else { // fewer out edges at the source
          // linear search the out_edges at the source
          foreach(edge_id_t eid, out_edges[source]) {
            ASSERT_LT(eid, edges.size());
            if(edges[eid].source() == source 
               && edges[eid].target() == target) {
              return std::make_pair(true, eid);
            }
          }
          return std::make_pair(false, -1);
        }
      } // End of else 
    } // end of find

    
    /** \brief A less safe version of find. 
        Returns the edge_id of an edge from src to target exists. 
        Assertion failure otherwise. */
    edge_id_t edge_id(vertex_id_t source, vertex_id_t target) const {
      std::pair<bool, edge_id_t> res = find(source, target);
      // The edge must exist
      ASSERT_TRUE(res.first);
      ASSERT_LT(res.second, edges.size());
      return res.second;
    } // end of edge_id

    
    /** \brief Returns the edge ID of the edge going in the opposite direction. 
        Assertion failure if such an edge is not found.  */
    edge_id_t rev_edge_id(edge_id_t eid) const {
      ASSERT_LT(eid, edges.size());
      vertex_id_t source = edges[eid].source();
      vertex_id_t target = edges[eid].target();    
      return edge_id(target, source);
    } // end of rev_edge_id

    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_t add_vertex(const VertexData& vdata = VertexData() ) {
      vertices.push_back(vdata);
      // Resize edge maps
      out_edges.resize(vertices.size());
      in_edges.resize(vertices.size());
      vcolors.resize(vertices.size());
      return (vertex_id_t)vertices.size() - 1;
    } // End of add vertex;


    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) {
      ASSERT_GE(num_vertices, vertices.size());
      vertices.resize(num_vertices);
      // Resize edge maps
      out_edges.resize(vertices.size());
      in_edges.resize(vertices.size());
      vcolors.resize(vertices.size());
    } // End of resize
    
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_t add_edge(vertex_id_t source, vertex_id_t target, 
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
          << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
          << "This operation is not permitted in GraphLab!" << std::endl;
        ASSERT_MSG(source != target, "Attempting to add self edge!");
      }

      // Add the edge to the set of edge data (this copies the edata)
      edges.push_back( edge( source, target, edata ) );

      // Add the edge id to in and out edge maps
      edge_id_t edge_id = (edge_id_t)edges.size() - 1;
      in_edges[target].push_back(edge_id);
      out_edges[source].push_back(edge_id);

      // Determine if the graph is still finalized A graph is
      // finalized if it was finalized and the newly added edge_id is
      // in the correct location in the in and out edge lists (which
      // is true if either the lists only contain a single element or
      // the last two elements are in the correct order).
      finalized = finalized &&
        ((in_edges[target].size() < 2) ||
         edge_id_less(*(in_edges[target].end()-2),
                      *(in_edges[target].end()-1))) &&
        ((out_edges[source].size() < 2) ||
         edge_id_less(*(out_edges[source].end()-2),
                      *(out_edges[source].end()-1)));
      return edge_id;
    } // End of add edge
        
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(vertex_id_t v) {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_t v) const {
      ASSERT_LT(v, vertices.size());
      return vertices[v];
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_t source, vertex_id_t target) {
      ASSERT_LT(source, vertices.size());
      ASSERT_LT(target, vertices.size());
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      ASSERT_TRUE(ans.first);
      // the edge id should be valid!
      ASSERT_LT(ans.second, edges.size());
      return edges[ans.second].data();
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_t source, vertex_id_t target) const {
      ASSERT_LT(source, vertices.size());
      ASSERT_LT(target, vertices.size());
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      ASSERT_TRUE(ans.first);
      // the edge id should be valid!
      ASSERT_LT(ans.second, edges.size());
      return edges[ans.second].data();
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_id_t edge_id) { 
      ASSERT_LT(edge_id, edges.size());
      return edges[edge_id].data();
    }
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(edge_id_t edge_id) const {
      ASSERT_LT(edge_id, edges.size());
      return edges[edge_id].data();
    }

    /** \brief Returns the source vertex of an edge. */
    vertex_id_t source(edge_id_t edge_id) const {
      //      ASSERT_LT(edge_id, edges.size());
      return edges[edge_id].source();
    }

    /** \brief Returns the destination vertex of an edge. */
    vertex_id_t target(edge_id_t edge_id) const {
      //      ASSERT_LT(edge_id, edges.size());
      return edges[edge_id].target();    
    }
    
    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_t vertex) const {
      ASSERT_LT(vertex, vertices.size());
      return vcolors[vertex];
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_t vertex) {
      ASSERT_LT(vertex, vertices.size());
      return vcolors[vertex];
    }


    /** \brief This function constructs a heuristic coloring for the 
    graph and returns the number of colors */
    size_t compute_coloring() {
      // Reset the colors
      for(vertex_id_t v = 0; v < num_vertices(); ++v) color(v) = 0;
      // construct a permuation of the vertices to use in the greedy
      // coloring. \todo Should probably sort by degree instead when
      // constructing greedy coloring.
      std::vector<std::pair<vertex_id_t, vertex_id_t> > 
	permutation(num_vertices());

      for(vertex_id_t v = 0; v < num_vertices(); ++v) 
        permutation[v] = std::make_pair(-num_in_neighbors(v), v);
      //      std::random_shuffle(permutation.begin(), permutation.end());
      std::sort(permutation.begin(), permutation.end());
      // Recolor
      size_t max_color = 0;
      std::set<vertex_color_type> neighbor_colors;
      for(size_t i = 0; i < permutation.size(); ++i) {
        neighbor_colors.clear();
        const vertex_id_t& vid = permutation[i].second;
        // Get the neighbor colors
        foreach(edge_id_t eid, in_edge_ids(vid)){
          const vertex_id_t& neighbor_vid = source(eid);
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          neighbor_colors.insert(neighbor_color);
        }
        foreach(edge_id_t eid, out_edge_ids(vid)){
          const vertex_id_t& neighbor_vid = target(eid);
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
      for(vertex_id_t vid = 0; vid < num_vertices(); ++vid) {
        const vertex_color_type& vertex_color = color(vid);
        edge_list in_edges = in_edge_ids(vid);
        // Get the neighbor colors
        foreach(edge_id_t eid, in_edges){
          const vertex_id_t& neighbor_vid = source(eid);
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          if(vertex_color == neighbor_color) return false;
        }
      }
      return true;
    }
    
    
    /** \brief Return the edge ids of the edges arriving at v */
    edge_list in_edge_ids(vertex_id_t v) const {
      ASSERT_LT(v, in_edges.size());
      return edge_list(in_edges[v]);
    } // end of in edges    

    /** \brief Return the edge ids of the edges leaving at v */
    edge_list out_edge_ids(vertex_id_t v) const {
      ASSERT_LT(v, out_edges.size());
      return edge_list(out_edges[v]);
    } // end of out edges
    
    /** \brief Get the set of in vertices of vertex v */
    std::vector<vertex_id_t> in_vertices(vertex_id_t v) const {
      std::vector<vertex_id_t> ret;
      foreach(edge_id_t eid, in_edges[v]) {
        ret.push_back(edges[eid].source());
      }
      return ret;
    }
    
    /** \brief Get the set of out vertices of vertex v */
    std::vector<vertex_id_t> out_vertices(vertex_id_t v) const {
      std::vector<vertex_id_t> ret;
      foreach(edge_id_t eid, out_edges[v]) {
        ret.push_back(edges[eid].target());
      }
      return ret;
    }
    
    
    /** \brief count the number of times the graph was cleared and rebuilt */
    size_t get_changeid() const {
      return changeid;
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();    
      // read the vertices and colors
      arc >> vertices
          >> edges
          >> in_edges
          >> out_edges
          >> vcolors
          >> finalized;
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << vertices
          << edges
          << in_edges
          << out_edges
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
      for(size_t i = 0; i < edges.size(); ++i) {
        fout << edges[i].source() << ", " << edges[i].target() << "\n";
        ASSERT_TRUE(fout.good());
      }          
      fout.close();
    }

    /**
     * \brief partition the graph with roughly the same number of edges for
     * each part. Equivalent to calling partition() with the 
     * partition_method::PARTITION_EDGE_NUM parameter.
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     */
    void edge_num_partition(size_t nparts, 
			    std::vector<vertex_id_t>& vertex2part){
      vertex2part.resize(num_vertices());
      size_t e = 2 * num_edges();
      size_t edge_per_part = e / nparts;

      size_t curpart = 0;
      std::vector<size_t> parts;
      parts.resize(nparts, 0);
      
      for (vertex_id_t i = 0; i< num_vertices(); i++){
        uint32_t ne = (uint32_t)(out_edge_ids(i).size() + in_edge_ids(i).size());
        vertex2part[i] = (uint32_t)curpart;
        parts[curpart]+= ne;
           
        if (parts[curpart] >= edge_per_part  && curpart < nparts-1){
          curpart++;
        }
      }
    }


    /**
     * \brief Randomly assign vertices to partitions.  This will assign
     * vertices evenly to each partition.
     * Equivalent to calling partition() with the 
     * partition_method::PARTITION_RANDOM parameter
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     */
    void random_partition(size_t nparts, std::vector<vertex_id_t>& vertex2part) {
      vertex2part.resize(num_vertices());
      for (vertex_id_t i = 0;i < num_vertices(); ++i) {
        vertex2part[i] = (vertex_id_t)(i % nparts);
      }
      random::shuffle(vertex2part.begin(), vertex2part.end());
    }


    /**
     * \brief Use a modified version of the METIS library to partition the
     * graph. Equivalent to calling partition() with the 
     * partition_method::PARTITION_METIS paramter
     *
     * \param numparts The number of parts to partition into
     * \param[out] ret_part A vector providing a vertex_id -> partition_id mapping
     *
     * The metis library is described in:
     *
     *   "A Fast and Highly Quality Multilevel Scheme for Partitioning
     *   Irregular Graphs”. George Karypis and Vipin Kumar. SIAM
     *   Journal on Scientific Computing, Vol. 20, No. 1, pp. 359—392,
     *   1999.
     *
     * We have modified an alpha version (5.0) to work with the
     * GraphLab framework.  Therefore users having trouble with this
     * function or the included Metis source should direct concerns to
     * the contact information provided at:
     *
     *   http://www.select.cs.cmu.edu/code
     */
    void metis_partition(size_t numparts , std::vector<vertex_id_t>& ret_part) {
      if (numparts == 1) {
        ret_part.assign(num_vertices(), 0);
        return;
      }
      // Determine parameters needed to construct the partitioning
      metis::idxtype numverts(num_vertices());
      ASSERT_GT(numverts, 0);
      // Compute the number of edges 
      metis::idxtype numedges(num_edges());

      // allocate metis data structures
      metis::idxtype* vweight = new metis::idxtype[numverts];
      ASSERT_NE(vweight, NULL);    
      metis::idxtype* xadj = new metis::idxtype[numverts + 1];
      ASSERT_NE(xadj, NULL);
      metis::idxtype* adjacency = new metis::idxtype[2 * numedges];
      ASSERT_NE(adjacency, NULL);
      metis::idxtype* eweight = NULL;
      //       if(weighted) {
      //         eweight = new idxtype[numedges];
      //         assert(eweigth != NULL);
      //       }
      metis::idxtype* res = new metis::idxtype[numverts];   
      ASSERT_NE(res, NULL);

      // Pass through vertices filling in the metis data structures
      size_t offset = 0;
      for(vertex_id_t u = 0; u < num_vertices(); ++u) {
        // Update vertex weight
        // Set weight
        vweight[u] = 1;
        // Update the offset
        xadj[u] = offset;
        // Fill the the adjacency data
      
        std::set<vertex_id_t> neighbors;
        foreach(edge_id_t eid, out_edge_ids(u)) {
          neighbors.insert(target(eid));
        }
        foreach(edge_id_t eid, in_edge_ids(u)) {
          neighbors.insert(source(eid));
        }
        foreach(vertex_id_t vid, neighbors) {
          if (vid == u) continue;
          adjacency[offset] = vid;
          ASSERT_GE(adjacency[offset], 0);
          offset++;
          ASSERT_GE(offset, 0);
        }
    
      } // end of data structure creation
      
      // Set the last entry in xadj to the end of the adjacency array
      xadj[numverts] = offset;
    
      // Set additional metis flags
      /**
       * 0 No weights (vwgts and adjwgt are NULL) 
       * 1 Weights on the edges only (vwgts = NULL) 
       * 2 Weights on the vertices only (adjwgt = NULL) 
       * 3 Weights both on vertices and edges. 
       */
      metis::idxtype weightflag = 2;
      // 0 for C-style numbering starting at 0 (1 for fortran style)
      metis::idxtype numflag = 0;
      // the number of parts to cut into
      metis::idxtype nparts = numparts;     
      // Options array (only care about first element if first element
      // is zero
      metis::idxtype options[5] = {0}; 
      // output argument number of edges cut
      metis::idxtype edgecut = 0;
    
      // Call kmetis
      metis::METIS_PartGraphKway(&(numverts), 
                                 xadj,
                                 adjacency,
                                 vweight,
                                 eweight,
                                 &(weightflag),
                                 &(numflag),
                                 &(nparts),
                                 options,
                                 &(edgecut),
                                 res);
    
      //     // Call pmetis
      //     metis::METIS_PartGraphRecursive(&(numverts), 
      //                                     xadj,
      //                                     adjacency,
      //                                     vweight,
      //                                     eweight,
      //                                     &(weightflag),
      //                                     &(numflag),
      //                                     &(nparts),
      //                                     options,
      //                                     &(edgecut),
      //                                     res);
    
      // destroy all unecessary data structures except res
      if(xadj != NULL) delete [] xadj;
      if(adjacency != NULL) delete [] adjacency;
      if(vweight != NULL) delete [] vweight;
      if(eweight != NULL) delete [] eweight;

      // Resize the partition
      ret_part.resize(num_vertices());
      // process the final results
      ASSERT_NE(res, NULL);
      for(vertex_id_t v = 0; v < num_vertices(); ++v) {
        ret_part[v] = (vertex_id_t)res[v];
      }    
      // Delete the result array
      if(res != NULL) delete [] res;
    } // end of metis partition


    /**
     * This function computes a weighted graph partition using METIS.
     * \param numparts The number of parts to partition into
     * \param[out] ret_part A vector providing a vertex_id -> partition_id mapping
     * \param vfunction A function of the type size_t (*)(const VertexData &v)
     *                  This function takes in the data on a vertex, and 
     *                  returns the weight of the vertex
     * \param wfunction A function of the type size_t (*)(const EdgeData &e)
     *                  This function takes in the data on an edge, and 
     *                  returns the weight of the edge
     * \param usemetisdefaults If set to true, uses Metis default parameter set.
     *                         defaults to false.
     *                         
     *
     * Use a modified version of the METIS library to partition the
     * graph using user provided edge and vertex weight functions.
     * The methis library is described in:
     *
     *   "A Fast and Highly Quality Multilevel Scheme for Partitioning
     *   Irregular Graphs”. George Karypis and Vipin Kumar. SIAM
     *   Journal on Scientific Computing, Vol. 20, No. 1, pp. 359—392,
     *   1999.
     *
     * We have modified an alpha version (5.0) to work with the
     * GraphLab framework.  Therefore users having trouble with this
     * function or the included Metis source should direct concerns to
     * the contact information provided at:
     *
     *   http://www.select.cs.cmu.edu/code
     *
     */
    template <typename EdgeWeightFunction, typename VertexWeightFunction>
    void metis_weighted_partition(size_t numparts ,
                                  std::vector<vertex_id_t>& ret_part,
                                  VertexWeightFunction vfunction,
                                  EdgeWeightFunction wfunction,
                                  bool usemetisdefaults = false) {
      if (numparts == 1) {
        ret_part.assign(num_vertices(), 0);
        return;
      }
      // Determine parameters needed to construct the partitioning
      metis::idxtype numverts(num_vertices());
      ASSERT_GT(numverts, 0);
      // Compute the number of edges 
      metis::idxtype numedges (num_edges());

      // allocate metis data structures
      metis::idxtype* vweight = new metis::idxtype[numverts];
      ASSERT_NE(vweight, NULL);
      metis::idxtype* xadj = new metis::idxtype[numverts + 1];
      ASSERT_NE(xadj, NULL);
      metis::idxtype* adjacency = new metis::idxtype[2 * numedges];
      ASSERT_NE(adjacency, NULL);
      metis::idxtype* eweight = NULL;
      eweight = new metis::idxtype[ 2 * numedges];
      ASSERT_NE(eweight, NULL);
      
      metis::idxtype* res = new metis::idxtype[numverts];   
      ASSERT_NE(res, NULL);

      // Pass through vertices filling in the metis data structures
      size_t offset = 0;
      for(vertex_id_t u = 0; u < num_vertices(); ++u) {
        // Update vertex weight
        // Set weight
        vweight[u] = double(vfunction(vertex_data(u)));
        // Update the offset
        xadj[u] = offset;
        // Fill the the adjacency data
      
        std::set<size_t> neighbors;
        std::map<size_t, double> nbrtoweight;
        foreach(edge_id_t eid, out_edge_ids(u)) {
          neighbors.insert(target(eid));
          nbrtoweight[target(eid)] = double(wfunction(edge_data(eid)));
        }
        foreach(edge_id_t eid, in_edge_ids(u)) {
          neighbors.insert(source(eid));
          nbrtoweight[source(eid)] = double(wfunction(edge_data(eid)));
        }
        foreach(vertex_id_t vid, neighbors) {
          if (vid == u) continue;
          adjacency[offset] = vid;
          eweight[offset] = nbrtoweight[vid];
          ASSERT_GE(adjacency[offset], 0);
          offset++;
          ASSERT_GE(offset, 0);
        }
    
      } // end of data structure creation
      
      // Set the last entry in xadj to the end of the adjacency array
      xadj[numverts] = offset;
    
      // Set additional metis flags
      /**
       * 0 No weights (vwgts and adjwgt are NULL) 
       * 1 Weights on the edges only (vwgts = NULL) 
       * 2 Weights on the vertices only (adjwgt = NULL) 
       * 3 Weights both on vertices and edges. 
       */
      metis::idxtype weightflag = 3;
      // 0 for C-style numbering starting at 0 (1 for fortran style)
      metis::idxtype numflag = 0;
      // the number of parts to cut into
      metis::idxtype nparts = numparts;     
      // Options array (only care about first element if first element
      // is zero
      metis::idxtype options[5] = {0};
      
      options[0] = 1;
      options[1]=3;
      options[2]=1;
      options[3]=2;
      options[4]=0;
      if (usemetisdefaults) options[0] = 0;
      // output argument number of edges cut
      metis::idxtype edgecut = 0;
    
      // Call kmetis
      metis::METIS_PartGraphKway(&(numverts), 
                                 xadj,
                                 adjacency,
                                 vweight,
                                 eweight,
                                 &(weightflag),
                                 &(numflag),
                                 &(nparts),
                                 options,
                                 &(edgecut),
                                 res);
    
      // Call pmetis
      /*   metis::METIS_PartGraphRecursive(&(numverts), 
           xadj,
           adjacency,
           vweight,
           eweight,
           &(weightflag),
           &(numflag),
           &(nparts),
           options,
           &(edgecut),
           res);*/
    
      // destroy all unecessary data structures except res
      if(xadj != NULL) delete [] xadj;
      if(adjacency != NULL) delete [] adjacency;
      if(vweight != NULL) delete [] vweight;
      if(eweight != NULL) delete [] eweight;

      // Resize the partition
      ret_part.resize(num_vertices());
      // process the final results
      ASSERT_NE(res, NULL);
      for(vertex_id_t v = 0; v < num_vertices(); ++v) {
        ret_part[v] = res[v];
      }    
      // Delete the result array
      if(res != NULL) delete [] res;
    } // end of metis partition
    
     /**
     * \brief Performs a breadth first search partitioning of the graph.
     * Equivalent to calling partition() with the 
     * partition_method::PARTITION_EDGE_NUM paramter
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     *
     * The algorithm works by picking up a random vertex and performing a breadth
     * first search until the number of vertices touched is |V|/nparts.
     * This will then be assigned as the first partition. This procedure repeats
     * until all partitions are filled.
     */
    void bfs_partition(size_t nparts, std::vector<vertex_id_t> &vertex2part) {
      // create a list of all unassigned variables
      std::set<vertex_id_t> unassigned;
      vertex2part.resize(num_vertices());
      // initialize the unassigned vertices
      for(vertex_id_t v = 0; v < num_vertices(); ++v) {
        unassigned.insert(v);
      }
      // Compute the partition size
      size_t maxpartsize = (size_t)(std::ceil(double(unassigned.size()) / (double)nparts));
      size_t partid = 0;
      while(!unassigned.empty()) {  
        std::list<vertex_id_t> queue;    // Breadth first queue 
        std::set<vertex_id_t>  visited;  // Set of visited vertices
        // While the task is still too small and their remains
        // unassigned vertices
        size_t curpartsize = 0;
        while(curpartsize < maxpartsize 
              && !unassigned.empty()) {
          if(queue.empty()) { 
            queue.push_front(*unassigned.begin());
            visited.insert(*unassigned.begin());
          }
          ASSERT_FALSE(queue.empty());
          // Pop the first element off the queue 
          vertex_id_t v = queue.front(); queue.pop_front();
          ASSERT_LT(partid, nparts);
          // Add the element to the task
          vertex2part[v] = (uint32_t)partid;
          ++curpartsize;
          // Remove the vertex from the set of unassigned vertices
          unassigned.erase(v); 
          // Add all its unassigned and unvisited neighbors to the queue
          foreach(edge_id_t eid, out_edge_ids(v)) {
            vertex_id_t u = target(eid);
            if(unassigned.find(u) != unassigned.end() &&
               visited.find(u) == visited.end()) {
              queue.push_back(u);
              visited.insert(u);
            }
          }
          foreach(edge_id_t eid, in_edge_ids(v)) {
            vertex_id_t u = source(eid);
            if(unassigned.find(u) != unassigned.end() &&
               visited.find(u) == visited.end()) {
              queue.push_back(u);
              visited.insert(u);
            }
          }
        } // End of block build foor loop
        // move to the next part
        partid++;
      }// end of outer while loop
    } // end of bfs partition



    /**
     * Partition the graph using one of the available partitioning
     * methods.
     * \param partmethod A partition method. \ref partition_method
     *
     * \param nparts The number of parts to partition into
     * \param[out] vertex2part A vector providing a vertex_id -> partition_id mapping
     * 
     */
    void partition(partition_method::partition_method_enum partmethod,
                   size_t nparts, std::vector<uint32_t>& vertex2part) {
      switch (partmethod) {
      case partition_method::PARTITION_METIS:
        return metis_partition(nparts, vertex2part);
      case partition_method::PARTITION_BFS:
        return bfs_partition(nparts, vertex2part);
      case partition_method::PARTITION_RANDOM:
        return random_partition(nparts, vertex2part);
      case partition_method::PARTITION_EDGE_NUM:
        return edge_num_partition(nparts, vertex2part);
      default:
        ASSERT_TRUE(false); //shoud never ever happen
      }
    }

    /**
     * builds a topological_sort of the graph returning it in topsort. 
     * 
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    bool topological_sort(std::vector<vertex_id_t>& topsort) const {
      topsort.clear();
      topsort.reserve(num_vertices());
    
      std::vector<size_t> indeg;
      indeg.resize(num_vertices());
      std::queue<vertex_id_t> q;
      for (size_t i = 0;i < num_vertices(); ++i) {
        indeg[i] = in_edge_ids(i).size();
        if (indeg[i] == 0) {
          q.push(i);
        }
      }
    
      while (!q.empty()) {
        vertex_id_t v = q.front();
        q.pop();
        topsort.push_back(v);
        foreach(edge_id_t eid, out_edge_ids(v)) {
          vertex_id_t destv = target(eid);
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
    class edge {
      vertex_id_t _source;
      vertex_id_t _target;
      EdgeData _data;
    public:
      edge() : _source(-1), _target(-1) { }
      edge(const edge& other) :
        _source(other.source()), _target(other.target()),
        _data(other.data()) { }
      edge(vertex_id_t source, vertex_id_t target) :
        _source(source), _target(target)  { }
      edge(vertex_id_t source, vertex_id_t target, EdgeData data) : 
        _source(source), _target(target), _data(data) {}

      bool operator<(const edge& other) const {
        return (_source < other._source) || 
          (_source == other._source && _target < other._target); 
      }
      
      inline vertex_id_t source() const { return _source; }
      inline vertex_id_t target() const { return _target; }   
      inline EdgeData& data() { return _data; }
      inline const EdgeData& data() const { return _data; }
      
      void load(iarchive& arc) {
        arc >> _source
            >> _target
            >> _data;
      }
      
      void save(oarchive& arc) const {
        arc << _source
            << _target
            << _data;
      }
    }; // end of edge_data


    
    struct edge_id_less_functor {
      graph* g_ptr;
      edge_id_less_functor(graph* g_ptr) : g_ptr(g_ptr) { }
      bool operator()(edge_id_t a, edge_id_t b) {
        return g_ptr->edge_id_less(a,b);
      }
    };
    

    
    /**
     * Used to order edge ids in the in and out edges vectors based on
     * the lexical ordering of the vertex ids of the corresponding
     * edge
     */
    inline bool edge_id_less(edge_id_t a, edge_id_t b) const {
      //      ASSERT_LT(a, edges.size());
      //      ASSERT_LT(b, edges.size());
      return edges[a] < edges[b];
    }

    
 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data */
    std::vector<VertexData> vertices;

    /** The edge data is a vector of edges where each edge stores its
        source, destination, and data. */
    std::vector<edge> edges;
    
    /** A map from src_vertex -> dest_vertex -> edge index */   
    std::vector< std::vector<edge_id_t> >  in_edges;
    
    /** A map from src_vertex -> dest_vertex -> edge index */   
    std::vector< std::vector<edge_id_t> >  out_edges;
    
    /** The vertex colors specified by the user. **/
    std::vector< vertex_color_type > vcolors;  
    
    /** Mark whether the graph is finalized.  Graph finalization is a
        costly procedure but it can also dramatically improve
        performance. */
    bool finalized;
    
    /** increments whenever the graph is cleared. Used to track the
     *  changes to the graph structure  */
    size_t changeid;

    // PRIVATE HELPERS =========================================================>
    /**
     * This function tries to find the edge in the vector.  If it
     * fails it returns size_t(-1)
     * TODO: switch to stl binary search
     */
    size_t binary_search(const std::vector<edge_id_t>& vec,
                         vertex_id_t source, vertex_id_t target) const {
      // Ensure that the graph is finalized before calling this function
      //      finalize();
      ASSERT_TRUE(finalized);
      // Compare to the middle of the list
      size_t first = 0;
      size_t last = vec.size() - 1;
      while(first <= last) {
        size_t mid = (first+last)/2;
        ASSERT_LT(mid, vec.size());
        vertex_id_t mid_source = edges[vec[mid]].source();
        vertex_id_t mid_target = edges[vec[mid]].target();
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
    for(vertex_id_t vid = 0; vid < graph.num_vertices(); ++vid) {
      foreach(edge_id_t eid, graph.out_edge_ids(vid))
        out << vid << ", " << graph.target(eid) << '\n';      
    }
    return out;
  }
  

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

