#ifndef GRAPH_LOCAL_STORE_HPP
#define GRAPH_LOCAL_STORE_HPP
#include <src/graphlab/graph/graph.hpp>
#include <src/graphlab/logger/assertions.hpp>
namespace graphlab {
namespace dist_graph_impl {
  


  
  template<typename VertexData, typename EdgeData> class graph;


  // CLASS GRAPH ==============================================================>
  /**
    \brief The local storage abstraction for the distributed graph data type.

   This implements local storage for a graph for a distributed graph data type
   and and is not meant to be used directly. It essentially replicates a 
   simplified version of the \ref graph datatype, but is
   modified to store vertex and edge data on disk.
   
   The local graph store only manages "local" vertex and edge ids and does
   not provide local <-> global mappings. This must be done at a higher
   level container */
  template<typename VertexData, typename EdgeData>
  class local_graph_store {
  public:

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
      
    
  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    local_graph_store(size_t maxv, size_t maxe) : finalized(true),changeid(0) {  }

    /**
     * \BUG: Should not reserve but instead directly create vertices.
     * Create a graph with space allocated for num_vertices and
     * num_edges.
     */
    local_graph_store(size_t maxv, size_t maxe, size_t nverts) : 
      vertices(nverts),
      in_edges(nverts), out_edges(nverts), vcolors(nverts),
      finalized(true),changeid(0) { }


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
      // check to see if the graph is already finalized
      if(finalized) return;

      // Assert that the graph is not finalized 
      typedef std::vector< edge_id_t > edge_set;

      edge_id_less_functor less_functor(this);
      
      // Sort all in edges sets
      foreach(edge_set& eset, in_edges) {
        std::sort(eset.begin(),
                  eset.end(),
                  less_functor);
        //                  boost::bind(&graph::edge_id_less, this, _1, _2) );
        // Test for duplicate edges
        if (eset.size() > 1) {
          for(size_t i = 0; i < eset.size()-1; ++i) {
            // Duplicate edge test
            assert(edge_id_less(eset[i], eset[i+1]));
          }
        }        
      }
      // Sort all out edges sets
      foreach(edge_set& eset, out_edges) {
        std::sort(eset.begin(),
                  eset.end(),
                  less_functor);
        //                 boost::bind(&graph::edge_id_less, this, _1, _2) );
        // Test for dupliate edges
        if (eset.size() > 1) {
          for(size_t i = 0; i < eset.size()-1; ++i) {
            // Duplicate edge test
            assert(edge_id_less(eset[i], eset[i+1]));
          }
        }
      }
      finalized = true;
    } // End of finalize
            
    /** \brief Get the number of vetices */
    size_t num_vertices() const {
      return vertices.size();
    } // end of num vertices

    /** \brief Get the number of edges */
    size_t num_edges() const {
      return edges.size();
    } // end of num edges


    /** \brief Get the number of in edges of a particular vertex */
    size_t num_in_neighbors(vertex_id_t v) const {
      assert(v < vertices.size());
      return in_edges[v].size();
    } // end of num vertices
    
    /** \brief Get the number of out edges of a particular vertex */
    size_t num_out_neighbors(vertex_id_t v) const  {
      assert(v < vertices.size());
      return out_edges[v].size();
    } // end of num vertices

    /** \brief Finds an edge.
    The value of the first element of the pair will be true if an 
    edge from src to target is found and false otherwise. If the 
    edge is found, the edge ID is returned in the second element of the pair. */
    std::pair<bool, edge_id_t>
    find(vertex_id_t source, vertex_id_t target) const {
      assert(source < in_edges.size());
      assert(target < out_edges.size());
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
            assert(eid < edges.size());
            if(edges[eid].source() == source 
               && edges[eid].target() == target) {
              return std::make_pair(true, eid);
            }
          }
          return std::make_pair(false, -1);
        } else { // fewer out edges at the source
          // linear search the out_edges at the source
          foreach(edge_id_t eid, out_edges[source]) {
            assert(eid < edges.size());
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
      assert(res.first);
      assert(res.second < edges.size());
      return res.second;
    } // end of edge_id

    
    /** \brief Returns the edge ID of the edge going in the opposite direction. 
        Assertion failure if such an edge is not found.  */
    edge_id_t rev_edge_id(edge_id_t eid) const {
      assert(eid < edges.size());
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
      return vertices.size() - 1;
    } // End of add vertex;


    /** 
     * \brief Add additional vertices up to provided num_vertices.  This will
     * fail if resizing down.
     */
    void resize(size_t num_vertices ) {
      assert(num_vertices >= vertices.size());
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
      edges.push_back( edge( source, target) );
      
      /** \todo: switch to mmap*/
      edgedata.push_back(edata);

      // Add the edge id to in and out edge maps
      edge_id_t edge_id = edges.size() - 1;
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
      assert(v < vertices.size());
      return vertices[v];
    } // end of data(v)
    
    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_t v) const {
      assert(v < vertices.size());
      return vertices[v];
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_t source, vertex_id_t target) {
      assert(source < vertices.size());
      assert(target < vertices.size());
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      assert(ans.first);
      // the edge id should be valid!
      assert(ans.second < edges.size());
      return edgedata[ans.second];
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the edge source->target */
    const EdgeData& edge_data(vertex_id_t source, vertex_id_t target) const {
      assert(source < vertices.size());
      assert(target < vertices.size());
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      assert(ans.first);
      // the edge id should be valid!
      assert(ans.second < edges.size());
      return edgedata[ans.second];
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_id_t edge_id) { 
      assert(edge_id < edges.size());
      return edgedata[edge_id];
    }
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(edge_id_t edge_id) const {
      assert(edge_id < edges.size());
      return edgedata[edge_id];
    }

    /** \brief Returns the source vertex of an edge. */
    vertex_id_t source(edge_id_t edge_id) const {
      //      assert(edge_id < edges.size());
      return edges[edge_id].source();
    }

    /** \brief Returns the destination vertex of an edge. */
    vertex_id_t target(edge_id_t edge_id) const {
      //      assert(edge_id < edges.size());
      return edges[edge_id].target();    
    }
    
    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_t vertex) const {
      assert(vertex < vertices.size());
      return vcolors[vertex];
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_t vertex) {
      assert(vertex < vertices.size());
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
        edge_list in_edges = in_edge_ids(vid);
        // Get the neighbor colors
        foreach(edge_id_t eid, in_edges){
          const vertex_id_t& neighbor_vid = source(eid);
          const vertex_color_type& neighbor_color = color(neighbor_vid);
          neighbor_colors.insert(neighbor_color);
        }
        // Find the lowest free color
        vertex_color_type& vertex_color = color(vid);
        foreach(vertex_color_type neighbor_color, neighbor_colors) {
          if(vertex_color != neighbor_color) break;
          else vertex_color++;
          // Ensure no wrap around
          assert(vertex_color != 0);                
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
    bool valid_coloring() {
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
      assert(v < in_edges.size());
      return edge_list(in_edges[v]);
    } // end of in edges    

    /** \brief Return the edge ids of the edges leaving at v */
    edge_list out_edge_ids(vertex_id_t v) const {
      assert(v < out_edges.size());
      return edge_list(out_edges[v]);
    } // end of out edges
    
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
      assert(fout.good());
      for(size_t i = 0; i < edges.size(); ++i) {
        fout << edges[i].source() << ", " << edges[i].target() << "\n";
        assert(fout.good());
      }          
      fout.close();
    }
    
    /**
     * Renumbers vertex id revmap[i] to vertex_id i.
     * revmap must be #vertices in length, and the range of revmap
     * must be [0, #vertices-1] where no element is duplicated.
     * The revmap vector will be modified by this function call
     */
    void renumber_vids(std::vector<vertex_id_t> &revmap) {
      // quick sanity check
      ASSERT_EQ(revmap.size(), vertices.size());

      // build the forward map
      std::vector<vertex_id_t> forwardmap(revmap.size());      
      for (size_t i = 0;i < revmap.size(); ++i) {
        forwardmap[revmap[i]] = i;
      }
      
      // forward map all the edges
      for (size_t i = 0; i < edges.size(); ++i) {
        edges[i]._source = forwardmap[edges[i]._source];
        edges[i]._target = forwardmap[edges[i]._target];
      }
      // remap the vertex data by swapping around the cycles
      for (size_t i = 0;i < revmap.size(); ++i) {
        // check if I need to remap
        if (remap[i] != i) {
          // yes I do!
          // begin a remapping cycle
          // remember the value of the first element in the cycle
          VertexData initialdata = vdata[i];
          size_t j = i;
          size_t prev;
          while (1) {
            if (remap[j] != i) {
              // if we are not back to the start of the cycle
              // move the data on cycle upwards the cycle
              vdata[j] = vdata[remap[j]];
              prev = j;
              // move down the cycle
              j = remap[j];
              // reset the remap value on this entry
              remap[prev] = prev;
            }
            else {
              // back at the start of the cycle!
              // make a sanity check
              vdata[j] = initialdata;
              // and we are done!
              remap[j] = j;
              break;
            }
          }
        }
      }
    }
    
  private:    
    /** Internal edge class  */   
    struct edge {
      vertex_id_t _source;
      vertex_id_t _target;

      edge() : _source(-1), _target(-1) { }
      edge(const edge& other) :
        _source(other.source()), _target(other.target()) { }
      edge(vertex_id_t source, vertex_id_t target) :
        _source(source), _target(target)  { }
      edge(vertex_id_t source, vertex_id_t target) : 
        _source(source), _target(target) {}

      bool operator<(const edge& other) const {
        return (_source < other._source) || 
          (_source == other._source && _target < other._target); 
      }
      
      inline vertex_id_t source() const { return _source; }
      inline vertex_id_t target() const { return _target; }   
      
      void load(iarchive& arc) {
        arc >> _source
            >> _target;
      }
      
      void save(oarchive& arc) const {
        arc << _source
            << _target;
      }
    }; // end of edge


    
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
      assert(a < edges.size());
      assert(b < edges.size());
      return edges[a] < edges[b];
    }

    
 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data 
     * \todo To switch to mmap
     */
    std::vector<VertexData> vertices;
    
    /** Vector of edge data 
     * \todo To switch to mmap*/
    std::vector<EdgeData> edgedata;
    
    /** The edge data is a vector of edges where each edge stores its
        source, destination. */
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
     */
    size_t binary_search(const std::vector<edge_id_t>& vec,
                         vertex_id_t source, vertex_id_t target) const {
      // Ensure that the graph is finalized before calling this function
      //      finalize();
      assert(finalized);
      // Compare to the middle of the list
      size_t first = 0;
      size_t last = vec.size() - 1;
      while(first <= last) {
        size_t mid = (first+last)/2;
        assert(mid < vec.size());
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
          assert(mid > 0);
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
                           const local_graph_store<VertexData, EdgeData>& graph) {
    for(vertex_id_t vid = 0; vid < graph.num_vertices(); ++vid) {
      foreach(edge_id_t eid, graph.out_edge_ids(vid))
        out << vid << ", " << graph.target(eid) << '\n';      
    }
    return out;
  }
  
  
  
}
}

#endif
