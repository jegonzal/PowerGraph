#ifndef GRAPH_LOCAL_STORE_HPP
#define GRAPH_LOCAL_STORE_HPP
#include <climits>
#include <graphlab/util/mmap_wrapper.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {
namespace dist_graph_impl {
  

#define PREFETCH_MERGE_LIMIT 4096
  
  template<typename VertexData, typename EdgeData> class graph_local_store;


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
  class graph_local_store {
  public:

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;

    struct vdata_store {
      VertexData data;
      uint64_t version;
    };

    struct edata_store {
      EdgeData data;
      uint64_t version;
    };
  public:

    /**
     * Build a basic graph
     */
    graph_local_store(): vertices(NULL), edgedata(NULL), finalized(true), changeid(0) {  }

    void create_store(size_t create_num_verts, size_t create_num_edges,
                std::string vertexstorefile, std::string edgestorefile) { 
      nvertices = create_num_verts;
      nedges = create_num_edges;
      
      edges.resize(nedges);
      in_edges.resize(nvertices);
      out_edges.resize(nvertices);
      vcolors.resize(nvertices);
      
      vertex_store_file = vertexstorefile;
      edge_store_file = edgestorefile;
      

      finalized = true;
      changeid = 0;
      
      // allocate the vdata and edata
      setup_mmap();
      
    }

    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
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
      return nvertices;
    } // end of num vertices

    /** \brief Get the number of edges */
    size_t num_edges() const {
      return nedges;
    } // end of num edges


    /** \brief Get the number of in edges of a particular vertex */
    size_t num_in_neighbors(vertex_id_t v) const {
      assert(v < nvertices);
      return in_edges[v].size();
    } // end of num vertices
    
    /** \brief Get the number of out edges of a particular vertex */
    size_t num_out_neighbors(vertex_id_t v) const  {
      assert(v < nvertices);
      return out_edges[v].size();
    } // end of num vertices

    /** \brief Finds an edge.
    The value of the first element of the pair will be true if an 
    edge from src to target is found and false otherwise. If the 
    edge is found, the edge ID is returned in the second element of the pair. */
    std::pair<bool, edge_id_t>
    find(vertex_id_t source, vertex_id_t target) const {


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
            assert(eid < nedges);
            if(edges[eid].source() == source 
               && edges[eid].target() == target) {
              return std::make_pair(true, eid);
            }
          }
          return std::make_pair(false, -1);
        } else { // fewer out edges at the source
          // linear search the out_edges at the source
          foreach(edge_id_t eid, out_edges[source]) {
            assert(eid < nedges);
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
      assert(res.second < nedges);
      return res.second;
    } // end of edge_id

    
    /** \brief Returns the edge ID of the edge going in the opposite direction. 
        Assertion failure if such an edge is not found.  */
    edge_id_t rev_edge_id(edge_id_t eid) const {
      assert(eid < nedges);
      vertex_id_t source = edges[eid].source();
      vertex_id_t target = edges[eid].target();    
      return edge_id(target, source);
    } // end of rev_edge_id

  
    
    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    void add_edge(edge_id_t edge_id, vertex_id_t source, vertex_id_t target) {
      if ( source >= nvertices
           || target >= nvertices ) {

        logstream(LOG_FATAL) 
          << "Attempting add_edge (" << source
          << " -> " << target
          << ") when there are only " << nvertices
          << " vertices" << std::endl;

        ASSERT_MSG(source < nvertices, "Invalid source vertex!");
        ASSERT_MSG(target < nvertices, "Invalid target vertex!");
      }

      if (edge_id >= nedges) {
        ASSERT_MSG(edge_id < nedges, "Invalid edge ID!");
      }
      if(source == target) {
        logstream(LOG_FATAL) 
          << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
          << "This operation is not permitted in GraphLab!" << std::endl;
        ASSERT_MSG(source != target, "Attempting to add self edge!");
      }

      // Add the edge to the set of edge data (this copies the edata)
      edges[edge_id] = edge(source, target);
      

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
    } // End of add edge
        
    
    /** \brief Returns a reference to the data stored on the vertex v. */
    VertexData& vertex_data(vertex_id_t v) {
      assert(v < nvertices);
      return vertices[v].data;
    } // end of data(v)

    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_t v) const {
      assert(v < nvertices);
      return vertices[v].data;
    } // end of data(v)

    /// Sets the vertex version. Setting the version also clears the modified flag
    void set_vertex_version(vertex_id_t v, uint64_t version) {
      vertices[v].version = version;
    }

    void increment_vertex_version(vertex_id_t v) {
      ++vertices[v].version;
    }

    uint64_t vertex_version(vertex_id_t v) const{
      return vertices[v].version & uint64_t(0x7FFFFFFFFFFFFFFF);
    }


    void set_vertex_modified(vertex_id_t v, bool modified) {
      if (modified) {
        vertices[v].version |= uint64_t(0x8000000000000000);
      }
      else {
        vertices[v].version = vertex_version(v);
      }
    }

    bool vertex_modified(vertex_id_t v) const{
      return vertices[v].version & uint64_t(0x8000000000000000);
    }

    
    /** \brief Returns a reference to the data stored on the edge source->target. */
    EdgeData& edge_data(vertex_id_t source, vertex_id_t target) {
      assert(source < nvertices);
      assert(target < nvertices);
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      assert(ans.first);
      // the edge id should be valid!
      assert(ans.second < nedges);
      return edgedata[ans.second].data;
    } // end of edge_data(u,v)
    
    /** \brief Returns a constant reference to the data stored on the edge source->target */
    const EdgeData& edge_data(vertex_id_t source, vertex_id_t target) const {
      assert(source < nvertices);
      assert(target < nvertices);
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      assert(ans.first);
      // the edge id should be valid!
      assert(ans.second < nedges);
      return edgedata[ans.second].data;
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(edge_id_t edge_id) { 
      assert(edge_id < nedges);
      return edgedata[edge_id].data;
    }
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(edge_id_t edge_id) const {
      assert(edge_id < nedges);
      return edgedata[edge_id].data;
    }

    void set_edge_version(edge_id_t edge_id, uint64_t version) {
      edgedata[edge_id].version = version;
    }

    void increment_edge_version(edge_id_t edge_id) {
      ++edgedata[edge_id].version;
    }
    
    uint64_t edge_version(edge_id_t edge_id) const{
      return edgedata[edge_id].version & uint64_t(0x7FFFFFFFFFFFFFFF);
    }

    void set_edge_modified(edge_id_t edge_id, bool modified) {
      if (modified) {
        edgedata[edge_id].version |= uint64_t(0x8000000000000000);
      }
      else {
        edgedata[edge_id].version = edge_version(edge_id);
      }
    }

    bool edge_modified(vertex_id_t v) const{
      return edgedata[v].version & uint64_t(0x8000000000000000);
    }

    size_t& edge_version(vertex_id_t source, vertex_id_t target) {
      assert(source < nvertices);
      assert(target < nvertices);
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      assert(ans.first);
      return edgedata[ans.second].version;
    } // end of edge_data(u,v)

    void increment_edge_version(vertex_id_t source, vertex_id_t target) {
      assert(source < nvertices);
      assert(target < nvertices);
      std::pair<bool, edge_id_t> ans = find(source, target);
      assert(ans.first);
      ++edgedata[ans.second].version;
    }

    size_t edge_version(vertex_id_t source, vertex_id_t target) const {
      assert(source < nvertices);
      assert(target < nvertices);
      std::pair<bool, edge_id_t> ans = find(source, target);
      // We must find the edge!
      assert(ans.first);
      return edgedata[ans.second].version;
    } // end of edge_data(u,v)

    /** \brief Returns the source vertex of an edge. */
    vertex_id_t source(edge_id_t edge_id) const {
      //      assert(edge_id < nedges);
      return edges[edge_id].source();
    }

    /** \brief Returns the destination vertex of an edge. */
    vertex_id_t target(edge_id_t edge_id) const {
      //      assert(edge_id < nedges);
      return edges[edge_id].target();    
    }
    
    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_t vertex) const {
      assert(vertex < nvertices);
      return vcolors[vertex];
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type& color(vertex_id_t vertex) {
      assert(vertex < nvertices);
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

      return edge_list(in_edges[v]);
    } // end of in edges    

    /** \brief Return the edge ids of the edges leaving at v */
    edge_list out_edge_ids(vertex_id_t v) const {

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
      arc >> nvertices
          >> nedges
          >> edges
          >> in_edges
          >> out_edges
          >> vcolors
          >> finalized;
      // rebuild the map
      delete vertexmmap;
      delete edgemmap;
      setup_mmap();
      deserialize(arc, vertices, sizeof(VertexData) * nvertices);
      deserialize(arc, edgedata, sizeof(EdgeData) * nedges);
      
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      // Write the number of edges and vertices
      arc << nvertices
          << nedges
          << edges
          << in_edges
          << out_edges
          << vcolors
          << finalized;

      serialize(arc, vertices, sizeof(VertexData) * nvertices);
      serialize(arc, edgedata, sizeof(EdgeData) * nedges);
          
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
      for(size_t i = 0; i < nedges; ++i) {
        fout << edges[i].source() << ", " << edges[i].target() << "\n";
        assert(fout.good());
      }          
      fout.close();
    }
    
    void flush() {
      vertexmmap->sync_all();
      edgemmap->sync_all();
    }
    
    void background_flush() {
      vertexmmap->background_sync_all();
      edgemmap->background_sync_all();
    }
    
    
    
    void compute_minimal_prefetch() {
      minimal_prefetch_vertex.resize(nvertices);
      minimal_prefetch_edge.resize(nvertices);
      for (size_t v = 0;v < nvertices; ++v) {
        std::map<void*, size_t> prefetchvertex;
        std::map<void*, size_t> prefetchedge;
        // first get a list of all the prefetch targets
        prefetchvertex[vertices + v] = sizeof(VertexData);
        for (size_t i = 0;i < in_edges[v].size(); ++i) {
          prefetchvertex[vertices + edges[in_edges[v][i]].source()] = sizeof(VertexData);
          prefetchedge[edgedata + in_edges[v][i]] = sizeof(EdgeData);
        }
        for (size_t i = 0;i < out_edges[v].size(); ++i) {
          prefetchvertex[vertices + edges[out_edges[v][i]].target()] = sizeof(VertexData);
          prefetchedge[edgedata + out_edges[v][i]] = sizeof(EdgeData);
        }
        reduce_prefetch_list(prefetchvertex);
        reduce_prefetch_list(prefetchedge);
        
        minimal_prefetch_vertex[v].clear();
        std::copy(prefetchvertex.begin(), prefetchvertex.end(),
                  std::back_inserter(minimal_prefetch_vertex[v]));
                  
        minimal_prefetch_edge[v].clear(); 
        std::copy(prefetchedge.begin(), prefetchedge.end(),
                  std::back_inserter(minimal_prefetch_edge[v]));
      }
    }
    void print_prefetch_list(vertex_id_t v) {
      std::cout << "Vertex " << v << " prefetch list:\n";
      std::cout << "Vertex: \n";
      for (size_t i = 0;i < minimal_prefetch_vertex[v].size(); ++i) {
        std::cout << "at idx: " << ((char*)minimal_prefetch_vertex[v][i].first - (char*)vertices) / sizeof(VertexData) 
                  << " span "
                  << minimal_prefetch_vertex[v][i].second / sizeof(VertexData) << "\n";
      }
      std::cout << "Edge: \n";
      for (size_t i = 0;i < minimal_prefetch_edge[v].size(); ++i) {
        std::cout << "at idx: " << ((char*)minimal_prefetch_edge[v][i].first - (char*)edgedata) / sizeof(EdgeData) 
                  << " span "
                  << minimal_prefetch_edge[v][i].second / sizeof(EdgeData) << "\n";      }
    }
    void prefetch_scope(vertex_id_t v) {
      for (size_t i = 0;i < minimal_prefetch_vertex[v].size(); ++i) {
        vertexmmap->prefetch(minimal_prefetch_vertex[v][i].first, minimal_prefetch_vertex[v][i].second);
      }
      for (size_t i = 0;i < minimal_prefetch_edge[v].size(); ++i) {
        edgemmap->prefetch(minimal_prefetch_edge[v][i].first, minimal_prefetch_edge[v][i].second);
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
      graph_local_store<VertexData, EdgeData> * g_ptr;
      edge_id_less_functor(graph_local_store<VertexData, EdgeData>* g_ptr) : g_ptr(g_ptr) { }
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
      assert(a < nedges);
      assert(b < nedges);
      return edges[a] < edges[b];
    }

    
 
    // PRIVATE DATA MEMBERS ===================================================>    
    /** The vertex data is simply a vector of vertex data 
     */
    vdata_store* vertices;
    
    /** Vector of edge data  */
    edata_store* edgedata;
    
    std::string vertex_store_file;
    std::string edge_store_file;
    
    mmap_wrapper *vertexmmap, *edgemmap;
    
    /** The edge data is a vector of edges where each edge stores its
        source, destination. */
    std::vector<edge> edges;
    
    /** A map from src_vertex -> dest_vertex -> edge index */   
    std::vector< std::vector<edge_id_t> >  in_edges;
    
    /** A map from src_vertex -> dest_vertex -> edge index */   
    std::vector< std::vector<edge_id_t> >  out_edges;
    
    /** The vertex colors specified by the user. **/
    std::vector< vertex_color_type > vcolors;  
    
    std::vector<std::vector<std::pair<void*, size_t> > > minimal_prefetch_vertex; 
    std::vector<std::vector<std::pair<void*, size_t> > > minimal_prefetch_edge;
    size_t nvertices;
    size_t nedges;
    
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

    void setup_mmap() {
      vertexmmap = new mmap_wrapper(vertex_store_file, sizeof(vdata_store) * nvertices);
      edgemmap = new mmap_wrapper(edge_store_file, sizeof(edata_store) * nedges);
      vertices = (vdata_store*)(vertexmmap->mapped_ptr());
      edgedata = (edata_store*)(edgemmap->mapped_ptr());
    }


    std::pair<void*, size_t> merge_targets(std::pair<void*, size_t> lower,
                                           std::pair<void*, size_t> higher) {
      char* lowleftptr = (char*)lower.first;
      char* lowrightptr = (char*)lower.first + lower.second;
      char* highleftptr = (char*)higher.first;
      char* highrightptr = (char*)higher.first + higher.second;
      if (lowrightptr >= highleftptr && lowrightptr >= highrightptr) {
        // new target intersects an existing target
        return lower;
      }
      else if (lowrightptr + PREFETCH_MERGE_LIMIT >= highrightptr && 
            lowrightptr + PREFETCH_MERGE_LIMIT >= highrightptr) {
        // see if we can extend the old target to include the new target.
        // don't extend by more than PREFETCH_MERGE_LIMIT
        // yes we do! extend
        lower.second = (char*)highrightptr - (char*)lowleftptr;
        return lower;
      }
      else {
        return std::make_pair<void*, size_t>(NULL, 0);
      }
    }
     
    void reduce_prefetch_list(std::map<void*, size_t> &current) {
      std::map<void*, size_t>::iterator iter = current.begin();
      while(iter != current.end()) {
        while(1) {
          std::map<void*, size_t>::iterator next = iter;
          next++;
          if (next == current.end()) break;    
          std::pair<void*, size_t> ret;
          ret = merge_targets(std::make_pair<void*, size_t>(iter->first, iter->second),
                              std::make_pair<void*, size_t>(next->first, next->second));
          
          if (ret.first != NULL) {
            iter->second = ret.second;
            current.erase(next);
          }
          else {
            break;
          }
        }
        ++iter;
      }
    }

   

  }; // End of graph

  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const graph_local_store<VertexData, EdgeData>& graph) {
  
    for(vertex_id_t vid = 0; vid < graph.num_vertices(); ++vid) {
      foreach(edge_id_t eid, graph.out_edge_ids(vid))
        out << vid << ", " << graph.target(eid) << '\n';      
    }
    return out;
  }
  
  
  
}
}
#include <graphlab/macros_undef.hpp>

#endif
