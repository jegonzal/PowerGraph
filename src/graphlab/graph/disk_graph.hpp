#ifndef GRAPHLAB_DISK_GRAPH_HPP
#define GRAPHLAB_DISK_GRAPH_HPP

#include <omp.h>

#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

  




template<typename VertexData, typename EdgeData> 
class disk_graph {
 public:
  /**
    * Create or open a graph using files
    * fbasename.0 , fbasename.1, etc till fbasename.[numfiles-1]
    * If a graph was created with N files, it must be opened 
    * with numfiles = N or behavior can be unpredictable
    */
  disk_graph(std::string fbasename, size_t numfiles){  
    atoms.resize(numfiles);
    numv.value = 0;
    nume.value = 0;
    for (size_t i = 0;i < numfiles; ++i) {
      atoms[i] = new disk_atom(fbasename + "." + tostr(i), i);
      numv.value += atoms[i]->num_local_vertices();
      nume.value += atoms[i]->num_local_edges();
    }
  }
  
  
  disk_graph<VertexData, EdgeData>&  operator=(const graph<VertexData, EdgeData> &g) {
    clear();
    size_t nv = g.num_vertices();
    #pragma omp parallel for 
    for (int i = 0;i < (int)nv; ++i) {
      vertex_id_t vid = i;
      size_t hashloc = (vid) % atoms.size();
      //place vertices sequentially
      uint16_t owner = (uint16_t)((i * atoms.size())/ nv);
      atoms[hashloc]->set_owner(vid, owner);
      atoms[owner]->add_vertex(vid, owner, g.vertex_data(i));
      atoms[owner]->set_color(vid, g.color(vid));
    }
    
    #pragma omp parallel for 
    for (int i = 0;i < (int)(g.num_edges()); ++i) {
      vertex_id_t target = g.target(i);
      vertex_id_t source = g.source(i);
      uint16_t targetowner = (uint16_t)((target * atoms.size())/ nv);
      uint16_t sourceowner = (uint16_t)((source * atoms.size())/ nv);
      atoms[targetowner]->add_edge(source, target, g.edge_data(i));
      atoms[targetowner]->inc_numlocale();
      // create ghosts
      if (sourceowner != targetowner) {
        atoms[sourceowner]->add_edge(source, target);
        if (atoms[sourceowner]->add_vertex_skip(target, targetowner)) {
          atoms[sourceowner]->set_color(target, g.color(target));
        }
        if (atoms[targetowner]->add_vertex_skip(source, sourceowner)) {
          atoms[targetowner]->set_color(source, g.color(source));
        }
      }
    }
    
    numv.value = g.num_vertices();
    nume.value = g.num_edges();
    
    return *this;
  }
  
  
  ~disk_graph() {
    for (size_t i = 0;i < atoms.size(); ++i) {
      delete atoms[i];
    }
    atoms.clear();
  }
  
  void clear() {
    for (size_t i = 0;i < atoms.size(); ++i) atoms[i]->clear();
    numv.value = 0;
    nume.value = 0;   
  }
  
  void finalize() { 
    // nothing to finalize
  }
  
  /** \brief Get the number of vertices */
  size_t num_vertices() {
    return numv.value;
  }
  
  /** \brief Get the number of vertices local to this machine */
  size_t local_vertices() {
    return numv.value;
  }
  
  /** \brief Get the number of edges */
  size_t num_edges() const {
    return nume.value;
  }


  /** \brief Get the number of in edges of a particular vertex */
  size_t num_in_neighbors(vertex_id_t v) const {
    uint16_t owner = atoms[v % atoms.size()]->get_owner(v);
    ASSERT_NE(owner, (uint16_t)(-1));
    return atoms[owner]->get_in_vertices(v).size();
  }
  
  /** \brief Get the number of out edges of a particular vertex */
  size_t num_out_neighbors(vertex_id_t v) const {
    uint16_t owner = atoms[v % atoms.size()]->get_owner(v);
    ASSERT_NE(owner, (uint16_t)(-1));
    return atoms[owner]->get_in_vertices(v).size();
  }
  

  /** 
  * \brief Creates a vertex containing the vertex data and returns the id
  * of the new vertex id. Vertex ids are assigned in increasing order with
  * the first vertex having id 0.
  * Vertices are placed in atom vid % num_atoms
  */
  vertex_id_t add_vertex(const VertexData& vdata = VertexData()) {
    vertex_id_t v = numv.inc_ret_last();
    uint16_t owner = v % atoms.size();
    atoms[owner]->add_vertex(v, owner, vdata);
    atoms[owner]->set_owner(v, owner);
    return v;
  }
  
  vertex_id_t add_vertex(const VertexData& vdata, 
                         uint16_t locationhint) {
    vertex_id_t v = numv.inc_ret_last();
    uint16_t owner = locationhint;
    atoms[owner]->add_vertex(v, owner, vdata);
    atoms[owner]->set_owner(v, owner);
    return v;
  }
  /**
   * Adds a collection vdata.size() 
   * The collection of vertices will be assigned the name 'collectionname'
   * This name must be unique and cannot be zero length.
   * The return value is the ID of the first vertex inserted.
   * The remaining vertices have the consecutive vertex IDs.
   */
  vertex_id_t add_vertex_collection(std::string collectionname,
                                    std::vector<VertexData> & vdata) {
    vertex_id_t ret = numv.inc_ret_last(vdata.size());
    collection_range[collectionname] = std::make_pair(ret, ret + vdata.size());

    #pragma omp parallel for
    for (int i = 0;i < (int)(vdata.size()); ++i) {
      vertex_id_t vid = i + ret;
      size_t hashloc = (vid) % atoms.size();
      //place vertices sequentially
      uint16_t owner = (uint16_t)((i * atoms.size())/ vdata.size());
      atoms[hashloc]->set_owner(vid, owner);
      atoms[owner]->add_vertex(vid, owner, vdata[i]);
    }
    return numv;
  }

  vertex_id_t add_vertex_collection(std::string collectionname,
                                    size_t numvertices,
                                    boost::function<VertexData (vertex_id_t)> generator) {
    vertex_id_t ret = numv.inc_ret_last(numvertices);
    collection_range[collectionname] = std::make_pair(ret, ret + numvertices);
    #pragma omp parallel for
    for (int i = 0;i < (int)numvertices; ++i) {
      vertex_id_t vid = i + ret;
      size_t hashloc = (vid) % atoms.size();
      //place vertices sequentially
      uint16_t owner = (uint16_t)((i * atoms.size())/ numvertices);
      atoms[hashloc]->set_owner(vid, owner);
      atoms[owner]->add_vertex(vid, owner, generator((vertex_id_t)i));
    }
    return ret;
  }

  
  
  /**
  * \brief Creates an edge connecting vertex source to vertex target. 
  *        Any existing data will be cleared.
  */
  void add_edge(vertex_id_t source, vertex_id_t target, 
                       const EdgeData& edata = EdgeData()) {
    nume.inc();
    uint16_t targetowner = atoms[target % atoms.size()]->get_owner(target);
    uint16_t sourceowner = atoms[source % atoms.size()]->get_owner(source);
    ASSERT_NE(targetowner, (uint16_t)(-1));
    ASSERT_NE(sourceowner, (uint16_t)(-1));
    atoms[targetowner]->add_edge(source, target, edata);
    atoms[targetowner]->inc_numlocale();
    // create ghosts
    if (sourceowner != targetowner) {
      atoms[sourceowner]->add_edge(source, target);
      atoms[sourceowner]->add_vertex_skip(target, targetowner);
      atoms[targetowner]->add_vertex_skip(source, sourceowner);
    }
  }


  
  /**
   * Defines the adjacency structure for a collection of vertices in
   * "collectionname". If collectionname is empty, it refers to all vertices
   * in the graph.
   * 
   * The provided user function efn will be called on each vertex in the collection,
   * passing the user function the data on the vertex (vdata) as well as the
   * ID of the vertex.
   * 
   * The function should then return in the remaining arguments, the 
   * set of outgoing vertices and their corresponding edge data, as well as the
   * set of incoming vertices and their corresponding edge data.
   * The prototype for the function is:
   * void efn(vertex_id_t vid, \n 
   *         const VertexData& vdata, \n
   *         std::vector<vertex_id_t>& inv,  \n
   *         std::vector<EdgeData>& inedata, \n
   *         std::vector<vertex_id_t>& outv,  \n
   *         std::vector<EdgeData>& outedata) \n
   */
  void add_edge_indirect(std::string collectionname,
                         boost::function<void (vertex_id_t,
                                               const VertexData&,
                                               std::vector<vertex_id_t>&, 
                                               std::vector<EdgeData>&,
                                               std::vector<vertex_id_t>&, 
                                               std::vector<EdgeData>&)> efn) {
    std::pair<vertex_id_t, vertex_id_t> vidrange = collection_range[collectionname];
    
    #pragma omp parallel for 
    for (vertex_id_t vid = vidrange.first; vid < vidrange.second; ++vid) {
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      VertexData vdata;
      ASSERT_TRUE(atoms[owner]->get_vertex(vid, owner, vdata));
      std::vector<vertex_id_t> inv;
      std::vector<EdgeData> inedata;
      std::vector<vertex_id_t> outv; 
      std::vector<EdgeData> outedata;
      efn(vid, vdata, inv, inedata, outv, outedata);
      for (size_t i = 0;i < inv.size(); ++i) add_edge(inv[i], vid, inedata[i]);
      for (size_t i = 0;i < outv.size(); ++i) add_edge(vid, outv[i], outedata[i]);
    }
  }
  
  
  void set_color(vertex_id_t vid, vertex_color_type color) {
    uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
    ASSERT_NE(owner, (uint16_t)(-1));
    atoms[owner]->set_color(vid, color);
  }

  /** \brief Returns the vertex color of a vertex.
        Coloring is only valid if compute_coloring() is called first.*/
  vertex_color_type get_color(vertex_id_t vid) {
    uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
    ASSERT_NE(owner, (uint16_t)(-1));
    return atoms[owner]->get_color(vid);
  }


  /** \brief This function constructs a heuristic coloring for the 
    graph and returns the number of colors */
  size_t compute_coloring() {
    // Reset the colors
    #pragma omp parallel for
    for(vertex_id_t v = 0; v < num_vertices(); ++v) set_color(v, 0);
    
    // Recolor
    size_t max_color = 0;
    std::set<vertex_color_type> neighbor_colors;
    for(vertex_id_t vid = 0; vid < num_vertices(); ++vid) {
      neighbor_colors.clear();
      // Get the neighbor colors
      foreach(vertex_id_t neighbor_vid, in_vertices(vid)){
        const vertex_color_type& neighbor_color = get_color(neighbor_vid);
        neighbor_colors.insert(neighbor_color);
      }
      foreach(vertex_id_t neighbor_vid, out_vertices(vid)){
        const vertex_color_type& neighbor_color = get_color(neighbor_vid);
        neighbor_colors.insert(neighbor_color);
      }

      vertex_color_type vertex_color = 0;
      foreach(vertex_color_type neighbor_color, neighbor_colors) {
        if(vertex_color != neighbor_color) break;
        else vertex_color++;
        // Ensure no wrap around
        ASSERT_NE(vertex_color, 0);                
      }
      set_color(vid, vertex_color);
      max_color = std::max(max_color, size_t(vertex_color) );
    }
    // Return the NUMBER of colors
    propagate_coloring();
    return max_color + 1;     
  }
  
  std::vector<vertex_id_t> in_vertices(vertex_id_t vid) const {
    uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
    ASSERT_NE(owner, (uint16_t)(-1));
    return atoms[owner]->get_in_vertices(vid);
  }
  
  std::vector<vertex_id_t> out_vertices(vertex_id_t vid) const {
    uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
    ASSERT_NE(owner, (uint16_t)(-1));
    return atoms[owner]->get_out_vertices(vid);
  }
  
  /** \brief Gets the edge data of the edge from vertex 'src' to vertex 'dest'*/
  EdgeData get_edge_data(vertex_id_t src, vertex_id_t dest) const {
    EdgeData ret;
    uint16_t owner = atoms[dest % atoms.size()]->get_owner(dest);
    ASSERT_NE(owner, (uint16_t)(-1));
    ASSERT_TRUE(atoms[owner]->get_edge(src, dest, ret));
    return ret;
  }
  
  /** \brief Sets the edge data of the edge from vertex 'src' to vertex 'dest'*/
  void set_edge_data(vertex_id_t src, vertex_id_t dest, const EdgeData &edata) {
    uint16_t owner = atoms[dest % atoms.size()]->get_owner(dest);
    ASSERT_NE(owner, (uint16_t)(-1));
    atoms[owner]->set_edge(src, dest, edata);
  }
  
  /** \brief Gets the vertex data of the vertex with id 'vid' */
  VertexData get_vertex_data(vertex_id_t vid) {
    VertexData ret;
    uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
    ASSERT_NE(owner, (uint16_t)(-1));
    ASSERT_TRUE(atoms[owner]->get_vertex(vid, owner, ret));
    return ret;
  }
  /** \brief Sets the vertex data of the vertex with id 'vid' */  
  void set_vertex_data(vertex_id_t vid, const VertexData& vdata) {
    uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
    ASSERT_NE(owner, (uint16_t)(-1));
    atoms[owner]->set_vertex(vid, owner, vdata);
  }
  
 private:
   
  // block copies
  disk_graph(const disk_graph&);
  disk_graph& operator=(const disk_graph&);
    
  mutable std::vector<disk_atom*> atoms;
  atomic<size_t> numv;
  atomic<size_t> nume;
  
  std::map<std::string, std::pair<vertex_id_t, vertex_id_t> > collection_range;
  
  void propagate_coloring() {
    #pragma omp parallel for
    for (size_t i = 0;i < atoms.size(); ++i) {
      foreach(vertex_id_t vid, atoms[i]->enumerate_vertices()) {
        uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
        if (owner != i) {
          atoms[i]->set_color(vid, atoms[owner]->get_color(vid));
        }
      }
    }
  }
};
  

  
} // namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif