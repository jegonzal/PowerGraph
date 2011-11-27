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



#ifndef GRAPHLAB_DISK_GRAPH_HPP
#define GRAPHLAB_DISK_GRAPH_HPP

#include <omp.h>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/graph_partitioner.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <graphlab/graph/memory_atom.hpp>
#include <graphlab/graph/write_only_disk_atom.hpp>
#include <graphlab/graph/atom_index_file.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

  
  struct disk_graph_atom_type {
    enum atom_type {
      DISK_ATOM,
      MEMORY_ATOM,
      WRITE_ONLY_ATOM
    };
  };



  /**
     \brief GraphLab disk based Graph container templatized over the vertex and edge types.
  
     This class more or less exposes the same interface as graphlab::graph
     but with a few critical differences.
  
     Firstly, the user will notice that edge IDs are no longer assigned. 
     This is an intentional design decision as there are many more edges
     in a graph than vertices, and as graph sizes grow large, it could become impractical
     to hold a single index table mapping each unique edge id to the data / location
     of the edge. Instead, we opt to switch to the more "natural" 
     method of identifying each edge simply by its source and destination 
     vertices.

     <b> Graph Creation </b>
  
     Vertices and edges can be added using the disk_graph::add_vertex()
     and graph::add_edge() member functions.
  
  
     \code
     vertex_it disk_graph::add_vertex(const VertexData& vdata = VertexData()) 
     void disk_graph::add_edge(vertex_id_t source, vertex_id_t target, const EdgeData& edata = EdgeData()) 
     \endcode
  
     There is also a "hinted" version of add_vertex which allows the user to 
     provide information about which atom should store the vertex.
  
     \code
     vertex_it disk_graph::add_vertex(const VertexData& vdata, size_t locationhint)
     \endcode
  
     However, these functions add one vertex or one edge at a time and
     could be slow. Vertex / edge creation can often be parallelized using the following 
     parallel creation functions.
  
     <b> Parallel Graph Creation </b>

     \code
     vertex_id_t add_vertex_collection(std::string collectionname,
     std::vector<VertexData> & vdata)

     vertex_id_t add_vertex_collection(std::string collectionname,
     size_t numvertices,
     boost::function<VertexData (vertex_id_t)> generator)
     \endcode

     The first form inserts a collection of vertices, assigning the collection a name
     and obtaining the vertex data from the vdata vector, 
     The function returns the index of the first vertex inserted. Remaining vertices
     are assigned the consecutive vertex IDs.
   
     The second form does the same but obtains the vertex data by calling a user provided
     function "VertexData generator(vertex_id_t)"
   
   
     Edges can then be inserted using the following function:
     \code
     void add_edge_indirect(std::string collectionname,
     boost::function<void (vertex_id_t,
     const VertexData&,
     std::vector<vertex_id_t>&, 
     std::vector<EdgeData>&,
     std::vector<vertex_id_t>&, 
     std::vector<EdgeData>&)> efn) 
     \endcode
     This rather complicated looking function prototype is actually quite straightforward.
   
     The provided user function will be called on each vertex in the collection,
     passing the user function the data on the vertex (vdata) as well as the
     ID of the vertex.

     The function should then return in the remaining arguments, the 
     set of outgoing vertices and their corresponding edge data, as well as the
     set of incoming vertices and their corresponding edge data.
   
     The prototype for the function is:
     \code
     void efn(vertex_id_t vid,          // input: vertex ID of the vertex being queried
     const VertexData& vdata,   // input: vertex data of the vertex being queried
     std::vector<vertex_id_t>& inv,    // output: Incoming vertices of the vertex being queried
     std::vector<EdgeData>& inedata,   // output: Data on the in-edges of the vertex. Paired with inv.
     std::vector<vertex_id_t>& outv,   // output: Outgoing vertices of the vertex being queried
     std::vector<EdgeData>& outedata); // output: Data on the out-edges of the vertex. Paired with outv.
     \endcode

     <b> Finalization </b>

     The finalize() function completes synchronization
     of disk IO as well as regenerate an atom index file.
     finalize() is automatically called in the disk_graph destructor.

     <b> Performance </b>

     Simple synthetic benchmarks have rated the edge insertion rate at about 100K to 200K edges
     per second. Vertex insertion rate can exceed 500K vertices per second,.
  */
  template<typename VertexData, typename EdgeData> 
  class disk_graph {
  public:
    //! \todo: Fix hack by defining proper types for this graph
    typedef graph<VertexData, EdgeData> base_graph_type;
    typedef typename base_graph_type::vertex_id_type    vertex_id_type;
    typedef typename base_graph_type::edge_id_type      edge_id_type;
    typedef typename base_graph_type::vertex_color_type vertex_color_type;
    typedef typename base_graph_type::edge_list_type    edge_list_type;
    /**
     * Creates or opens a disk graph with base name 'fbasename' using 
     * 'numfiles' of atoms. The atoms will be named fbasename.0, fbasename.1
     * fbasename.2, etc. An atom index file will be created in fbasename.idx.
     * This is useful when either creating a new
     * disk graph, or to open a disk graph without an atom index.
     * 
     * \note finalize() will always rebuild an atom index file irregardless of
     * which constructor is used.
     */
    disk_graph(std::string fbasename, size_t numfiles, 
               disk_graph_atom_type::atom_type atype = disk_graph_atom_type::DISK_ATOM) {  
      atoms.resize(numfiles);
      atomtype = atype;
      numv.value = 0;
      nume.value = 0;
      
      for (size_t i = 0;i < numfiles; ++i) {
        if (atomtype == disk_graph_atom_type::DISK_ATOM) {
          atoms[i] = new disk_atom(fbasename + "." + tostr(i), i);
          numv.value += atoms[i]->num_vertices();
          nume.value += atoms[i]->num_edges();
        }
        else if (atomtype == disk_graph_atom_type::MEMORY_ATOM) {
          atoms[i] = new memory_atom(fbasename + "." + tostr(i) + ".fast", i);
          numv.value += atoms[i]->num_vertices();
          nume.value += atoms[i]->num_edges();
        }
        else if (atomtype == disk_graph_atom_type::WRITE_ONLY_ATOM) {
          atoms[i] = new write_only_disk_atom(fbasename + "." + tostr(i) + ".dump", i, true);
        }
      }
      indexfile = fbasename + ".idx";
      ncolors = vertex_color_type(-1);
    }
      
    /**
     * Opens a graph using 'atomindex' as the disk graph index file.
     * If the file 'atomindex' does not exist, an assertion failure will be issued
     * Use the other constructor to create a new disk graph.
     */    
    disk_graph(disk_graph_atom_type::atom_type atype, std::string atomindex) {
      atomtype = atype;
      indexfile = atomindex;
    
      atom_index_file idxfile;
      idxfile.read_from_file(atomindex);

      ncolors = idxfile.ncolors;

      atoms.resize(idxfile.atoms.size());
      numv.value = 0;
      nume.value = 0;

      for (size_t i = 0;i < idxfile.atoms.size(); ++i) {
        if (atomtype == disk_graph_atom_type::DISK_ATOM) {
          atoms[i] = new disk_atom(idxfile.atoms[i].file, i);
          numv.value += atoms[i]->num_vertices();
          nume.value += atoms[i]->num_edges();
        }
        else if (atomtype == disk_graph_atom_type::MEMORY_ATOM) {
          atoms[i] = new memory_atom(idxfile.atoms[i].file + ".fast", i);
          numv.value += atoms[i]->num_vertices();
          nume.value += atoms[i]->num_edges();
        }
        else if (atomtype == disk_graph_atom_type::WRITE_ONLY_ATOM) {
          atoms[i] = new write_only_disk_atom(idxfile.atoms[i].file + ".dump", i, true);
        }
      }
      if (atomtype != disk_graph_atom_type::WRITE_ONLY_ATOM) {
        ASSERT_EQ(numv.value, idxfile.nverts);
        ASSERT_EQ(nume.value, idxfile.nedges);
      }
    }
  
    void create_from_graph(const graph<VertexData, EdgeData> &g,
                           const std::vector<vertex_id_type> &partids) {
      clear();
      size_t nv = g.num_vertices();
      logger(LOG_WARNING, "storing vertices...");
#pragma omp parallel for 
      for (int i = 0;i < (int)nv; ++i) {
        vertex_id_type vid = i;
        size_t hashloc = (vid) % atoms.size();
        //place vertices sequentially
        uint16_t owner = partids[i] % atoms.size();
        atoms[hashloc]->set_owner(vid, owner);
        atoms[owner]->add_vertex(vid, owner, g.vertex_data(i));
        atoms[owner]->set_color(vid, g.color(vid));
      }
      logger(LOG_WARNING, "storing edges...");
#pragma omp parallel for 
      for (int i = 0;i < (int)(g.num_edges()); ++i) {
        vertex_id_type target = g.target(i);
        vertex_id_type source = g.source(i);
        uint16_t sourceowner = partids[source] % atoms.size();
        uint16_t targetowner = partids[target] % atoms.size();
        if (sourceowner != targetowner) atoms[sourceowner]->add_edge(source, sourceowner, target, targetowner);
        atoms[targetowner]->add_edge(source, sourceowner, target, targetowner, g.edge_data(i));
      }
    
      numv.value = g.num_vertices();
      nume.value = g.num_edges();
      finalize();
    }
  
    disk_graph<VertexData, EdgeData>& operator=(const graph<VertexData, EdgeData> &g) {
      std::vector<vertex_id_type> partids(g.num_vertices());
      size_t nv = g.num_vertices();
      for (size_t i = 0;i < nv; ++i) {
        partids[i] = ((i * atoms.size())/ nv);
      }
      create_from_graph(g, partids);
      return *this;
    }
  
  
    ~disk_graph() {
      finalize();
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
  
    /**
       Synchronizes all disk atoms. This also recomputes the graph coloring
       if it is not already available. An atom index file is also regenerated.
    */
    void finalize() { 
      // synchronize all atoms
      for (size_t i = 0;i < atoms.size(); ++i) {
        atoms[i]->synchronize();
      }
      // compute ncolors if missing
      if (ncolors == vertex_color_type(-1)) {
        ncolors = 0;
        for (size_t i = 0;i < atoms.size(); ++i) {
          ncolors = std::max(ncolors, atoms[i]->max_color());
        }
        ++ncolors;
      }
      // generate an atom index file
      atom_index_file idx;
      idx.nverts = num_vertices();
      idx.nedges = num_edges();
      idx.natoms = atoms.size();
      idx.ncolors = ncolors;
      idx.atoms.resize(atoms.size());
      for (size_t i = 0;i < atoms.size(); ++i) {
        idx.atoms[i].protocol = "file";
        idx.atoms[i].file = atoms[i]->get_filename();
        // if end with .fast, strip it out
        if (idx.atoms[i].file.length() >= 5 && 
            (idx.atoms[i].file.substr(idx.atoms[i].file.length() - 5, 5) == ".fast" ||
            idx.atoms[i].file.substr(idx.atoms[i].file.length() - 5, 5) == ".dump")) {
          idx.atoms[i].file = idx.atoms[i].file.substr(0, idx.atoms[i].file.length() - 5);
        }
        idx.atoms[i].nverts = atoms[i]->num_vertices();
        idx.atoms[i].nedges = atoms[i]->num_edges();
        std::map<uint16_t, uint32_t> adj = atoms[i]->enumerate_adjacent_atoms();
        std::map<uint16_t, uint32_t>::iterator iter = adj.begin();
        while (iter != adj.end()) {
          idx.atoms[i].adjatoms.push_back(iter->first);
          idx.atoms[i].optional_weight_to_adjatoms.push_back(iter->second);
          ++iter;
        }
      }
      idx.write_to_file(indexfile);
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
    size_t num_in_neighbors(vertex_id_type v) const {
      uint16_t owner = atoms[v % atoms.size()]->get_owner(v);
      ASSERT_NE(owner, (uint16_t)(-1));
      return atoms[owner]->get_in_vertices(v).size();
    }
  
    /** \brief Get the number of out edges of a particular vertex */
    size_t num_out_neighbors(vertex_id_type v) const {
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
    vertex_id_type add_vertex(const VertexData& vdata = VertexData()) {
      vertex_id_type v = numv.inc_ret_last();
      uint16_t owner = v % atoms.size();
      atoms[owner]->add_vertex(v, owner, vdata);
      atoms[owner]->set_owner(v, owner);
      return v;
    }
  
    /** 
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     * Vertex is stored in atom locationhint
     */
    vertex_id_type add_vertex(const VertexData& vdata, 
                              uint16_t locationhint) {
      vertex_id_type v = numv.inc_ret_last();
      uint16_t owner = locationhint;
      atoms[owner]->add_vertex(v, owner, vdata);
      atoms[v % atoms.size()]->set_owner(v, owner);
      return v;
    }
  
    /**
     * Adds a vertex with vertex ID vid.
     * This function should be used with care as it will conflict
     * with the other add_vertex functions which automatically assign vertex_ids
     * Generally, this function is not safe to use together with the other
     * add_vertex functions.
     */
    void add_vertex_unsafe(vertex_id_type vid,
                           const VertexData& vdata, 
                           uint16_t locationhint) {
      uint16_t owner = locationhint;
      atoms[owner]->add_vertex(vid, owner, vdata);
      atoms[vid % atoms.size()]->set_owner(vid, owner);
    }
  
    /**
     * Adds a collection vdata.size() 
     * The collection of vertices will be assigned the name 'collectionname'
     * This name must be unique and cannot be zero length.
     * The return value is the ID of the first vertex inserted.
     * The remaining vertices have the consecutive vertex IDs.
     */
    vertex_id_type add_vertex_collection(std::string collectionname,
                                         std::vector<VertexData> & vdata) {
      vertex_id_type ret = numv.inc_ret_last(vdata.size());
      collection_range[collectionname] = std::make_pair(ret, ret + vdata.size());

#pragma omp parallel for
      for (int i = 0;i < (int)(vdata.size()); ++i) {
        vertex_id_type vid = i + ret;
        size_t hashloc = (vid) % atoms.size();
        //place vertices sequentially
        uint16_t owner = (uint16_t)((i * atoms.size())/ vdata.size());
        atoms[hashloc]->set_owner(vid, owner);
        atoms[owner]->add_vertex(vid, owner, vdata[i]);
      }
      return numv;
    }

    vertex_id_type add_vertex_collection(std::string collectionname,
                                         size_t numvertices,
                                         boost::function<VertexData (vertex_id_type)> generator) {
      vertex_id_type ret = numv.inc_ret_last(numvertices);
      collection_range[collectionname] = std::make_pair(ret, ret + numvertices);
#pragma omp parallel for
      for (int i = 0;i < (int)numvertices; ++i) {
        vertex_id_type vid = i + ret;
        size_t hashloc = (vid) % atoms.size();
        //place vertices sequentially
        uint16_t owner = (uint16_t)((i * atoms.size())/ numvertices);
        atoms[hashloc]->set_owner(vid, owner);
        atoms[owner]->add_vertex(vid, owner, generator((vertex_id_type)i));
      }
      return ret;
    }

  
  
    /**
     * \brief Creates an edge connecting vertex source to vertex target. 
     *        Any existing data will be cleared.
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData()) {
      nume.inc();
      uint16_t sourceowner = atoms[source % atoms.size()]->get_owner(source);
      ASSERT_NE(sourceowner, (uint16_t)(-1));
      uint16_t targetowner = atoms[target % atoms.size()]->get_owner(target);
      ASSERT_NE(targetowner, (uint16_t)(-1));
      if (sourceowner != targetowner) atoms[sourceowner]->add_edge(source, sourceowner, target, targetowner);
      atoms[targetowner]->add_edge(source, sourceowner, target, targetowner, edata);
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
     * void efn(vertex_id_type vid, \n 
     *         const VertexData& vdata, \n
     *         std::vector<vertex_id_type>& inv,  \n
     *         std::vector<EdgeData>& inedata, \n
     *         std::vector<vertex_id_type>& outv,  \n
     *         std::vector<EdgeData>& outedata) \n
     */
    void add_edge_indirect(std::string collectionname,
                           boost::function<void (vertex_id_type,
                                                 const VertexData&,
                                                 std::vector<vertex_id_type>&, 
                                                 std::vector<EdgeData>&,
                                                 std::vector<vertex_id_type>&, 
                                                 std::vector<EdgeData>&)> efn) {
      std::pair<vertex_id_type, vertex_id_type> vidrange = collection_range[collectionname];
    
#pragma omp parallel for 
      for (int vid = (int)vidrange.first; vid < (int)vidrange.second; ++vid) {
        uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
        ASSERT_NE(owner, (uint16_t)(-1));
        VertexData vdata;
        ASSERT_TRUE(atoms[owner]->get_vertex(vid, owner, vdata));
        std::vector<vertex_id_type> inv;
        std::vector<EdgeData> inedata;
        std::vector<vertex_id_type> outv; 
        std::vector<EdgeData> outedata;
        efn(vid, vdata, inv, inedata, outv, outedata);
        for (size_t i = 0;i < inv.size(); ++i) add_edge(inv[i], vid, inedata[i]);
        for (size_t i = 0;i < outv.size(); ++i) add_edge(vid, outv[i], outedata[i]);
      }
    }
    
    /**
     * Inserts an edge using explicitly specified locations for the source
     * and target. Faster, but potentially unsafe.
     */
    void add_edge_explicit(vertex_id_type source, uint16_t sourceowner,
                           vertex_id_type target, uint16_t targetowner,
                           const EdgeData& edata = EdgeData()) {
      nume.inc(); 
      if (sourceowner != targetowner) atoms[sourceowner]->add_edge(source, sourceowner, target, targetowner);
      atoms[targetowner]->add_edge(source, sourceowner, target, targetowner, edata);
    }
    
    
    void set_color(vertex_id_type vid, vertex_color_type color) {
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      atoms[owner]->set_color(vid, color);
      // reset ncolors. we will need to recompute it on save
      ncolors = vertex_color_type(-1); 
    }

    void set_color_unsafe(vertex_id_type vid, vertex_color_type color, uint16_t locationhint) {
      uint16_t owner = locationhint;
      ASSERT_NE(owner, (uint16_t)(-1));
      atoms[owner]->set_color(vid, color);
      // reset ncolors. we will need to recompute it on save
      ncolors = vertex_color_type(-1); 
    }

    
    /** \brief Returns the vertex color of a vertex.
        Coloring is only valid if compute_coloring() is called first.*/
    vertex_color_type get_color(vertex_id_type vid) {
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      uint32_t vc = atoms[owner]->get_color(vid);
      // color is not set in the file! return 0
      if (vc == (uint32_t)(-1)) {
        return 0; 
      }
      else return vc;
      
    }
    
    
    /** \brief This function constructs a heuristic coloring for the 
        graph and returns the number of colors */
    size_t compute_coloring() {
      // Reset the colors
#pragma omp parallel for
      for(int v = 0; v < (int)num_vertices(); ++v) set_color(v, 0);
      
      // Recolor
      size_t max_color = 0;
      std::set<vertex_color_type> neighbor_colors;
      for(vertex_id_type vid = 0; vid < num_vertices(); ++vid) {
        neighbor_colors.clear();
        // Get the neighbor colors
      foreach(vertex_id_type neighbor_vid, in_vertices(vid)){
        const vertex_color_type& neighbor_color = get_color(neighbor_vid);
        neighbor_colors.insert(neighbor_color);
      }
      foreach(vertex_id_type neighbor_vid, out_vertices(vid)){
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
      ncolors = max_color + 1;
    return max_color + 1;     
    }
  
    std::vector<vertex_id_type> in_vertices(vertex_id_type vid) const {
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      return atoms[owner]->get_in_vertices(vid);
    }
    
    std::vector<vertex_id_type> out_vertices(vertex_id_type vid) const {
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      return atoms[owner]->get_out_vertices(vid);
    }
  
    /** \brief Gets the edge data of the edge from vertex 'src' to vertex 'dest'*/
    EdgeData get_edge_data(vertex_id_type src, vertex_id_type dest) const {
      EdgeData ret;
      uint16_t owner = atoms[dest % atoms.size()]->get_owner(dest);
      ASSERT_NE(owner, (uint16_t)(-1));
      ASSERT_TRUE(atoms[owner]->get_edge(src, dest, ret));
      return ret;
    }
    
    /** \brief Sets the edge data of the edge from vertex 'src' to vertex 'dest'*/
    void set_edge_data(vertex_id_type src, vertex_id_type dest, const EdgeData &edata) {
      uint16_t owner = atoms[dest % atoms.size()]->get_owner(dest);
      ASSERT_NE(owner, (uint16_t)(-1));
      atoms[owner]->set_edge(src, dest, edata);
    }
    
    /** \brief Gets the vertex data of the vertex with id 'vid' */
    VertexData get_vertex_data(vertex_id_type vid) {
      VertexData ret;
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      ASSERT_TRUE(atoms[owner]->get_vertex(vid, owner, ret));
      return ret;
    }

    /** \brief Sets the vertex data of the vertex with id 'vid' */  
    void set_vertex_data(vertex_id_type vid, const VertexData& vdata) {
      uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
      ASSERT_NE(owner, (uint16_t)(-1));
      atoms[owner]->set_vertex(vid, owner, vdata);
    }
    
    size_t num_atoms() const {
      return atoms.size();
    }
    
    void make_memory_atoms() {
//      #pragma omp parallel for
      for (int i = 0;i < (int)atoms.size(); ++i) {
        std::string fname = atoms[i]->get_filename();
        // Make sure that this is not already a fast file
        if (fname.length() < 5 || fname.substr(fname.length() - 5, 5) != ".fast") {
          if (typeid(*atoms[i]) == typeid(disk_atom)) {
            dynamic_cast<disk_atom*>(atoms[i])->build_memory_atom(atoms[i]->get_filename() + ".fast");
          }
          else {
            std::string mfile = fname.substr(0, fname.length() - 5) + ".fast";
            memory_atom matom(mfile, atoms[i]->atom_id());
            matom.clear();
            dynamic_cast<write_only_disk_atom*>(atoms[i])->play_back(&matom);
          }
        }
      }
    }
    
  private:
    
    // block copies
  disk_graph(const disk_graph&);
    disk_graph& operator=(const disk_graph&);
    
    mutable std::vector<graph_atom*> atoms;
    atomic<size_t> numv;
    atomic<size_t> nume;
    
    std::map<std::string, std::pair<vertex_id_type, vertex_id_type> > collection_range;
    
    std::string indexfile;
    vertex_color_type ncolors; // this is (-1) if it is not set
    
    disk_graph_atom_type::atom_type atomtype;
    
    void propagate_coloring() {
#pragma omp parallel for
    for (int i = 0;i < (int)atoms.size(); ++i) {
      foreach(vertex_id_type vid, atoms[i]->enumerate_vertices()) {
        uint16_t owner = atoms[vid % atoms.size()]->get_owner(vid);
        if (owner != i) {
          atoms[i]->set_color(vid, atoms[owner]->get_color(vid));
        }
      }
    }
  }
    
  }; // end of Disk Graph
  























  /**
   * Quick utility function to quickly turn a graph into a disk graph
   * parts is a vector assigning each vertex to an atom ID. atom IDs begin 
   * at 0. The total number of atoms created is 1 + the largest 
   * value in the parts vector.
   *
   * The result graph will be saved using the basename provided. The
   * atom index will have file name "[basename].idx" and the atoms 
   * will be stored in "[basename].0", "[basename].1", etc.
   */
  template<typename VertexData, typename EdgeData>
  void graph_partition_to_atomindex(graph<VertexData, EdgeData> &g,
                                    std::vector<graph_partitioner::part_id_type>& parts,
                                    std::string basename) {
    // get the number of atoms
    size_t nparts = 0;
    for (size_t i = 0;i < parts.size(); ++i) {
      nparts = std::max<size_t>(nparts, parts[i]);
    }
    ++nparts;
    disk_graph<VertexData, EdgeData> dg(basename, nparts + 1, 
                              graphlab::disk_graph_atom_type::WRITE_ONLY_ATOM);
    dg.create_from_graph(g, parts);
    dg.finalize();
  }                           
  
  
} // namespace graphlab

#include <graphlab/macros_undef.hpp>
#endif

