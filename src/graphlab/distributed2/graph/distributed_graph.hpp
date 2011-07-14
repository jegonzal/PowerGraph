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


#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_HPP

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/caching_dht.hpp>
#include <graphlab/rpc/lazy_dht.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/metrics/metrics.hpp>
#include <graphlab/graph/atom_index_file.hpp>
#include <graphlab/graph/disk_atom.hpp>
#include <graphlab/distributed2/graph/dgraph_edge_list.hpp>
#include <graphlab/distributed2/graph/graph_local_store.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {


  // forward definition of a friend
  template <typename GraphType> class graph_lock;


  /**
   * \brief Distributed Graph Implementation.
   * 
   * A fully distributed implementation of the \ref graph object.
   * Vertices are partitioned across machines. Each vertex is owned by a
   * unique machine. Each edge is owned by its destination vertex.  Each
   * machine stores all vertex data for vertices within its partition,
   * as well as vertex data for vertices/edges on the boundary of the
   * partition.  Each vertex data instance is therefore replicated as
   * many times as the number of distinct machines owning neighbors of
   * the vertex in question.
   * 
   * The distributed graph structure is <b> not mutable </b>. Only the graph
   * data is mutable.
   * 
   * Formally, where \f$ \Gamma(v)\f$ is the set of neighbors of \f$
   * v\f$ and \f$o(v)\f$ is the owner of vertex v, vertex v is
   * replicated \f$ \mbox{DISTINCT} \left( \left\{ o(u), u \in \left\{
   * v,\Gamma(v) \right\} \right\} \right) \f$ times.
   * 
   * Each edge is replicated a maximum of 2 times.
   * 
   * To standardize on terminology, we call the set of vertices and
   * edges owned by a machine, the machine's <b> partition </b>. We call
   * the set of vertices and edges adjacent to the partition (but not in
   * the partition), the <b> boundary </b>. Finally, we will call a
   * machine's local copy of the partition + boundary, the machine's <b>
   * fragment </b>.
   * 
   * 
   *
   * Vertex IDs: 
   * 
   * Every vertex in the graph has a uniquely assigned global
   * vertex. 
   * 
   * Each machine has a local representation for its fragment of the
   * graph. Within the local fragment, each vertex/edge has a local
   * vertex/edge ID. This local ID is hidden and abstracted from the
   * user. Implementors however should keep in mind the following
   * requirements for the local representation:
   *
   * <ul>
   * <li> Local vertex / edge IDs are unique and sequentially assigned </li>
   * <li> Sorting all vertices/edges in the local fragment must
   *      produce the same sequence whether or not we sort by global IDs 
   *      or Local IDs. </li>
   * </ul>
   * 
   * 
   * Graph Data:
   * 
   * Accessing and modifying graph data through get_edge_data_from_pair(), 
   * get_vertex_data(), set_edge_data_from_pair() and set_vertex_data() functions
   * are fully parallel and can be called simultaneously from any machine.
   * Remote accesses are automatically managed.
   * 
   * The regular vertex_data() and edge_data() functions which return references
   * are still available. However, as a matter of necessity, only local vertices /
   * local edges are accessible through these functions. An assertion failure
   * will be raised if remote vertices / remote edges are accessed through 
   * vertex_data() or edge_data().
   * 
   * Consistency: 
   * 
   * Consistency of graph data, is not managed and must be done manually
   * through the various synchronize() operations. All data reads will
   * be accessed through the local fragment if the local fragment
   * contains the data. Otherwise, it will be requested from the owner
   * of the data. All data writes will be sent to the owner of the
   * data. The writes may not however, update all fragments unless
   * explicitly requested.
   * 
   * 
   */
  template<typename VertexData, typename EdgeData> 
  class distributed_graph {
 
  public:
    typedef VertexData vertex_data_type;
    typedef EdgeData edge_data_type;


    typedef dist_graph_impl::graph_local_store<VertexData, EdgeData>
    graph_local_store_type;
    
    typedef typename graph_local_store_type::vertex_id_type 
    vertex_id_type;
    typedef typename graph_local_store_type::vertex_color_type
    vertex_color_type;
    
    typedef typename graph_local_store_type::edge_id_type   
    edge_id_type;
    typedef typename graph_local_store_type::edge_list_type 
    local_edge_list_type;

  


    /** This class defines a set of edges */
    class edge_list {
    public:

      typedef typename boost::function<edge_id_type (edge_id_type)> 
      TransformType;
      typedef typename boost::transform_iterator<TransformType, const edge_id_type*> 
      iterator;
      typedef typename boost::transform_iterator<TransformType, const edge_id_type*> 
      const_iterator;
      typedef edge_id_type value_type;
    private:
      
      template<typename T>
      inline static T identity(T x) { return x; }
      
      // this class can wrap two methods of describing an edge list
      // method 1: a regular edge_list which requires some conversion
      // method 2: a vector
      
      // method 1:
      local_edge_list_type elist;
      // method 2:
      std::vector<edge_id_type> edgeidvec;
      size_t numel;
      bool useelist;
      
    public:    
      // an edge list which wraps a regular elist. Method 1.
      inline edge_list(local_edge_list_type elist) : 
        elist(elist), numel(elist.size()), useelist(true){ }
      
      // an edge list which wraps a vector. Method 2.
      inline edge_list(const std::vector<edge_id_type> &edgeidvec) :
        edgeidvec(edgeidvec), numel(edgeidvec.size()), useelist(false) { }

      /** \brief Get the size of the edge list */
      inline size_t size() const { return numel; }

      /** \brief Get the ith edge in the edge list */
      inline edge_id_type operator[](size_t i) const {
        assert(i < size());
        return *(begin() + i);
      }

      /** \brief Returns a pointer to the start of the edge list */
      inline iterator begin() const {
        if (useelist) {
          return boost::
            make_transform_iterator(elist.begin(),
                                    identity<edge_id_type>);
        } else {
          return boost::
            make_transform_iterator(&(edgeidvec[0]),
                                    identity<edge_id_type>);
        }
      }

      /** \brief Returns a pointer to the end of the edge list */
      inline iterator end() const {
        if (useelist) {
          return boost::
            make_transform_iterator(elist.end(),
                                    identity<edge_id_type>);
        } else {
          return boost::
            make_transform_iterator(&(edgeidvec[edgeidvec.size()]),
                                                identity<edge_id_type>);

        }
      }

      /** \brief Fill a vector with edge id list */
      inline void fill_vector(std::vector<edge_id_type>& lvalue) const {
        lvalue.clear();
        foreach(edge_id_type eid, *this) lvalue.push_back(eid);
      }

      /** \brief test if the edge list is empty */
      inline bool empty() const { return size() == 0; }
    }; // end edge_list


    typedef edge_list edge_list_type;




    /**
     * Constructs a distributed graph loading the graph from the atom index
     * 'indexfilename'
     */
    distributed_graph(distributed_control &dc, 
                      std::string indexfilename, 
                      bool do_not_load_data = false,
                      bool do_not_mmap = true,
                      bool sliced_partitioning = false):
      rmi(dc, this),
      indexfilename(indexfilename),
      globalvid2owner(dc, 65536),
      pending_async_updates(true, 0),
      pending_push_updates(true, 0),
      graph_metrics("distributed_graph"){

      if (do_not_mmap == false) {
        logstream(LOG_WARNING) << "Using MMAP for local storage is highly experimental" << std::endl;
      }
                                
      cur_proc_vector.push_back(rmi.procid());
      // read the atom index.
      atom_index_file atomindex;
      atomindex.read_from_file(indexfilename);
      // store the graph size
      numglobalverts = atomindex.nverts;
      numglobaledges = atomindex.nedges;
      numcolors = atomindex.ncolors;
      // machine 0 partitions it
      if (dc.procid() == 0) {
        if (sliced_partitioning == false) {
          atompartitions = partition_atoms(atomindex, dc.numprocs());
        }
        else {
          atompartitions = partition_atoms_sliced(atomindex, dc.numprocs());
        }

        for (size_t i = 0;i < atompartitions.size(); ++i) {
          logstream(LOG_DEBUG) << i << ":";
          for (size_t j = 0; j < atompartitions[i].size(); ++j) {
            logstream(LOG_DEBUG) << atompartitions[i][j] << ", ";
          }
          logstream(LOG_DEBUG) << std::endl;
        }
      }
      dc.barrier();
      rmi.broadcast(atompartitions, dc.procid() == 0);
      construct_local_fragment(atomindex, atompartitions, rmi.procid(), do_not_load_data, do_not_mmap);
      rmi.barrier();
    }

    ~distributed_graph() {
      rmi.barrier();
    }

  
  
    /// Saves the local fragment back to disk
    void save() {
      atom_index_file atomindex;
      atomindex.read_from_file(indexfilename);
    
      std::vector<size_t>& atoms_in_curpart = atompartitions[rmi.procid()];

#pragma omp parallel for
      for (int i = 0;i < (int)atoms_in_curpart.size(); ++i) {
        size_t atomid = atoms_in_curpart[i];
        disk_atom atom(atomindex.atoms[atomid].file, atomid);
        // go through the set of owned vertices
        foreach(vertex_id_type vid, atom.enumerate_vertices()) {
          // update the data on the vertex. I should own this vertex
          atom.set_vertex(vid, atomid, vertex_data(vid));
          // go through the set of owned edges
          foreach(vertex_id_type src, atom.get_in_vertices(vid)) {
            // update the data on the edge. I should own this edge
            atom.set_edge(src, vid, edge_data(src, vid));
          }
        }
      }
    }
  

    /// Returns the number of vertices in the graph.
    size_t num_vertices() const{
      return numglobalverts;
    }
  
    size_t local_vertices() const {
      return owned_vertices().size();
    }

    /// Returns the number of edges in the graph.
    size_t num_edges() const{
      return numglobaledges;
    }

    /// Returns the number of incoming edges of vertex 'vid'
    size_t num_in_neighbors(vertex_id_type vid) const {
      typename global2localvid_type::const_iterator iter = 
        global2localvid.find(vid);
      // if I have this vertex in my fragment
      if (iter != global2localvid.end()) {
        // and if I own it (it is interior)
        vertex_id_type localvid = iter->second;
        if (localvid2owner[localvid] == rmi.procid()) {
          return localstore.num_in_neighbors(localvid);
        }
      }
      std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
      assert(vidowner.first);
      // otherwise I need to ask the owner
      return rmi.remote_request(vidowner.second,
                                &distributed_graph<VertexData, EdgeData>::
                                num_in_neighbors,
                                vid);
    }

    /// Returns the number of outgoing edges of vertex 'vid'
    size_t num_out_neighbors(vertex_id_type vid) const {
      typename global2localvid_type::const_iterator iter = 
        global2localvid.find(vid);
      // if I have this vertex in my fragment
      if (iter != global2localvid.end()) {
        // and if I own it (it is interior)
        vertex_id_type localvid = iter->second;
        if (localvid2owner[localvid] == rmi.procid()) {
          return localstore.num_out_neighbors(localvid);
        }
      }

      std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
      assert(vidowner.first);

      // otherwise I need to ask the owner
      return rmi.remote_request(vidowner.second,
                                &distributed_graph<VertexData, EdgeData>::
                                num_out_neighbors,
                                vid);
    }


    /** \brief Finds an edge return the edge id
        The value of the first element of the pair will be true if an 
        edge from src to target is found and false otherwise. If the 
        edge is found, the edge ID is returned in the second element of the pair. 
        The target vertex must be owned by this machine. */
    std::pair<bool, edge_id_type>
    find(vertex_id_type source, vertex_id_type target) const {
      std::pair<bool, edge_id_type> ret;
      // hmm. surprisingly tricky
      typename global2localvid_type::const_iterator itersource = 
        global2localvid.find(source);
      typename global2localvid_type::const_iterator itertarget = 
        global2localvid.find(target);
      // both are local, I can find it
      if (itersource != global2localvid.end() && 
          itertarget != global2localvid.end()) {
        ret = localstore.find(itersource->second, itertarget->second);
        // convert to global edge ids
        if (ret.first) ret.second = ret.second;
        return ret;
      }
      // if the edge exists, the owner of either the source or target must have it
      // lets use the target
      std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(target);
      assert(vidowner.first);
      procid_t targetowner = vidowner.second;
    
      // if I am the owner, then this edge can't possibly exist
      if (targetowner == rmi.procid()) {
        ret.first = false; ret.second = 0;
        return ret;
      }
      else {
        ASSERT_MSG(false, "edge is not in this fragment. eid not available");
      }
    }

    /// unsafe version of find
    edge_id_type edge_id(vertex_id_type source, vertex_id_type target) const {
      std::pair<bool, edge_id_type> res = find(source, target);
      // The edge must exist
      assert(res.first);
      return res.second;
    }

    /// Reverses an existing edge id
    edge_id_type rev_edge_id(edge_id_type eid) const {
      // do I have this edge in the fragment?
      return localstore.rev_edge_id(eid);
    } // end of rev_edge_id


    /** \brief Returns the source vertex of an edge. */
    vertex_id_type source(edge_id_type eid) const {
      // do I have this edge in the fragment?
      return local2globalvid[localstore.source(eid)];
    }

    /** \brief Returns the destination vertex of an edge. */
    vertex_id_type target(edge_id_type eid) const {
      return local2globalvid[localstore.target(eid)];
    }

    /// Returns the set of incoming vertices of vertex v
    std::vector<vertex_id_type> in_vertices(vertex_id_type v) const {
      if (is_owned(v)) {
        // switch to localeid
        vertex_id_type localvid = globalvid_to_localvid(v);
        //get the local edgelist
        local_edge_list_type elist = localstore.in_edge_ids(localvid);
     
        std::vector<vertex_id_type> ret(elist.size());
        // get the source of each edge and return it. remember to convert back to global
        for (size_t i = 0;i < elist.size(); ++i) {
          ret[i] = localvid_to_globalvid(localstore.source(elist[i]));
        }
        return ret;
      }
      else {
        return rmi.remote_request(globalvid2owner.get_cached(v).second,
                                  &distributed_graph<VertexData,EdgeData>::in_vertices,
                                  v);
      }
    }
  
    /// Returns the set of outgoing vertices of vertex v
    std::vector<vertex_id_type> out_vertices(vertex_id_type v) const {
      if (is_owned(v)) {
        // switch to localeid
        vertex_id_type localvid = globalvid_to_localvid(v);
        //get the local edgelist
        local_edge_list_type elist = localstore.out_edge_ids(localvid);
     
        std::vector<vertex_id_type> ret(elist.size());
        // get the source of each edge and return it. remember to convert back to global
        for (size_t i = 0;i < elist.size(); ++i) {
          ret[i] = localvid_to_globalvid(localstore.target(elist[i]));
        }
        return ret;
      }
      else {
        return rmi.remote_request(globalvid2owner.get_cached(v).second,
                                  &distributed_graph<VertexData,EdgeData>::out_vertices,
                                  v);
      }
    }

    /** \brief Return the edge ids of the edges arriving at v */
    edge_list_type in_edge_ids(vertex_id_type v) const {
      typename global2localvid_type::const_iterator iter = 
        global2localvid.find(v);
      // if I have the vertex in my fragment
      // and if it is interior
      if (iter != global2localvid.end()) {
        vertex_id_type localvid = iter->second;
        if (localvid2owner[localvid]  == rmi.procid()) {
          return edge_list_type(localstore.in_edge_ids(localvid));
        }
      }
      // ok I need to construct a vector
      return edge_list_type(in_edge_id_as_vec(v));
    } // end of in edges

    /** \brief Return the edge ids of the edges arriving at v as a vector */
    std::vector<edge_id_type> in_edge_id_as_vec(vertex_id_type v) const {
      typename global2localvid_type::const_iterator iter = 
        global2localvid.find(v);
      // if I have the vertex in my fragment
      // and if it is interior
      std::vector<edge_id_type> ret;
      if (iter != global2localvid.end()) {
        vertex_id_type localvid = iter->second;
        if (localvid2owner[localvid]  == rmi.procid()) {
          foreach(edge_id_type eid, localstore.in_edge_ids(localvid)) {
            ret.push_back(eid);
          }
          return ret;
        }
      }
      ASSERT_MSG(false, "Edge IDs of non-local vertex cannot be obtained");
    } // end of in edges

    /** \brief Return the edge ids of the edges leaving at v */
    edge_list_type out_edge_ids(vertex_id_type v) const {
      typename global2localvid_type::const_iterator iter = 
        global2localvid.find(v);
      // if I have the vertex in my fragment
      // and if it is interior
      if (iter != global2localvid.end()) {
        vertex_id_type localvid = iter->second;
        if (localvid2owner[localvid]  == rmi.procid()) {
          return edge_list_type(localstore.out_edge_ids(localvid));
        }
      }
      // ok I need to construct a vector
      return edge_list_type(out_edge_id_as_vec(v));
    } // end of out edges



    /** \brief Return the edge ids of the edges leaving at v as a vector*/
    std::vector<edge_id_type> out_edge_id_as_vec(vertex_id_type v) const {
      typename global2localvid_type::const_iterator iter = 
        global2localvid.find(v);
      // if I have the vertex in my fragment
      // and if it is interior
      std::vector<edge_id_type> ret;
      if (iter != global2localvid.end()) {
        vertex_id_type localvid = iter->second;
        if (localvid2owner[localvid]  == rmi.procid()) {
          foreach(edge_id_type eid, localstore.out_edge_ids(localvid)) {
            ret.push_back(eid);
          }
          return ret;
        }
      }
    
      std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(v);
      assert(vidowner.first);
    
      return rmi.remote_request(vidowner.second,
                                &distributed_graph<VertexData, EdgeData>::
                                out_edge_id_as_vec,
                                v);
    } // end of in edges


    /// Prints the local store
    void print(std::ostream &out) const {
      for (size_t i = 0;i < localstore.num_edges(); ++i) {
        std::cout << local2globalvid[localstore.source(i)] << ", " 
                  << local2globalvid[localstore.target(i)] << "\n";
      }
    }

    vertex_id_type globalvid_to_localvid(vertex_id_type vid) const {
      typename global2localvid_type::
        const_iterator iter = global2localvid.find(vid);
      assert(iter != global2localvid.end());
      return iter->second;
    }

    vertex_id_type localvid_to_globalvid(vertex_id_type vid) const {
      return local2globalvid[vid];
    }
  
    procid_t globalvid_to_owner(vertex_id_type vid) const {
      if (vertex_is_local(vid)) {
        typename global2localvid_type::
          const_iterator iter = global2localvid.find(vid);
        if (iter == global2localvid.end()) return false;
        vertex_id_type localvid = iter->second;
        return localvid2owner[localvid];
      }
      else {
        std::pair<bool, procid_t> ret = globalvid2owner.get_cached(vid);
        assert(ret.first);
        return ret.second;
      }
    }

    bool vertex_is_local(vertex_id_type vid) const{
      return global_vid_in_local_fragment(vid);
    }
  
  
    bool is_ghost(vertex_id_type vid) const{
      typename global2localvid_type::
        const_iterator iter = global2localvid.find(vid);
      if (iter == global2localvid.end()) return false;
      return localvid2owner[iter->second] != rmi.procid();
    }
  
  
    bool localvid_is_ghost(vertex_id_type localvid) const {
      return localvid2owner[localvid] != rmi.procid();
    }

    bool is_owned(vertex_id_type vid) const{
      return vertex_is_local(vid) && !is_ghost(vid);
    }
  
    const std::vector<vertex_id_type>& owned_vertices() const{
      return ownedvertices;
    }

    const std::vector<vertex_id_type>& boundary_scopes() const{
      return boundaryscopes;
    }

    const boost::unordered_set<vertex_id_type>& boundary_scopes_set() const{
      return boundaryscopesset;
    }
  
    bool on_boundary(vertex_id_type vid) const{
      return boundaryscopesset.find(vid) != boundaryscopesset.end();
    }
 
    const std::vector<vertex_id_type>& ghost_vertices() const{
      return ghostvertices;
    }
  
    /** returns a vector of all processors having a replica of this globalvid
     *  This vector is guaranteed to be in sorted order of processors.
     */
    const std::vector<procid_t>& localvid_to_replicas(vertex_id_type localvid) const {
      if (localvid2ghostedprocs[localvid].empty()) {
        return cur_proc_vector;
      }
      else {
        return localvid2ghostedprocs[localvid];
      }
    }
  
    /// returns a vector of all processors having a replica of this globalvid
    const std::vector<procid_t>& globalvid_to_replicas(vertex_id_type globalvid) const {
      vertex_id_type localvid = globalvid_to_localvid(globalvid);
      return localvid_to_replicas(localvid);
    }
    /**
     * Returns a reference to the edge data on the edge source->target
     * assertion failure if the edge is not within the current fragment
     */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target) {
      ASSERT_TRUE(global_vid_in_local_fragment(source) && global_vid_in_local_fragment(target));
      return localstore.edge_data(global2localvid.find(source)->second,
                                  global2localvid.find(target)->second);
    }

    /**
     * Returns a constant reference to the edge data on the edge source->target
     * assertion failure if the edge is not within the current fragment
     */
    const EdgeData& edge_data(vertex_id_type source, vertex_id_type target) const{
      ASSERT_TRUE(global_vid_in_local_fragment(source) && global_vid_in_local_fragment(target));

      return localstore.edge_data(global2localvid.find(source)->second,
                                  global2localvid.find(target)->second);
    }

    /**
     * Returns a reference to the edge data on the edge eid
     * assertion failure if the edge is not within the current fragment
     */
    EdgeData& edge_data(edge_id_type eid) {
      return localstore.edge_data(eid);
    }

    /**
     * Returns a constant reference to the edge data on the edge eid
     * assertion failure if the edge is not within the current fragment
     */
    const EdgeData& edge_data(edge_id_type eid) const{
      return localstore.edge_data(eid);
    }

    /**
     * Returns a reference to the vertex data on vertex vid
     * assertion failure if the vertex is not within the current fragment
     */
    VertexData& vertex_data(vertex_id_type vid) {
      ASSERT_TRUE(global_vid_in_local_fragment(vid));
      return localstore.vertex_data(global2localvid[vid]);
    }

    void __attribute__((__deprecated__)) increment_vertex_version(vertex_id_type vid) {
   
    }

    void __attribute__((__deprecated__)) increment_edge_version(vertex_id_type eid) {

    }

    void vertex_is_modified(vertex_id_type vid) {
      ASSERT_TRUE(global_vid_in_local_fragment(vid));
      // If this is a ghost, this just marks as modified
      vertex_id_type localvid = global2localvid[vid];
      if (localvid2owner[localvid] == rmi.procid()) {
        localstore.increment_vertex_version(localvid);
      }
      else {    
        localstore.set_vertex_modified(localvid, true);
      }
    }

    void edge_is_modified(edge_id_type eid) {
      // If the edge is a ghost, this just marks as modified
      vertex_id_type localtargetvid = localstore.target(eid);
      if (localvid2owner[localtargetvid] == rmi.procid()) {
        localstore.increment_edge_version(eid);
      }
      else {    
        localstore.set_edge_modified(eid, true);
      }
    }

    void __attribute__((__deprecated__)) vertex_clear_modified(vertex_id_type vid) {
      //    localstore.set_vertex_modified(global2localvid[vid], false);
    }

    void __attribute__((__deprecated__)) edge_clear_modified(edge_id_type eid) {
      //   localstore.set_edge_modified(eid, false);  
    }


    bool __attribute__((__deprecated__)) is_vertex_modified(vertex_id_type vid) {
      return false;
    }

    bool __attribute__((__deprecated__)) is_edge_modified(edge_id_type eid) {
      return false;
    }

    /**
     * Returns a constant reference to the vertex data on vertex vid
     * assertion failure if the vertex is not within the current fragment
     */
    const VertexData& vertex_data(vertex_id_type vid) const{
      assert(global_vid_in_local_fragment(vid));
      return localstore.vertex_data(global2localvid.find(vid)->second);
    }


    /**
     * Returns a copy of the edge data on the edge source->target
     * If the edge is not on this fragment, the request is sent
     * to a remote machine.
     */
    EdgeData get_edge_data_from_pair(vertex_id_type source, 
                                     vertex_id_type target) const {
      if (global_vid_in_local_fragment(source) && 
          is_owned(target)) {
        return edge_data(source, target);
      }
      else {
        std::pair<bool, procid_t> vidowner = 
          globalvid2owner.get_cached(target);
        assert(vidowner.first);

        return rmi.remote_request(vidowner.second,
                                  &distributed_graph<VertexData,EdgeData>::
                                  get_edge_data_from_pair,
                                  source,
                                  target);
      }
    }

    /**
     * Returns a copy of the edge data on the edge eid
     * If the edge is not on this fragment, the request is sent
     * to a remote machine.
     */
    EdgeData get_edge_data_from_eid(edge_id_type eid) const{
      return edge_data(eid);
    }

    /**
     * Returns a copy of the edge data on the edge source->target
     * If the edge is not on this fragment, the request is sent
     * to a remote machine.
     */
    EdgeData get_edge_data(vertex_id_type source, vertex_id_type target) const {
      return get_edge_data_from_pair(source, target);
    }

    /**
     * Returns a copy of the edge data on the edge eid
     * If the edge is not on this fragment, the request is sent
     * to a remote machine.
     */
    EdgeData get_edge_data(edge_id_type eid) const{
      return get_edge_data_from_eid(eid);
    }

    /**
     * Returns a copy of the vertex data on the vertex vid
     * If the vertex is not on this fragment, the request is sent
     * to a remote machine.
     */
    VertexData get_vertex_data(vertex_id_type vid) const{
      if (global_vid_in_local_fragment(vid)) {
        return vertex_data(vid);
      }
      else {
        std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
        assert(vidowner.first);

        return rmi.remote_request(vidowner.second,
                                  &distributed_graph<VertexData,EdgeData>::
                                  get_vertex_data,
                                  vid);
      }
    }



    /**
     * Sets the data on the edge source->target If the edge is not on
     * this fragment, the request is sent to a remote machine. If async
     * is true, the function returns immediately without waiting for
     * confirmation from the remote machine.
     */
    void set_edge_data_from_pair(vertex_id_type source, vertex_id_type target,
                                 const EdgeData edata, bool async) {
      // sets must go straight to the owner
      typename global2localvid_type::const_iterator targetiter = 
        global2localvid.find(target);
      // if I own the target vertex, then I own the edge
      if (targetiter != global2localvid.end()) {
        if (localvid2owner[targetiter->second] == rmi.procid()) {
          edge_data(source, target) = edata;
          vertex_id_type localsourcevid = global2localvid[source];
          vertex_id_type localtargetvid = global2localvid[target];
          // get the localeid
          // edge must exist. I don't need to check the return of find
          edge_id_type eid = localstore.find(localsourcevid, localtargetvid).second;
          edge_is_modified(eid);
          push_owned_edge_to_replicas(eid,
                                      async,  // async  
                                      async); // untracked

          return;
        }
      }
      std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(target);
      assert(vidowner.first);

      if (async) {
        rmi.remote_call(vidowner.second,
                        &distributed_graph<VertexData,EdgeData>::
                        set_edge_data_from_pair,
                        source,
                        target,
                        edata,
                        async);
      }
      else {
        rmi.remote_request(vidowner.second,
                           &distributed_graph<VertexData,EdgeData>::
                           set_edge_data_from_pair,
                           source,
                           target,
                           edata,
                           async);
      }
    }

    /**
     * Sets the data on the edge eid
     * If the edge is not on this fragment, the request is sent
     * to a remote machine. If async is true, the function returns immediately
     * without waiting for confirmation from the remote machine.
     */
    void set_edge_data_from_eid(edge_id_type eid, 
                                const EdgeData edata, bool async){
      if (localvid2owner[localstore.target(eid)] == rmi.procid()) {
        // if I do. then I must own the edge.
        edge_data(eid) = edata;
        edge_is_modified(eid);
        push_owned_edge_to_replicas(eid, 
                                    async,  // async  
                                    async); // untracked

        return;
      }
  
      ASSERT_MSG(false,
                 "Remote edge set impossible due to use of canonical edge numbering");

    }

    /**
     * Sets the data on the edge source->target
     * If the edge is not on this fragment, the request is sent
     * to a remote machine. This operation is performed synchronously.
     * It will wait for the remote machine to complete the modification before
     * returning control.
     */
    void set_edge_data(vertex_id_type source, vertex_id_type target, 
                       const EdgeData edata) {
      set_edge_data_from_pair(source, target, edata, false);
    }

    /**
     * Sets the data on the edge eid
     * If the edge is not on this fragment, the request is sent
     * to a remote machine. This operation is performed synchronously.
     * It will wait for the remote machine to complete the modification before
     * returning control.
     */
    void set_edge_data(edge_id_type eid, const EdgeData edata){
      set_edge_data_from_eid(eid, edata, false);
    }

    /**
     * Sets the data on the vertex vid
     * If the vertex is not on this fragment, the request is sent
     * to a remote machine. This operation is performed synchronously.
     * It will wait for the remote machine to complete the modification before
     * returning control.
     */
    void set_vertex_data(vertex_id_type vid, const VertexData vdata){
      if (is_owned(vid)) {
        vertex_data(vid) = vdata;
        // increment version
        vertex_is_modified(vid);
        push_owned_vertex_to_replicas(vid, 
                                      false,  // async  
                                      false); // untracked
      }
      else {
        std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
        assert(vidowner.first);
        rmi.remote_request(vidowner.second,
                           &distributed_graph<VertexData,EdgeData>::set_vertex_data,
                           vid,
                           vdata);
      }
    }

    /**
     * Sets the data on the edge source->target
     * If the edge is not on this fragment, the request is sent
     * to a remote machine. This modification is performed asynchronously.
     */
    void set_edge_data_async(vertex_id_type source, 
                             vertex_id_type target, 
                             const EdgeData edata) {
      set_edge_data_from_pair(source, target, edata, true);
    }

    /**
     * Sets the data on the edge eid
     * If the edge is not on this fragment, the request is sent
     * to a remote machine. This modification is performed asynchronously.
     */
    void set_edge_data_async(edge_id_type eid, const EdgeData edata){
      set_edge_data_from_eid(eid, edata, true);
    }

    /**
     * Sets the data on the vertex vid.
     * If the vertex is not on this fragment, the request is sent
     * to a remote machine. This modification is performed asynchronously.
     */
    void set_vertex_data_async(vertex_id_type vid, const VertexData vdata){
      if (global_vid_in_local_fragment(vid)) {
        vertex_data(vid) = vdata;      
        vertex_is_modified(vid);
        push_owned_vertex_to_replicas(vid, 
                                      true,  // async  
                                      true); // untracked

      }
      else {
        std::pair<bool, procid_t> vidowner = 
          globalvid2owner.get_cached(vid);
        assert(vidowner.first);
        rmi.remote_call(vidowner.second,
                        &distributed_graph<VertexData,EdgeData>::set_vertex_data_async,
                        vid,
                        vdata);
      }
    }


    /** 
     * Get (and cache) the number of colors
     */
    size_t recompute_num_colors() {
      vertex_color_type max_color(0);
      for(size_t i = 0; i < localstore.num_vertices(); ++i) {
        max_color = std::max(max_color, localstore.color(i));
      }
      std::vector<vertex_color_type> proc2colors(rmi.numprocs());
      proc2colors[rmi.procid()] = max_color +1;
      rmi.all_gather(proc2colors);
      numcolors = 0;
      for(size_t i = 0; i < proc2colors.size(); ++i)  
        numcolors += proc2colors[i];
      return numcolors;
    }

    /**
     * This function should only be called after the number of colors
     * has be computed
     */
    size_t num_colors() const {
      return numcolors;
    }

    /**
     * Gets a reference to the color on vertex vid.
     * Assertion failure if vid is not on this machine.
     */
    vertex_color_type& color(vertex_id_type vid) {
      assert(global_vid_in_local_fragment(vid));
      return localstore.color(global2localvid[vid]);
    }

    /**
     * Gets a constant reference to the color on vertex vid.
     * Assertion failure if vid is not on this machine.
     */
    const vertex_color_type& color(vertex_id_type vid) const {
      assert(global_vid_in_local_fragment(vid));
      return localstore.color(global2localvid[vid]);
    }

    /**
     * Gets the color on vertex vid.
     * If the vertex is not on this fragment, the request is sent
     * to a remote machine.
     */
    vertex_color_type get_color(vertex_id_type vid) const{
      if (global_vid_in_local_fragment(vid)) {
        return localstore.color(globalvid_to_localvid(vid));
      }
      else {
        std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
        assert(vidowner.first);
        return rmi.remote_request(vidowner.second,
                                  &distributed_graph<VertexData,EdgeData>::get_color,
                                  vid);
      }
    }

    /**
     * Sets the color on vertex vid.
     * If the vertex is not on this fragment, the request is sent
     * to a remote machine. This operation is performed synchronously.
     * It will wait for the remote machine to complete the modification before
     * returning control.
     */
    void set_color(vertex_id_type vid, vertex_color_type color) {
      if (global_vid_in_local_fragment(vid)) {
        localstore.color(global2localvid[vid]) = color;
      }
      else {
        std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
        assert(vidowner.first);
        return rmi.remote_request(vidowner.second,
                                  &distributed_graph<VertexData,EdgeData>::set_color,
                                  vid,
                                  color);
      }
    }

    /**
     * Sets the color on vertex vid.
     * If the vertex is not on this fragment, the request is sent
     * to a remote machine. This operation is performed asynchronously.
     */
    void set_color_async(vertex_id_type vid, vertex_color_type color) const{
      if (global_vid_in_local_fragment(vid)) {
        localstore.color(global2localvid[vid]) = color;
      }
      else {
        std::pair<bool, procid_t> vidowner = globalvid2owner.get_cached(vid);
        assert(vidowner.first);
        return rmi.remote_call(vidowner.second,
                               &distributed_graph<VertexData,EdgeData>::set_color_async,
                               vid,
                               color);
      }
    }

    /**
       Collects all the vertex data onto one machine.
       The target machine will be returned a vector containing all the vertex data
       while all machines will returned an empty vector.
       All machines must call this function together and must have the same
       value for 'targetmachine'
       \warning This is only useful when the graph is small enough that
       all the vertex data fits in memory.
    */
    std::vector<VertexData> collect_vertices(procid_t targetmachine) {
      std::vector<VertexData> ret;
      typename std::vector<std::map<vertex_id_type, VertexData> > gather(rmi.numprocs());
      foreach(vertex_id_type vid, owned_vertices()) {
        gather[rmi.procid()][vid]  = vertex_data(vid);
      }
    
      rmi.gather(gather, targetmachine);
    
      if (rmi.procid() == 0) {
        ret.resize(num_vertices());
        for (size_t i = 0;i < gather.size(); ++i) {
          typename std::map<vertex_id_type, VertexData>::const_iterator iter = gather[i].begin();
          while (iter != gather[i].end()) {
            ret[iter->first] = iter->second;
            ++iter;
          }
        }
      }
      return ret;
    }



    /**
     * Like collects a subset of vertex data onto the current machine, 
     * 'vids' is a list of vertex IDs to collect. Result will be returned in the return value
     * as an associative map from vertex ID to vertex Data.
     * This function can be called from one machine.
     */
    std::map<vertex_id_type, VertexData> collect_vertex_subset_one_way(const std::vector<vertex_id_type> &vids) {
      dc_impl::stored_ret_type ret(rmi.numprocs());
      for (size_t i = 0;i < rmi.numprocs(); ++i) {
        rmi.remote_call(i,
                        &distributed_graph<VertexData,EdgeData>::send_vertex_subset,
                        vids,
                        rmi.procid(),
                        reinterpret_cast<size_t>(&ret));
      }
      ret.wait();
    
      std::map<vertex_id_type, VertexData> final_ret;
      // unpack
      std::map<procid_t, dc_impl::blob>::iterator iter = ret.val.begin();
      while (iter != ret.val.end()) {
        // read the blob
        boost::iostreams::stream<boost::iostreams::array_source> retstrm(iter->second.c, iter->second.len); 
        iarchive iarc(retstrm);
        std::map<vertex_id_type, VertexData> v_to_data;
        iarc >> v_to_data;
        iter->second.free();
      
      
        // insert into final_ret
        typename std::map<vertex_id_type, VertexData>::const_iterator miter = v_to_data.begin();
        while (miter != v_to_data.end()) {
          final_ret[miter->first] = miter->second;
          ++miter;
        }
        ++iter;
      }
      return final_ret;
    }
  
    // synchronzation calls. These are called from the ghost side
    // to synchronize against the owner.
  
    /**
     * synchronize the data on vertex with global id vid
     * vid must be a ghost
     */
    void synchronize_vertex(vertex_id_type vid, bool async = false);
  
    /**
     * synchronize the data on edge with global id eid
     * target of edge must be a ghost
     */

    void synchronize_edge(edge_id_type eid, bool async = false);
  
    /**
     * synchronize the entire scope for vertex vid
     * vid must be owned by the current machine. 
     * This function synchronizes all ghost data on this scope
     * with their owners.
     */
    void synchronize_scope(vertex_id_type vid, bool async = false);   
  
    void allocate_scope_callbacks() {
      scope_callbacks.resize(local2globalvid.size());
      for (size_t i = 0;i < local2globalvid.size(); ++i) {
        scope_callbacks[i].callback = NULL;
      }
    }
    /**
       Synchonizes the ghost data of the scope around vid.
       When done, the callback is issued with the tag as the first argument.
       Allocate vertex callbacks must be called first. There can be at most one
       callback associated with any vertex. an assertion failure will be thrown otherwise.
       vid must be on the boundary of the fragment. An assertion will be thrown otherwise.
    */
    void async_synchronize_scope_callback(vertex_id_type vid, 
                                          boost::function<void (void)>);

    /**
     * Waits for all asynchronous data synchronizations to complete
     */
    void wait_for_all_async_syncs();
    /**
       Synchronize all ghost vertices
    */
    void synchronize_all_vertices(bool async = false);
  
    /**
       Synchronize all ghost edges
    */
    void synchronize_all_edges(bool async = false);
  
    /**
       Synchronize all ghost scopes.
    */
    void synchronize_all_scopes(bool async = false);

  
    /** Called from the owner side to synchronize the owner against ghosts. 
     * Pushes an owned vertex to all ghosts.
     */
    void push_owned_vertex_to_replicas(vertex_id_type vid, bool async = false, bool untracked = false);
  
    /** Called from the owner side to synchronize the owner against ghosts. 
     * Pushes an owned edge to all ghosts.
     */
    void push_owned_edge_to_replicas(edge_id_type eid, bool async = false, bool untracked = false);
  
    /** Called from the owner side to synchronize the owner against ghosts. 
     * Pushes an owned scope to all ghosts.
     */
    void push_owned_scope_to_replicas(vertex_id_type vid, 
                                      bool modified_only, 
                                      bool clear_modified, 
                                      bool async = false, 
                                      bool untracked = false);

    /** Updates all ghosted vertices. */
    void push_all_owned_vertices_to_replicas();
  
    /** Updates all ghosted edges . */
    void push_all_owned_edges_to_replicas();
  
    /** Waits for all asynchronous push requests to complete */
    void wait_for_all_async_pushes();
  public:
  
    // extra types
    template <typename DataType>
    struct conditional_store{
      bool hasdata;
      DataType data;
      void save(oarchive &oarc) const {
        oarc << hasdata;
        if (hasdata) oarc << data;
      }
      void load(iarchive &iarc) {
        iarc >>  hasdata;
        if (hasdata) iarc >> data;
      }
    };
  
    typedef conditional_store<std::pair<VertexData, uint64_t> >  vertex_conditional_store;
    typedef conditional_store<std::pair<EdgeData, uint64_t> >  edge_conditional_store;

  
    struct block_synchronize_request2 {
      std::vector<vertex_id_type> vid;
      std::vector<uint64_t > vidversion;
      std::vector<vertex_conditional_store> vstore;
      std::vector<std::pair<vertex_id_type, vertex_id_type> > srcdest;
      std::vector<uint64_t > edgeversion;
      std::vector<edge_conditional_store> estore;

      void clear() {
        vid.clear();
        vidversion.clear();
        vstore.clear();
        srcdest.clear();
        edgeversion.clear();
        estore.clear();
      }
 
      void save(oarchive &oarc) const{
        oarc << vid << vidversion << vstore
             << srcdest << edgeversion << estore;
      }

      void load(iarchive &iarc) {
        iarc >> vid >> vidversion >> vstore
             >> srcdest >> edgeversion >> estore;
      }
    };

    typedef std::pair<block_synchronize_request2, 
                      typename std::vector<vertex_id_type>::iterator>
    request_veciter_pair_type;


    /// RMI object
    mutable dc_dist_object<distributed_graph<VertexData, EdgeData> > rmi;

  private:
  
    std::string indexfilename;
    std::vector<std::vector<size_t> > atompartitions;
   
    /** Protects structural modifications of the graph.
     * Modifications to the data store and to the local<->global mappings
     * must lock this.
     */
    mutex alldatalock;
  
    /// stores the local fragment of the graph
    graph_local_store_type localstore;


    typedef boost::unordered_map<vertex_id_type, vertex_id_type> 
    global2localvid_type;

    /** all the mappings requried to move from global to local vid/eids
     *  We only store mappings if the vid/eid is in the local fragment
     */
    global2localvid_type global2localvid;

    std::vector<vertex_id_type> local2globalvid;

    /// collection of vertices I own
    std::vector<vertex_id_type> ownedvertices;

    /// collection of ghost vertices
    std::vector<vertex_id_type> ghostvertices;

    /// collection of owned vertices which have a ghost neighbor
    boost::unordered_set<vertex_id_type> boundaryscopesset;
    std::vector<vertex_id_type> boundaryscopes;
  
    std::vector<std::vector<procid_t> > localvid2ghostedprocs;
    std::vector<procid_t> cur_proc_vector;  // vector containing only 1 element. the current proc
  
    /** To avoid requiring O(V) storage on each maching, the 
     * global_vid -> owner mapping cannot be stored in its entirely locally
     * instead, we store it in a DHT. \see globaleid2owner
     */
    lazy_dht<vertex_id_type, procid_t> globalvid2owner;
  
    /** This provides a fast mapping from the local vids in the fragment
     * to its owner. Since this operation is quite frequently needed.
     */
    std::vector<procid_t> localvid2owner;
  
    /**
     * The number of vertices and edges in the entire graph so far.
     * Currently only consistent on machine 0 since machine 0 manages 
     * the allocation of global VIDs and local VIDs.
     */
    size_t numglobalverts, numglobaledges, numcolors;

    dc_impl::reply_ret_type pending_async_updates;
    dc_impl::reply_ret_type pending_push_updates;
  
  
    struct async_scope_callback {
      atomic<unsigned short> counter;
      boost::function<void (void)> callback;
    };
    std::vector<async_scope_callback> scope_callbacks;

    metrics graph_metrics;

    /**
     * Returns true if the global vid is in the local fragment
     */
    bool global_vid_in_local_fragment(vertex_id_type globalvid) const{
      // easiest way to check is to see if it is in the global2localvid mapping
      return global2localvid.find(globalvid) != global2localvid.end();
    }
  

  
    /**
     * From the atoms listed in the atom index file, construct the local store
     * using all the atoms in the current partition.
     */
    void construct_local_fragment(const atom_index_file &atomindex,
                                  std::vector<std::vector<size_t> > partitiontoatom,
                                  size_t curpartition,
                                  bool do_not_load_data,
                                  bool do_not_mmap) {
      timer loadtimer;
      loadtimer.start();
      // first make a map mapping atoms to machines
      // we will need this later
      std::vector<procid_t> atom2machine;
      for (size_t i = 0 ;i< partitiontoatom.size(); ++i) {
        for (size_t j = 0 ; j < partitiontoatom[i].size(); ++j) {
          if (atom2machine.size() <= partitiontoatom[i][j]) 
            atom2machine.resize(partitiontoatom[i][j] + 1);
          atom2machine[partitiontoatom[i][j]] = i;
        }
      }

    
    
      // the atomfiles for the local fragment
      std::vector<disk_atom*> atomfiles;
      // for convenience take a reference to the list of atoms in this partition
      std::vector<size_t>& atoms_in_curpart = partitiontoatom[curpartition];
      dense_bitset atoms_in_curpart_set(atomindex.atoms.size()); // make a set vertion for quick lookup
      atoms_in_curpart_set.clear();
      std::vector<std::vector<vertex_id_type> > vertices_in_atom;
      // create the atom file readers
      // open the atoms
      atomfiles.resize(atoms_in_curpart.size());
      logstream(LOG_INFO) << "Atoms on this machine: " << atoms_in_curpart.size() << std::endl;
      // open the atoms we are assigned to
      vertices_in_atom.resize(atoms_in_curpart.size());
      for (int i = 0;i < (int)(atoms_in_curpart.size()); ++i) {
        atoms_in_curpart_set.set_bit(atoms_in_curpart[i]);
        atomfiles[i] = new disk_atom(atomindex.atoms[atoms_in_curpart[i]].file, 
                                     atoms_in_curpart[i]);
        atomfiles[i]->precache();
        vertices_in_atom[i] = atomfiles[i]->enumerate_vertices();
      }
    
    
      logger(LOG_INFO, "Generating mappings");
      /******************* get all the owned vertices  *************************/
      // Takes all the open atoms, and returns the union of all the
      // vertices in each atom (owned and ghosts). 
      //
      // cat all the globalvids and globaleids into a single big list
      // and sort it. Make sure that my owned vertices come first before the ghost vertices
      //
      // globalvid_notowned_zip[i].second is the vertex ID
      // and globalvid_notowned_zip[i].first == true if we do not own the vertex
      //
      // sorting this array will therefore put all the owned vertices at the start
      std::vector<std::pair<bool, vertex_id_type> > globalvid_notowned_zip;
      for (size_t i = 0;i < atomfiles.size(); ++i) {
        std::vector<vertex_id_type>& vertices = vertices_in_atom[i];
        for (size_t j = 0;j < vertices.size() ; ++j) {
          uint16_t owneratom;
          ASSERT_TRUE(atomfiles[i]->get_vertex(vertices[j], owneratom));
          globalvid_notowned_zip.push_back(std::make_pair(atom2machine[owneratom] != rmi.procid(),
                                                          vertices[j]));
        }
      }
    
      // Find only unique occurances of each vertex, by sorting, unique,
      // and resize
      std::sort(globalvid_notowned_zip.begin(), 
                globalvid_notowned_zip.end());
      typename std::vector<std::pair<bool, vertex_id_type> >::
        iterator uviter = 
        std::unique(globalvid_notowned_zip.begin(), 
                    globalvid_notowned_zip.end());
              
      globalvid_notowned_zip.resize(uviter - globalvid_notowned_zip.begin());
      std::sort(globalvid_notowned_zip.begin(), globalvid_notowned_zip.end());
    
    
      
      local2globalvid.resize(globalvid_notowned_zip.size());
      for (size_t i = 0; i < globalvid_notowned_zip.size(); ++i) {
        // this is a sanity checks that all the owned vertices come first
        if (i > 0 && globalvid_notowned_zip[i].first == false) ASSERT_EQ(globalvid_notowned_zip[i-1].first,  false); 
        local2globalvid[i] = globalvid_notowned_zip[i].second;
      } 
    
      //construct the reverse maps
      for (size_t i = 0; i < local2globalvid.size(); ++i) {
        global2localvid[local2globalvid[i]] = i;
      }
      global2localvid.rehash(2 * global2localvid.size());




      logger(LOG_INFO, "Counting Edges");
      /****** figure out how many edges I need to instantiate from each atom ****/
      // This will let me isntantiate the edges in parallel later

      size_t nedges_to_create = 0;
      std::vector<size_t> acc_edges_created_in_this_atom(atomfiles.size(), 0);
#pragma omp parallel for reduction(+ : nedges_to_create)
      for (int i = 0;i < (int)(atomfiles.size()); ++i) {
        std::vector<vertex_id_type> vertices = vertices_in_atom[i];
        foreach(vertex_id_type dest, vertices) {
          uint16_t destowneratom;
          ASSERT_TRUE(atomfiles[i]->get_vertex(dest, destowneratom));
          // the atom owns the edge if the target is within the atom
          bool newedge = (destowneratom == atomfiles[i]->atom_id());
        
          // otherwise this must be an edge out connecting to the ghost of an atom
          // Now if I also own the target atom then we do not need to instantiate the edge
          // since the target atom will instantiate it.
          //
          // on the other hand, we do need to instantiate the edge if we do not
          // own the target since this means that it is a true ghosted edge
          newedge = newedge || (!atoms_in_curpart_set.get(destowneratom)); 
          if (newedge) {
            size_t ne = atomfiles[i]->get_in_vertices(dest).size();
            nedges_to_create += ne;
            acc_edges_created_in_this_atom[i] += ne;
          }
        }
      }
      // Now we compute the prefix sum of the acc_edges_created_in_this_atom 
      // which will give us the first edge id to be instanted by each atom file
      std::vector<size_t> atom_file_edge_first_id(atomfiles.size(), 0);
   
      for (size_t i = 1;i < atomfiles.size(); ++i) {
        acc_edges_created_in_this_atom[i] += acc_edges_created_in_this_atom[i - 1];
        atom_file_edge_first_id[i] = acc_edges_created_in_this_atom[i - 1];
      } 
    
    
      logstream(LOG_INFO) << "Creating " << nedges_to_create << " edges locally." << std::endl;


      // open the local store
      logger(LOG_INFO, "Creating local store");
      // now lets construct the graph structure
      localstore.create_store(local2globalvid.size(), nedges_to_create,
                              "vdata." + tostr(curpartition),
                              "edata." + tostr(curpartition),
                              do_not_mmap);
      localstore.zero_all();
      // create a course grained lock for vertex IDs
      // edge insertions should not touch the same vertex at the same time
      std::vector<mutex> hashvlock;
      hashvlock.resize((1 << 14));
    
#pragma omp parallel for
      for (int i = 0;i < (int)atomfiles.size(); ++i) {
        std::cerr << ".";
        std::cerr.flush();
      
        size_t nextedgeid = atom_file_edge_first_id[i];

        // iterate through all the edges in this atom
        // in the same order as we did the counting
        foreach(vertex_id_type dest, vertices_in_atom[i]) {
          uint16_t destowneratom;
          ASSERT_TRUE(atomfiles[i]->get_vertex(dest, destowneratom));
          foreach(vertex_id_type src, atomfiles[i]->get_in_vertices(dest)) {

            std::pair<vertex_id_type, vertex_id_type> globaledge = std::make_pair(src, dest);
          
            bool newedge = (destowneratom == atomfiles[i]->atom_id());
            newedge = newedge || (!atoms_in_curpart_set.get(destowneratom)); 
          
            if (newedge == false) continue;
            edge_id_type eid = nextedgeid;
            nextedgeid++;
          
            std::pair<vertex_id_type, vertex_id_type> localedge = global_edge_to_local_edge(globaledge);
            size_t k1 = localedge.first & ((1 << 14) - 1);
            size_t k2 = localedge.second & ((1 << 14) - 1);
            if (k1 < k2) {
              hashvlock[k1].lock();
              hashvlock[k2].lock();
            }
            else if (k2 < k1) {
              hashvlock[k2].lock();
              hashvlock[k1].lock();
            }
            else {
              hashvlock[k1].lock();
            }
            localstore.add_edge(eid, localedge.first, localedge.second);
            if (k1 != k2) {
              hashvlock[k1].unlock();
              hashvlock[k2].unlock();
            }
            else {
              hashvlock[k1].unlock();
            }
          }
        }
      }
    

      std::cerr << std::endl;
    
      logstream(LOG_INFO) << "Local structure creation complete." << std::endl;
    
      // construct localvid2owner
      localvid2owner.resize(local2globalvid.size());
    
#pragma omp parallel for
      for (int i = 0;i < (int)atomfiles.size(); ++i) {
        // set the color and localvid2owner mappings
        foreach (vertex_id_type globalvid, vertices_in_atom[i]) {
          // get the localvid
          vertex_id_type localvid = global2localvid[globalvid];
          uint16_t owneratom;
          ASSERT_TRUE(atomfiles[i]->get_vertex(globalvid, owneratom));
          localvid2owner[localvid] = atom2machine[owneratom];
          localstore.color(localvid) = atomfiles[i]->get_color(globalvid);
          // if I own this vertex, set the global ownership to me
          if (localvid2owner[localvid] == rmi.procid()) {
            globalvid2owner.set(globalvid, rmi.procid());
          }
        }
      }
    
      // check for contiguousness of localvid
      for (size_t i = 1; i < local2globalvid.size(); ++i) {
        if (localvid2owner[i] == rmi.procid()) ASSERT_EQ(localvid2owner[i-1], rmi.procid());
      }
      logstream(LOG_INFO) << "vid -> Owner DHT set complete" << std::endl;
      logstream(LOG_INFO) << "Constructing auxiliary datastructures..." << std::endl;

      // fill the ownedvertices list
      // create thread local versions of ownedvertices, ghostvertices and boundaryscopeset
    
      std::vector<std::vector<vertex_id_type> > __ownedvertices__(omp_get_max_threads());
      std::vector<std::vector<vertex_id_type> > __ghostvertices__(omp_get_max_threads());
      std::vector<boost::unordered_set<vertex_id_type> > __boundaryscopesset__(omp_get_max_threads());
      std::vector<dense_bitset> __ghostownerset__(omp_get_max_threads());
    
      for (size_t i = 0;i < __ghostownerset__.size(); ++i) {
        __ghostownerset__[i].resize(rmi.numprocs());
      }
      spinlock ghostlock, boundarylock;
      // construct the vid->replica mapping
      // for efficiency reasons this only contains the maps for ghost vertices
      // if the corresponding vector is empty, the only replica is the current machine
      // use localvid_to_replicas() instead to access this since it provides 
      // more consistent behavior (it checks. if the array is empty, it will return
      // a vector containing only the current procid)
    
      localvid2ghostedprocs.resize(localvid2owner.size());
    
#pragma omp parallel for
      for (long i = 0;i < (long)localvid2owner.size(); ++i) {
        int thrnum = omp_get_thread_num();
        if (localvid2owner[i] == rmi.procid()) {
          __ownedvertices__[thrnum].push_back(local2globalvid[i]);
          // loop through the neighbors and figure out who else might
          // have a ghost of me. Fill the vertex2ghostedprocs vector
          // those who have a ghost of me are the owners of my ghost vertices
          __ghostownerset__[thrnum].clear();
          __ghostownerset__[thrnum].set_bit_unsync(rmi.procid());  // must always have the current proc
          foreach(edge_id_type ineid, localstore.in_edge_ids(i)) {
            vertex_id_type localinvid = localstore.source(ineid);
            if (localvid2owner[localinvid] != rmi.procid()) {
              __ghostownerset__[thrnum].set_bit_unsync(localvid2owner[localinvid]);
            }
          }
          foreach(edge_id_type outeid, localstore.out_edge_ids(i)) {
            vertex_id_type localoutvid = localstore.target(outeid);
            if (localvid2owner[localoutvid] != rmi.procid()) {
              __ghostownerset__[thrnum].set_bit_unsync(localvid2owner[localoutvid]);
            }
          }
          uint32_t b;
          ASSERT_TRUE(__ghostownerset__[thrnum].first_bit(b));
          do {
            localvid2ghostedprocs[i].push_back(b);
          } while(__ghostownerset__[thrnum].next_bit(b));
        }
        else {
          __ghostvertices__[thrnum].push_back(local2globalvid[i]);
          // if any of my neighbors are not ghosts, they are a boundary scope
          foreach(edge_id_type ineid, localstore.in_edge_ids(i)) {
            vertex_id_type localinvid = localstore.source(ineid);
            if (localvid2owner[localinvid] == rmi.procid()) {
              __boundaryscopesset__[thrnum].insert(local2globalvid[localinvid]);
            }
          }
          foreach(edge_id_type outeid, localstore.out_edge_ids(i)) {
            vertex_id_type localoutvid = localstore.target(outeid);
            if (localvid2owner[localoutvid] == rmi.procid()) {
              __boundaryscopesset__[thrnum].insert(local2globalvid[localoutvid]);
            }
          }
        }
      }
    
    
      // join all the thread local datastructures
      for (size_t i = 0;i < __ghostvertices__.size(); ++i) {
        std::copy(__ownedvertices__[i].begin(), __ownedvertices__[i].end(), 
                  std::back_inserter(ownedvertices));
        std::copy(__ghostvertices__[i].begin(), __ghostvertices__[i].end(), 
                  std::back_inserter(ghostvertices));
        boundaryscopesset.insert(__boundaryscopesset__[i].begin(), __boundaryscopesset__[i].end());
      }
    
      std::sort(ownedvertices.begin(), ownedvertices.end());
      std::sort(ghostvertices.begin(), ghostvertices.end());
      std::copy(boundaryscopesset.begin(), boundaryscopesset.end(),
                std::back_inserter(boundaryscopes));
      __ownedvertices__.clear();
      __ghostvertices__.clear();
      __boundaryscopesset__.clear();
    
    
      if (do_not_load_data == false) {
        logger(LOG_INFO, "Loading data");
        // done! structure constructed!  now for the data!  
      
#pragma omp parallel for
        for (int i = 0;i < (int)(atomfiles.size()); ++i) {
          std::cerr << ".";
          std::cerr.flush();
        
          // loop through the vertices
          foreach(vertex_id_type globalvid, vertices_in_atom[i]) {
            uint16_t owneratom;
            ASSERT_TRUE(atomfiles[i]->get_vertex(globalvid, owneratom));

            // if the atomfile contains the data.
            if (owneratom == atomfiles[i]->atom_id()) {
              size_t localvid = global2localvid[globalvid];
              ASSERT_TRUE(atomfiles[i]->get_vertex<VertexData>(globalvid, 
                                                               owneratom, 
                                                               localstore.vertex_data(localvid)));
              localstore.set_vertex_version(localvid, 1);
            }
          }

    
          size_t nextedgeid = atom_file_edge_first_id[i];
          // loop through the edges
        
          // iterate through all the edges in this atom
          // in the same order as we did the counting
          foreach(vertex_id_type dest, vertices_in_atom[i]) {
            uint16_t destowneratom;
            ASSERT_TRUE(atomfiles[i]->get_vertex(dest, destowneratom));
            foreach(vertex_id_type src, atomfiles[i]->get_in_vertices(dest)) {

              bool newedge = (destowneratom == atomfiles[i]->atom_id());
              newedge = newedge || (!atoms_in_curpart_set.get(destowneratom)); 

              if (newedge == false) {
                continue;
              }
              edge_id_type eid = nextedgeid;
              nextedgeid++;    
              // get the local edge
              EdgeData temp;
              if (atomfiles[i]->get_edge(src, dest, temp)) {
                localstore.edge_data(eid) = temp;
                localstore.set_edge_version(eid, 1);
              }
            }
          }
        }



        rmi.barrier();
        logger(LOG_INFO, "Synchronizing ghost data...");
        // shuffle for all the ghost data
        push_all_owned_vertices_to_replicas();
        rmi.dc().full_barrier();
        logger(LOG_INFO, "vertices synchronized.");
      
        push_all_owned_edges_to_replicas();
        rmi.dc().full_barrier();
        logger(LOG_INFO, "edges synchronized.");

        logger(LOG_INFO, "Synchronization complete.");
        rmi.dc().barrier();
        logger(LOG_INFO, "Performing data verification.");
        for (size_t i = 0;i < localstore.num_vertices(); ++i) {
          ASSERT_EQ(localstore.vertex_version(i), 1);
        }
        for (size_t i = 0;i < localstore.num_edges(); ++i) {
          ASSERT_EQ(localstore.edge_version(i), 1);
        }
      }
      // flush the store
      logger(LOG_INFO, "Finalize");
      localstore.finalize();
      logger(LOG_INFO, "Flush");
      localstore.flush();
      if (do_not_mmap == false) {
        logger(LOG_INFO, "Prefetch computation"); 
        localstore.compute_minimal_prefetch();
      }
      logger(LOG_INFO, "Load complete.");
      rmi.comm_barrier();
      std::cout << "Load complete in " << loadtimer.current_time() << std::endl;
    }

    vertex_conditional_store get_vertex_if_version_less_than(vertex_id_type vid, 
                                                             uint64_t vertexversion,
                                                             vertex_conditional_store &vdata);
                                                       
    edge_conditional_store get_edge_if_version_less_than2(vertex_id_type source, 
                                                          vertex_id_type target, 
                                                          uint64_t edgeversion,
                                                          edge_conditional_store &edata);
  
  
    void async_get_vertex_if_version_less_than(procid_t srcproc, 
                                               vertex_id_type vid, 
                                               uint64_t vertexversion,
                                               vertex_conditional_store &vdata);
                                           
    void async_get_edge_if_version_less_than2(procid_t srcproc, 
                                              vertex_id_type source, 
                                              vertex_id_type target, 
                                              uint64_t edgeversion,
                                              edge_conditional_store &edata);
                                   
                                   
  
    void reply_vertex_data_and_version(vertex_id_type vid,
                                       vertex_conditional_store &estore);
  
  
    void reply_edge_data_and_version(edge_id_type eid,
                                     edge_conditional_store &estore);
  
  
    void reply_edge_data_and_version2(vertex_id_type source, 
                                      vertex_id_type target, 
                                      edge_conditional_store &estore);
  
  
    std::pair<vertex_id_type, vertex_id_type> 
    local_edge_to_global_edge(std::pair<vertex_id_type, vertex_id_type> e) const{
      return std::make_pair(local2globalvid[e.first], local2globalvid[e.second]);
    }
  
    std::pair<vertex_id_type, vertex_id_type> 
    global_edge_to_local_edge(std::pair<vertex_id_type, vertex_id_type> e) const{
      typename global2localvid_type::const_iterator iter1 = global2localvid.find(e.first);
      typename global2localvid_type::const_iterator iter2 = global2localvid.find(e.second);
      assert(iter1 != global2localvid.end());
      assert(iter2 != global2localvid.end()); 
      return std::make_pair(iter1->second, iter2->second);
    }
  
 
    /**
       Constructs the request set for a scope synchronization
       the request type is a little strange , but it is for efficiency reasons.
       the second component of the pair can be ignored
    */
    void synchronize_scope_construct_req(vertex_id_type vid, 
                                         std::map<procid_t, request_veciter_pair_type > &requests);   
 
 
    void update_vertex_data_and_version(vertex_id_type vid,
                                        vertex_conditional_store &estore);
  
  
    void update_edge_data_and_version(edge_id_type eid,
                                      edge_conditional_store &estore);
  
  
    void update_edge_data_and_version2(vertex_id_type source, 
                                       vertex_id_type target, 
                                       edge_conditional_store &estore);
  
  
    block_synchronize_request2& get_alot2(block_synchronize_request2 &request);
  
    void async_get_alot2(procid_t srcproc, 
                         block_synchronize_request2 &request, 
                         size_t replytarget, size_t tag);
  
    void update_alot2(block_synchronize_request2 &request) ;
  
    void reply_alot2(block_synchronize_request2 &request, size_t replytarget, size_t tag);


    void update_vertex_data_and_version_and_reply(
                                                  vertex_id_type vid, 
                                                  distributed_graph<VertexData, EdgeData>::vertex_conditional_store &vstore,
                                                  procid_t srcproc,
                                                  size_t reply);


    void update_edge_data_and_version_and_reply2(
                                                 vertex_id_type source, 
                                                 vertex_id_type target, 
                                                 edge_conditional_store &estore,
                                                 procid_t srcproc, size_t reply);
                          
    friend class graph_lock<distributed_graph<VertexData, EdgeData> >;
 
  
  
    /***********************************************************************
     *    End of Synchronization Ops. 
     *   Start of Miscellenous stuff
     ***********************************************************************/
  
  public: 
  
    void fill_metrics() {
      std::map<std::string, size_t> ret = rmi.gather_statistics();
    
      std::vector<size_t> procpartitionsize(rmi.numprocs(), 0);
      std::vector<size_t> procghosts(rmi.numprocs(), 0);
      procpartitionsize[rmi.procid()] = local_vertices();
      procghosts[rmi.procid()] = ghost_vertices().size();
    
      rmi.gather(procpartitionsize, 0);
      rmi.gather(procghosts, 0);
    
      if (rmi.procid() == 0) {
        graph_metrics.set("num_vertices", num_vertices(), INTEGER);
        graph_metrics.set("num_edges", num_edges(), INTEGER);
        graph_metrics.set("total_calls_sent", ret["total_calls_sent"], INTEGER);
        graph_metrics.set("total_bytes_sent", ret["total_bytes_sent"], INTEGER);
      
        for(int i=0; i<rmi.numprocs(); i++) {
          graph_metrics.set_vector_entry("local_part_size", i, procpartitionsize[i]);
          graph_metrics.set_vector_entry("ghosts_size", i, procghosts[i]);
        
        }
      } 
    }

    metrics get_metrics() {
      return graph_metrics;
    }


    void reset_metrics() {
      graph_metrics.clear();
    }

    void report_metrics(imetrics_reporter &reporter) {
      graph_metrics.report(reporter);
    }


    /**
       Receiving side of the collect_vertex_subset_one_way() function. Not to be used directly.
    */
    void send_vertex_subset(const std::vector<vertex_id_type> &vids, procid_t srcproc, size_t id) {
      std::map<vertex_id_type, VertexData> ret;
      for (size_t i = 0;i < vids.size(); ++i) {
        if (is_owned(vids[i])) {
          ret[vids[i]] = vertex_data(vids[i]);
        }
      }
      boost::iostreams::stream<resizing_array_sink> retstrm(128); 
      oarchive oarc(retstrm);
      oarc << ret;
      retstrm.flush();
    
      rmi.dc().remote_call(srcproc, stored_increment_counter, id, 
                           dc_impl::blob(retstrm->str, retstrm->len));
      free(retstrm->str);
    }
  
  
  }; // End of class distributed Graph



























  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ////                                                                 ////
  ////                      Implementation                             ////
  ////                                                                 ////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////



  ///=====================================================================>
  // Synchronization 

  // #define FROM_DISTRIBUTED_GRAPH_INCLUDE
  // #include <graphlab/distributed2/graph/distributed_graph_synchronization.hpp>
  // #undef FROM_DISTRIBUTED_GRAPH_INCLUDE





  /**
   * synchronize the data on vertex with global id vid
   * vid must be a ghost
   */
  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  synchronize_vertex(vertex_id_type vid, bool async) {
    vertex_id_type localvid = global2localvid[vid];
    if (is_ghost(vid)) {
      vertex_conditional_store out;
      out.hasdata = localstore.vertex_modified(localvid);
      if (out.hasdata) {
        localstore.set_vertex_modified(localvid, false);
        out.data.first = localstore.vertex_data(localvid);
      }
      if (async == false) {
        vertex_conditional_store v;
        v = rmi.remote_request(localvid2owner[localvid],
                               &distributed_graph<VertexData, EdgeData>::
                               get_vertex_if_version_less_than,
                               vid,
                               localstore.vertex_version(localvid),
                               out);
        if (v.hasdata) {
          update_vertex_data_and_version(vid, v);
        }
      } else {
        pending_async_updates.flag.inc();
        rmi.remote_call(localvid2owner[localvid],
                        &distributed_graph<VertexData, EdgeData>::
                        async_get_vertex_if_version_less_than,
                        rmi.procid(),
                        vid,
                        localstore.vertex_version(localvid),
                        out);
      }
    }
  } // end of sycnhronize vertex


  /**
   * synchronize the data on edge with global id eid
   * target of edge must be a ghost
   */
  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  synchronize_edge(edge_id_type eid, bool async) {
    vertex_id_type localtargetvid = localstore.target(eid);
    vertex_id_type targetvid = local2globalvid[localtargetvid];
    vertex_id_type localsourcevid = localstore.source(eid);
    vertex_id_type sourcevid = local2globalvid[localsourcevid];

    if (is_ghost(targetvid)) {

      edge_conditional_store out;
      out.hasdata = localstore.edge_modified(eid);
      if (out.hasdata) {
        localstore.set_edge_modified(eid, false);
        out.data.first = localstore.edge_data(eid);
      }
      if (async == false) {
        edge_conditional_store e = 
          rmi.remote_request(localvid2owner[localtargetvid],
                             &distributed_graph<VertexData, EdgeData>::
                             get_edge_if_version_less_than2,
                             sourcevid,
                             targetvid,
                             localstore.edge_version(eid),
                             out);
        if (e.hasdata) {
          update_edge_data_and_version2(sourcevid, targetvid, e);
        }
      } else {
        pending_async_updates.flag.inc();
        rmi.remote_call(localvid2owner[localtargetvid],
                        &distributed_graph<VertexData, EdgeData>::
                        async_get_edge_if_version_less_than2,
                        rmi.procid(),
                        sourcevid,
                        targetvid,
                        localstore.edge_version(eid),
                        out);
      }
    }
  } // end of synchronize edge





  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  synchronize_scope(vertex_id_type vid, bool async) {
    // construct he requests
    typedef request_veciter_pair_type pair_type;
    typedef std::map<procid_t, pair_type> map_type;
    map_type requests;
    synchronize_scope_construct_req(vid, requests);
  
    if (async) {
      // if asynchronous, the reply goes to pending_async_updates
      typename map_type::iterator iter;
      iter = requests.begin();
      size_t replytarget = reinterpret_cast<size_t>(&pending_async_updates);
      pending_async_updates.flag.inc(requests.size());

      while(iter != requests.end()) {
        rmi.remote_call(iter->first,
                        &distributed_graph<VertexData, EdgeData>::
                        async_get_alot2,
                        rmi.procid(),
                        iter->second.first,
                        replytarget,
                        0);
        ++iter;
      }
    } else {
      // otherwise we collect it into a local reply ret tye
      dc_impl::reply_ret_type reply(true, 0);
      typename map_type::iterator iter;
      iter = requests.begin();
      size_t replytarget = reinterpret_cast<size_t>(&reply);
      reply.flag.inc(requests.size());

      while(iter != requests.end()) {
        rmi.remote_call(iter->first,
                        &distributed_graph<VertexData, EdgeData>::
                        async_get_alot2,
                        rmi.procid(),
                        iter->second.first,
                        replytarget,
                        0);
        ++iter;
      }

      reply.wait();
    }  
  } // end of synchronize scope





  /**
   * Constructs the request for synchronizing the scope for vertex vid
   * vid must be owned by the current machine. 
   */
  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  synchronize_scope_construct_req(vertex_id_type vid, 
                                  std::map<procid_t, 
                                           request_veciter_pair_type > &requests) {
    ASSERT_FALSE(is_ghost(vid));
    if (boundaryscopesset.find(vid) == boundaryscopesset.end()) return;

    // now this is very annoying. A significant amount of code is identical here.
    // whether with edge canonical numbering on or not. But I cannot refactor it
    // easily because the types are different and I don't really want to introduce
    // templates here for something so trivial.
  
    vertex_id_type localvid = global2localvid[vid];
    requests.clear();
    // I should have all the in-edges. but I need the in vertices.
    // need to track the vertices added so I don't add duplicate vertices
    // if the vertex has both in-out edges to this vertex.
    // trick! vids are ordered!

    foreach(edge_id_type localineid, localstore.in_edge_ids(localvid)) {
      vertex_id_type localsourcevid = localstore.source(localineid);
      if (localvid_is_ghost(localsourcevid)) {
        // need to synchronize incoming vertex
        procid_t sourceowner = localvid2owner[localsourcevid];
        block_synchronize_request2 &req = requests[sourceowner].first;
        req.vid.push_back(local2globalvid[localsourcevid]);
        req.vidversion.push_back(localstore.vertex_version(localsourcevid));
        vertex_conditional_store vs;
        vs.hasdata = localstore.vertex_modified(localsourcevid);
        if (vs.hasdata) {
          localstore.set_vertex_modified(localsourcevid, false);
          vs.data.first = localstore.vertex_data(localsourcevid);
          vs.data.second = localstore.vertex_version(localsourcevid);
        }
        requests[sourceowner].second=req.vid.end();
        req.vstore.push_back(vs);
      }
    }
    // now for the out edges
    foreach(edge_id_type localouteid, localstore.out_edge_ids(localvid)) {
      vertex_id_type localtargetvid = localstore.target(localouteid);
      procid_t targetowner = localvid2owner[localstore.target(localouteid)];

      if (localvid_is_ghost(localtargetvid)) {
        block_synchronize_request2 &req = requests[targetowner].first;
        // need to synchronize outgoing vertex and outgoing edge
        // do outgoing vertex first
        if (std::binary_search(req.vid.begin(), requests[targetowner].second, 
                               local2globalvid[localtargetvid]) == false) {
          req.vid.push_back(local2globalvid[localtargetvid]);
          req.vidversion.push_back(localstore.vertex_version(localtargetvid));
          vertex_conditional_store vs;
          vs.hasdata = localstore.vertex_modified(localtargetvid);
          if (vs.hasdata) {
            localstore.set_vertex_modified(localtargetvid, false);
            vs.data.first = localstore.vertex_data(localtargetvid);
            vs.data.second = localstore.vertex_version(localtargetvid);
          }
          req.vstore.push_back(vs);
        }
        // now for the outgoing edge
      
        req.srcdest.push_back(std::make_pair(vid, local2globalvid[localtargetvid]));
        req.edgeversion.push_back(localstore.edge_version(localouteid));
        edge_conditional_store es;
        es.hasdata = localstore.edge_modified(localouteid);
        if (es.hasdata) {
          localstore.set_edge_modified(localouteid, false);
          es.data.first = localstore.edge_data(localouteid);
          es.data.second = localstore.edge_version(localouteid);
        }
        req.estore.push_back(es);
      }
    }
  }



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  async_synchronize_scope_callback(vertex_id_type vid, 
                                   boost::function<void (void)> callback){
    vertex_id_type localvid = global2localvid[vid];
    ASSERT_TRUE(scope_callbacks[localvid].callback == NULL);
    ASSERT_TRUE(boundary_scopes_set().find(vid) != boundary_scopes_set().end());
    // register the callback
    scope_callbacks[localvid].callback = callback;
    // construct the requests
    typedef std::map<procid_t, request_veciter_pair_type > map_type; 
    map_type requests;
    synchronize_scope_construct_req(vid, requests);
  
    //send the stuff
    typename map_type::iterator iter;
    iter = requests.begin();
  
    ASSERT_EQ(scope_callbacks[localvid].counter.value, 0);
    scope_callbacks[localvid].counter.inc(requests.size());
    while(iter != requests.end()) {

      // the reply target is 0. see reply_alot2
      rmi.remote_call(iter->first,
                      &distributed_graph<VertexData, EdgeData>::
                      async_get_alot2,
                      rmi.procid(),
                      iter->second.first,
                      0,
                      localvid);
      ++iter;
    }
  } // end of async synchronize scope callback




  /**
   * Waits for all asynchronous data synchronizations to complete
   */
  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  wait_for_all_async_syncs() {  pending_async_updates.wait(); }





  template <typename VertexData, typename EdgeData> 
  typename distributed_graph<VertexData, EdgeData>::vertex_conditional_store 
  distributed_graph<VertexData, EdgeData>::
  get_vertex_if_version_less_than(vertex_id_type vid,
                                  uint64_t  vertexversion,
                                  vertex_conditional_store &vdata) {
    vertex_conditional_store ret;
    size_t localvid = global2localvid[vid];
    uint64_t local_vertex_version = localstore.vertex_version(localvid);
    // Now I must the the owner of this vertex
    ASSERT_EQ(localvid2owner[localvid], rmi.procid());
#ifdef DGRAPH_DEBUG
    logstream(LOG_DEBUG) 
      << "get vertex: " << vid << ":" << vertexversion 
      << " vs " << local_vertex_version << ". " << vdata.hasdata 
      << std::endl;
#endif
    ret.hasdata = false;
  
    if (local_vertex_version  > vertexversion) {
      // send if I have a later version
      ret.hasdata = true;
      ret.data.first = localstore.vertex_data(localvid);
      ret.data.second = local_vertex_version;
    }
    else if (local_vertex_version == vertexversion) {
      // if version is the same and there is data, store and increment the version    
      if (vdata.hasdata) {
        localstore.increment_and_update_vertex(localvid, vdata.data.first);
      }
    }
    else {
      logstream(LOG_FATAL) << "Remote node attempted to update vertex " 
                           << vid << " with a newer version than the owner" << std::endl;
    }
    return ret;
  } // end of get_vertex_if_version_than




  template <typename VertexData, typename EdgeData> 
  typename distributed_graph<VertexData, EdgeData>::edge_conditional_store 
  distributed_graph<VertexData, EdgeData>::
  get_edge_if_version_less_than2(vertex_id_type source,
                                 vertex_id_type target,
                                 uint64_t  edgeversion,
                                 edge_conditional_store &edata) {
    edge_conditional_store ret;
    size_t localsource = global2localvid[source];
    size_t localtarget = global2localvid[target];


    std::pair<bool, edge_id_type> findret = 
      localstore.find(localsource, localtarget);
    assert(findret.first);
    edge_id_type localeid = findret.second;
  
    uint64_t  local_edge_version = localstore.edge_version(localeid);
    ret.hasdata = false;
#ifdef DGRAPH_DEBUG
    logstream(LOG_DEBUG) 
      << "get edge2: " << "(" << source << "->" << target << ")" 
      << ":" << edgeversion << " vs " << local_edge_version 
      << ". " << edata.hasdata << std::endl;
#endif
    if (local_edge_version > edgeversion) {
      ret.hasdata = true;
      ret.data.first = localstore.edge_data(localeid);
      ret.data.second = local_edge_version;
    }
    else if (local_edge_version == edgeversion) {
      if (edata.hasdata) {
        localstore.increment_and_update_edge(localeid, edata.data.first);
      }
    }
    else {
      logstream(LOG_FATAL) 
        << "Remote node attempted to update edge (" 
        << source <<  ", " << target 
        << ") with a newer version than the owner" << std::endl;
    }
    return ret;
  } // end of get_edge_if_version_less_than2



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  async_get_vertex_if_version_less_than(
                                        procid_t srcproc, 
                                        vertex_id_type vid, 
                                        uint64_t vertexversion,
                                        vertex_conditional_store &vdata) {
    rmi.remote_call(srcproc,
                    &distributed_graph<VertexData, EdgeData>::
                    reply_vertex_data_and_version,
                    vid,
                    get_vertex_if_version_less_than(vid, 
                                                    vertexversion, 
                                                    vdata));
  } // end of async_get_vertex_if_version_less_than


  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  async_get_edge_if_version_less_than2(
                                       procid_t srcproc, 
                                       vertex_id_type source, 
                                       vertex_id_type target, 
                                       uint64_t edgeversion,
                                       edge_conditional_store &edata) {
    rmi.remote_call(srcproc,
                    &distributed_graph<VertexData, EdgeData>::
                    reply_edge_data_and_version2,
                    source,
                    target,
                    get_edge_if_version_less_than2(source, 
                                                   target, 
                                                   edgeversion, 
                                                   edata));
  } // end of async_get_vertex_if_version_less_than


  template <typename VertexData, typename EdgeData> 
  typename distributed_graph<VertexData, EdgeData>::block_synchronize_request2& 
  distributed_graph<VertexData, EdgeData>::
  get_alot2(distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request) {
    std::vector<vertex_conditional_store> vresponse(request.vid.size());
    std::vector<edge_conditional_store> eresponse(request.srcdest.size());
    for (size_t i = 0;i < request.vid.size(); ++i) {
      request.vstore[i] = get_vertex_if_version_less_than(request.vid[i], 
                                                          request.vidversion[i], 
                                                          request.vstore[i]);
    }
    for (size_t i = 0;i < request.srcdest.size(); ++i) {
      request.estore[i] = get_edge_if_version_less_than2(request.srcdest[i].first, 
                                                         request.srcdest[i].second, 
                                                         request.edgeversion[i], 
                                                         request.estore[i]);
    }
    request.vidversion.clear();
    request.edgeversion.clear();
    return request;
  } // end of get_alot2



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  async_get_alot2(procid_t srcproc,
                  distributed_graph<VertexData, EdgeData>::
                  block_synchronize_request2 &request,
                  size_t replytarget,
                  size_t tag) {
    get_alot2(request);
    rmi.remote_call(srcproc,
                    &distributed_graph<VertexData, EdgeData>::reply_alot2,
                    request,
                    replytarget, 
                    tag);
  } // end of distributed_graph



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  reply_vertex_data_and_version(
                                vertex_id_type vid, 
                                distributed_graph<VertexData, EdgeData>::
                                vertex_conditional_store &vstore) {
    if (vstore.hasdata) update_vertex_data_and_version(vid, vstore);
    // the dc and procid are meaningless. Just pass something
    reply_increment_counter(rmi.dc(), 0, 
                            reinterpret_cast<size_t>(&pending_async_updates), dc_impl::blob());
  } // end of reply_vertex_data_and_version



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  reply_edge_data_and_version2(vertex_id_type source, 
                               vertex_id_type target, 
                               distributed_graph<VertexData, EdgeData>::edge_conditional_store &estore) {
    if (estore.hasdata) update_edge_data_and_version2(source, target, estore);
    reply_increment_counter(rmi.dc(), 0, 
                            reinterpret_cast<size_t>(&pending_async_updates), dc_impl::blob());

  } // end of reply_edge_data_and_version2



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  update_vertex_data_and_version(vertex_id_type vid, 
                                 distributed_graph<VertexData, EdgeData>::vertex_conditional_store &vstore) {
    vertex_id_type localvid = global2localvid[vid];
    // this must be a ghost
    ASSERT_NE(localvid2owner[localvid], rmi.procid());
#ifdef DGRAPH_DEBUG
    logstream(LOG_DEBUG) << "Receiving vertex " << vid << "(" << localvid << ")  "  << ". "
                         << vstore.data.second << " vs " << localstore.vertex_version(localvid) << std::endl;
#endif
    localstore.conditional_update_vertex(localvid, vstore.data.first, vstore.data.second);
  }

  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  update_edge_data_and_version2(
                                vertex_id_type source, 
                                vertex_id_type target, 
                                distributed_graph<VertexData, EdgeData>::edge_conditional_store &estore) {
    if (estore.hasdata) {
      vertex_id_type localsourcevid = global2localvid[source];
      vertex_id_type localtargetvid = global2localvid[target];
      ASSERT_NE(localvid2owner[localtargetvid], rmi.procid());
      std::pair<bool, edge_id_type> findret = localstore.find(localsourcevid, localtargetvid);
    
      assert(findret.first);
#ifdef DGRAPH_DEBUG
      logstream(LOG_DEBUG) << "Receiving edge (" << source << ","<<target << ")  " << ". "
                           << estore.data.second << " vs " << localstore.edge_version(findret.second) << std::endl;
#endif

      localstore.conditional_update_edge(findret.second, estore.data.first, estore.data.second);
    }
  }


  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::update_alot2(
                                                             distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request) {
    for (size_t i = 0;i < request.vid.size(); ++i) {
      update_vertex_data_and_version(request.vid[i], request.vstore[i]);
    }

    for (size_t i = 0;i < request.srcdest.size(); ++i) {
      update_edge_data_and_version2(request.srcdest[i].first, request.srcdest[i].second, request.estore[i]);
    }
  }


  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::reply_alot2(
                                                            distributed_graph<VertexData, EdgeData>::block_synchronize_request2 &request,
                                                            size_t replytarget,
                                                            size_t tag) {
    update_alot2(request);
  
    // special handling for callbacks
    if (replytarget != 0)  {
      reply_increment_counter(rmi.dc(), 0, 
                              replytarget, dc_impl::blob());  
    }
    else {
      // tag is local vid
      vertex_id_type localvid = tag;

      ASSERT_TRUE(scope_callbacks[localvid].callback != NULL);
      if (scope_callbacks[localvid].counter.dec() == 0) {
        // make a copy of it and clear the callbacks entry.
        boost::function<void (void)> tmp = scope_callbacks[localvid].callback;
        scope_callbacks[localvid].callback = NULL;
        tmp();
      }
    } 
  }


  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::synchronize_all_vertices(bool async) {
    foreach(vertex_id_type vid, ghostvertices) {
      synchronize_vertex(vid, async);
    }
  }

  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::synchronize_all_edges(bool async) {
    foreach(vertex_id_type vid, ghostvertices) {
      foreach(edge_id_type eid, localstore.in_edge_ids(global2localvid[vid])) {
        synchronize_edge(eid, async);
      }
    }
  }

  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::synchronize_all_scopes(bool async) {
    foreach(vertex_id_type vid, boundaryscopes) {
      synchronize_scope(vid, async);
    }
  }



  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  update_vertex_data_and_version_and_reply(vertex_id_type vid, 
                                           distributed_graph<VertexData, EdgeData>::vertex_conditional_store &vstore,
                                           procid_t srcproc,
                                           size_t reply) {
                        
    update_vertex_data_and_version(vid, vstore);

    if (srcproc != procid_t(-1)) {
      rmi.dc().remote_call(srcproc, reply_increment_counter, reply, dc_impl::blob());
    }
  }




  template <typename VertexData, typename EdgeData> 
  void distributed_graph<VertexData, EdgeData>::
  update_edge_data_and_version_and_reply2(vertex_id_type source, 
                                          vertex_id_type target, 
                                          distributed_graph<VertexData, EdgeData>::edge_conditional_store &estore,
                                          procid_t srcproc, size_t reply) {
    update_edge_data_and_version2(source, target, estore);
    if (srcproc != procid_t(-1)) {
      rmi.dc().remote_call(srcproc, reply_increment_counter, reply, dc_impl::blob());
    }
  }




  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  push_owned_vertex_to_replicas(vertex_id_type vid, bool async, bool untracked) {
    vertex_id_type localvid = global2localvid[vid];
    // get the replicas
    const std::vector<procid_t>& replicas = localvid_to_replicas(localvid);
    // owner is a replica too. if there are no other replicas quit
    if (replicas.size() <= 1) return;
  
    dc_impl::reply_ret_type ret(true, replicas.size() - 1);
  
    // if async, set the return reply to go to the global pending push updates
    size_t retptr;
    if (async == false) {
      retptr = reinterpret_cast<size_t>(&ret);
    }
    else {
      retptr = reinterpret_cast<size_t>(&pending_push_updates);
      pending_push_updates.flag.inc(replicas.size() - 1);
    }
  
    /**
       If untracked, set the reply procid to -1. That will mean that
       no reply will be returned at all
    */
    procid_t srcprocid;
    if (untracked == false) srcprocid = rmi.procid();
    else srcprocid = procid_t(-1);
    // build the store
    vertex_conditional_store vstore;
    vstore.hasdata = true;
    vstore.data.first = localstore.vertex_data(localvid);
    vstore.data.second = localstore.vertex_version(localvid);
  
    foreach(procid_t proc, replicas) {
      if (proc != rmi.procid()) {
#ifdef DGRAPH_DEBUG
        logger(LOG_DEBUG, "Pushing vertex %d to proc %d", vid, proc);
#endif
        rmi.remote_call(proc,
                        &distributed_graph<VertexData, EdgeData>::update_vertex_data_and_version_and_reply,
                        vid,
                        vstore,
                        srcprocid,
                        retptr);
      }
    }
    if (async == false && untracked == false)  ret.wait();
  }


  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::push_owned_edge_to_replicas(edge_id_type eid, bool async, bool untracked) {
    // get the replicas
    vertex_id_type globalsource = source(eid);
    vertex_id_type globaltarget = target(eid);
    // firstly I must own this edge
    if (!is_owned(globaltarget)) return;
  
    // Now, there are 2 cases. 
    // Case 1, I own both source and target.
    //         in that case there is nothing to sync. Any other
    //         machine at most has ghosts of source and target, but will
    //         not have the end itself
    // Case 2, I own the target, but not the source
    //         then the only replica is the owner of the source

    procid_t sendto;
  
    if (is_owned(globalsource)) {
      return;
    }
    else {
      // otherwise, it is the owner of the source
      sendto = localvid2owner[localstore.source(eid)];
    }

    dc_impl::reply_ret_type ret(true, 1);  
    // if async, set the return reply to go to the global pending push updates
    size_t retptr;
    if (async == false) {
      retptr = reinterpret_cast<size_t>(&ret);
    }
    else {
      retptr = reinterpret_cast<size_t>(&pending_push_updates);
    }
  
    /**
       If untracked, set the reply procid to -1. That will mean that
       no reply will be returned at all
    */
    procid_t srcprocid;
    if (untracked == false) srcprocid = rmi.procid();
    else srcprocid = procid_t(-1);
    // build the store
    edge_conditional_store estore;
    estore.hasdata = true;
    estore.data.first = localstore.edge_data(eid);
    estore.data.second = localstore.edge_version(eid);
  
#ifdef DGRAPH_DEBUG
    logstream(LOG_DEBUG) << "Pushing edge (" << globalsource << ", " << globaltarget << ") to proc " << sendto << std::endl;
#endif
    rmi.remote_call(sendto,
                    &distributed_graph<VertexData, EdgeData>::update_edge_data_and_version_and_reply2,
                    globalsource,
                    globaltarget,
                    estore,
                    srcprocid,
                    retptr);

    if (async == false && untracked == false)  ret.wait();
  }


  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::push_all_owned_vertices_to_replicas() {
    // perform block collective
    std::vector<block_synchronize_request2> blockpushes(rmi.numprocs()); 

    foreach(vertex_id_type vid, boundary_scopes()) {

      vertex_id_type localvid = global2localvid[vid];
      // get the replicas
      const std::vector<procid_t>& replicas = localvid_to_replicas(localvid);
      // owner is a replica too. if there are no other replicas quit
      if (replicas.size() <= 1) continue;
  
 
      // build the store
      vertex_conditional_store vstore;
      vstore.hasdata = true;
      vstore.data.first = localstore.vertex_data(localvid);
      vstore.data.second = localstore.vertex_version(localvid);
  
      foreach(procid_t proc, replicas) {
        if (proc != rmi.procid()) {
          blockpushes[proc].vid.push_back(vid);
          blockpushes[proc].vidversion.push_back(localstore.vertex_version(localvid));
          blockpushes[proc].vstore.push_back(vstore);
          if (blockpushes[proc].vid.size() >= 1024*1024/sizeof(VertexData)) {
            rmi.remote_call(proc,
                            &distributed_graph<VertexData, EdgeData>::update_alot2,
                            blockpushes[proc]);
            blockpushes[proc].clear();
          }
        }
      }
    }

    for(size_t proc = 0; proc < rmi.numprocs(); ++proc) {
      if (blockpushes[proc].vid.size() > 0) {
        assert(proc != rmi.procid());
        rmi.remote_call(proc,
                        &distributed_graph<VertexData, EdgeData>::update_alot2,
                        blockpushes[proc]);
        blockpushes[proc].clear();
      }
    }
  }

  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::push_all_owned_edges_to_replicas() {
    std::vector<std::vector<block_synchronize_request2> > blockpushes(omp_get_max_threads()); 
    for (size_t i = 0;i < blockpushes.size(); ++i) blockpushes[i].resize(rmi.numprocs());
  
  
#pragma omp parallel for
    for (long i = 0;i < (long)ghostvertices.size(); ++i) {
      int thrid = omp_get_thread_num();
    
      vertex_id_type vid = ghost_vertices()[i];
      vertex_id_type localvid = global2localvid[vid];
      procid_t proc = localvid2owner[localvid];
      foreach(edge_id_type localeid, localstore.out_edge_ids(localvid)) {
        vertex_id_type targetvid = local2globalvid[localstore.target(localeid)];
        edge_conditional_store estore;
        estore.hasdata = true;
        estore.data.first = localstore.edge_data(localeid);
        estore.data.second = localstore.edge_version(localeid);
  

        blockpushes[thrid][proc].srcdest.push_back(std::make_pair<vertex_id_type, vertex_id_type>(vid, targetvid));
        blockpushes[thrid][proc].edgeversion.push_back(localstore.edge_version(localeid));
        blockpushes[thrid][proc].estore.push_back(estore);
        if (blockpushes[thrid][proc].srcdest.size() >= 1*1024*1024/sizeof(EdgeData)) {
          rmi.remote_call(proc,
                          &distributed_graph<VertexData, EdgeData>::update_alot2,
                          blockpushes[thrid][proc]);
          blockpushes[thrid][proc].clear();
        }
      }
    }
    for (size_t i = 0;i < blockpushes.size(); ++i) {
      for(size_t proc = 0; proc < rmi.numprocs(); ++proc) {
        if (blockpushes[i][proc].srcdest.size() > 0) {
          assert(proc != rmi.procid());
          rmi.remote_call(proc,
                          &distributed_graph<VertexData, EdgeData>::update_alot2,
                          blockpushes[i][proc]);
          blockpushes[i][proc].clear();
        }
      }
    }
  }


  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::push_owned_scope_to_replicas(vertex_id_type vid, 
                                                                             bool onlymodified, 
                                                                             bool clearmodified, 
                                                                             bool async,
                                                                             bool untracked) {
    // fast exit if this is not on a boundary
    if (boundaryscopesset.find(vid) == boundaryscopesset.end()) return;
    if (0) {
      if (is_owned(vid)) {
        vertex_id_type localvid = global2localvid[vid];
        if (localstore.vertex_modified(localvid)) {
          localstore.set_vertex_modified(localvid, false);
          push_owned_vertex_to_replicas(vid, async, untracked);
        
        }
        foreach(edge_id_type eid, in_edge_ids(vid)) {
          if (localstore.edge_modified(eid)) {
            localstore.set_edge_modified(eid, false);
            push_owned_edge_to_replicas(eid, async, untracked);
          }
        }
      }
    }
    else {
      if (is_owned(vid)) {
        push_owned_vertex_to_replicas(vid, async, untracked);
        foreach(edge_id_type eid, in_edge_ids(vid)) {
          push_owned_edge_to_replicas(eid, async, untracked);
        }
      }
    }
  }

  template <typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::wait_for_all_async_pushes() {
    pending_push_updates.wait();
  }








































  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const distributed_graph<VertexData, EdgeData>& graph) {
    graph.print(out);
    return out;
  }


}

#include <graphlab/macros_undef.hpp>
#endif

