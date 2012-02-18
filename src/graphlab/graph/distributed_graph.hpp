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

#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_HPP

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
#include <functional>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>


#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/mpi_tools.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/graph/graph2.hpp>
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_batch_ingress.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab { 


  // CLASS GRAPH ==============================================================>  
  template<typename VertexData, typename EdgeData>
  class distributed_graph {
  public:

    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;
    /// Type for vertex colors 
    typedef graphlab::vertex_color_type vertex_color_type;


    /// The type of the local graph used to store the graph data 
    typedef graphlab::graph2<VertexData, EdgeData> local_graph_type;

    typedef idistributed_ingress<VertexData, EdgeData> 
    idistributed_ingress_type;


    /** 
     * The type of the local vertex id and local edge id.
     * While this is the same as the
     * vertex_id_type giving it a separate name will make method calls
     * clearer.
     */
    typedef typename local_graph_type::vertex_id_type  lvid_type;
    typedef typename local_graph_type::edge_id_type    leid_type;
    typedef typename local_graph_type::edge_type       local_edge_type;
    typedef typename local_edge_type::edge_dir         ledge_dir_type;

    
    /** This class represents an edge with source() and target()*/
    class edge_type {
    private:
      vertex_id_type _src;
      vertex_id_type _tar;
      edge_id_type _eid;
      procid_t _owner;
      bool _empty;      
      inline edge_id_type edge_id() const { return _eid; }
      friend class distributed_graph;
    public:
      edge_type (vertex_id_type src, vertex_id_type tar,
                 edge_id_type eid, procid_t owner):
        _src(src), _tar(tar), _eid(eid), _owner(owner), _empty(false) { }
      edge_type () : _src(-1), _tar(-1), _eid(-1), _owner(-1), _empty(true) { }
      inline vertex_id_type source() const { ASSERT_FALSE(empty()); return _src; }
      inline vertex_id_type target() const { ASSERT_FALSE(empty()); return _tar; }

      procid_t owner() const { return _owner; }
      bool empty() const { return _empty; }
    }; // end of class edge_type.



    /** This class represents an edge list stored on local machine*/
    class local_edge_list_type {
      typedef typename local_graph_type::edge_list_type lgraph_edge_list_type;
      typedef typename lgraph_edge_list_type::const_iterator 
      local_const_iterator_type;
      
      struct edge_functor : 
        public std::unary_function<local_edge_type, edge_type> {

        edge_functor(const distributed_graph* graph_ptr = NULL) : 
          graph_ptr(graph_ptr) { }

        edge_type operator()(const local_edge_type& edge) const {
          if (edge.empty()) {
            return edge_type();
          } else {
            ASSERT_TRUE(graph_ptr != NULL);
            const vertex_id_type source = graph_ptr->global_vid(edge.source());
            const vertex_id_type target = graph_ptr->global_vid(edge.target());
            const vertex_id_type eid = 
              graph_ptr->global_eid(graph_ptr->local_eid(edge));
            const procid_t owner = graph_ptr->rpc.procid();
            return edge_type(source, target, eid, owner);
          }
        } // end of operator()
        const distributed_graph* graph_ptr;
      }; // end of edge_functor

    public:
      typedef boost::
      transform_iterator<edge_functor, local_const_iterator_type> iterator;
      typedef iterator const_iterator;
      typedef edge_type value_type;
    private:
      const_iterator begin_iter, end_iter;      
    public:
      local_edge_list_type(const distributed_graph* graph_ptr = NULL,
                           const lgraph_edge_list_type& lgraph_edge_list =
                           lgraph_edge_list_type()) :
        begin_iter(lgraph_edge_list.begin(), edge_functor(graph_ptr)),
        end_iter(lgraph_edge_list.end(), edge_functor(graph_ptr)) { }
      size_t size() const { return end_iter - begin_iter; }
      edge_type operator[](size_t i) const {
        ASSERT_LT(i, size()); return *(begin_iter + i);
      }
      const_iterator begin() const { return begin_iter; }
      const_iterator end() const { return end_iter; }
      bool empty() const { return size() == 0; }  
    }; // end of local_edge_list_type


    /** This class represents an iterator over edge list stored on remote machine*/
    class remote_edge_iterator {
    public:
      typedef typename local_edge_type::edge_dir iterator_type;
      typedef local_edge_type reference;
    public:
      // Cosntructors
      remote_edge_iterator(procid_t proc) : 
        offset(-1), empty(true), proc(proc) { }

      remote_edge_iterator(procid_t proc, vertex_id_type _center, 
                            size_t _offset, iterator_type _itype) : 
        proc(proc), center(_center), offset(_offset), 
        itype(_itype), empty(false) { }

      remote_edge_iterator(const remote_edge_iterator& it) :
        proc(it.proc), center(it.center), offset(it.offset), 
        itype(it.itype), empty(it.empty) { }

      inline edge_type operator*() const  {
        ASSERT_TRUE(!empty);
        return make_value();
      }

      typedef boost::detail::
      operator_arrow_result<edge_type, edge_type, edge_type*> arrow_type;
      inline typename arrow_type::type operator->() const {
        return arrow_type::make(make_value());
      }


      inline bool operator==(const remote_edge_iterator& it) const {
        return (proc == it.proc && empty && it.empty) || 
          (proc == it.proc && empty == it.empty &&
           itype == it.itype && center == it.center && 
           offset == it.offset);
      }

      inline bool operator!=(const remote_edge_iterator& it) const { 
        return !(*this == it);
      }

      inline remote_edge_iterator& operator++() {
        ASSERT_TRUE(!empty);
        ++offset;
        return *this;
      }

      inline remote_edge_iterator operator++(int) {
        ASSERT_TRUE(!empty);
        const remote_edge_iterator copy(*this);
        operator++();
        return copy;
      }

      inline int operator-(const remote_edge_iterator& it) const {
        ASSERT_TRUE(!empty && itype == it.itype && center == it.center);
        return offset - it.offset;
      }

      inline remote_edge_iterator operator+(size_t i) const {
        remote_edge_iterator retval(proc, center, offset+i, itype);
        return retval;
      }

      // Generate the ret value of the iterator.
      inline edge_type make_value() const {
        // Not implemented;
        logstream(LOG_WARNING) << "make_value not implemented. " << std::endl;
        ASSERT_TRUE(false);
        return edge_type();
      }

    private:
      procid_t proc;
      vertex_id_type center;
      size_t offset;
      iterator_type itype;
      bool empty;
    }; // end of remote_edge_iterator



    class edge_iterator : 
      public std::iterator<std::forward_iterator_tag, edge_type> {
      typedef typename local_edge_list_type::iterator l_iterator;
      typedef remote_edge_iterator r_iterator;
    public:
      // Constructors
      edge_iterator () : total_counts(0), proc(-1), 
                         offset(-1), global_offset(-1) { }

      edge_iterator(l_iterator l_iter, size_t local_counts, 
                    std::vector<r_iterator> r_iters, 
                    std::vector<size_t> remote_counts, 
                    size_t total_counts, size_t global_offset) :
        l_iter(l_iter), local_counts(local_counts), r_iters(r_iters), 
        remote_counts(remote_counts), total_counts(total_counts), 
        proc(-1), offset(0),  global_offset(global_offset) { }

      edge_iterator (const edge_iterator& it) :
        l_iter(it.l_iter), local_counts(local_counts), 
        r_iters(it.r_iters), remote_counts(it.remote_counts), 
        total_counts(it.total_counts), proc(it.proc), 
        offset(it.offset), global_offset(it.global_offset) { }

      edge_type operator*() const {
        ASSERT_TRUE(global_offset < total_counts);
        return proc > remote_counts.size() ? 
          *l_iter : *(r_iters[proc]);
      }
      
      bool operator==(const edge_iterator& it) const {
        return (l_iter==it.l_iter && global_offset == it.global_offset);
        /* Because the local iterator has all information about this
           edge list the rest fields should be equal assuming correct
           construction of this iterator. For example, (proc ==
           it.proc) && (offset == it.offset) && (total_counts ==
           it.total_counts) && (local_counts == it.local_counts);
        */
      }

      bool operator!=(const edge_iterator& it) const { 
        return (l_iter != it.l_iter || global_offset != it.global_offset);
      }

      edge_iterator& operator++() {
        if (proc > remote_counts.size()) { // Still in local iterator.
          if (offset < local_counts - 1) {
            ++offset;
          } else { // Find new proc and start remote iterator.
            proc = 0; offset= 0;
            while (proc < remote_counts.size() && remote_counts[proc] == 0) 
              ++proc;
          }
        } else { // Still in same proc iterator.
          if (offset < remote_counts[proc] -1) {
            ++offset;
          } else { // Find new proc iterator.
            offset = 0;
            while (proc < remote_counts.size() && remote_counts[proc] == 0) 
              ++proc;
          }
        }
        ++global_offset;
        ASSERT_TRUE(global_offset <= total_counts);
      } // end of operator++

      edge_iterator operator++(int) {
        const edge_iterator copy(*this);
        operator++();
        return copy;
      }

      int operator-(const edge_iterator& it) const {
        ASSERT_TRUE(l_iter == it.l_iter);
        return (global_offset - it.global_offset);
      }

      edge_iterator operator+(size_t i) const {
        edge_iterator copy(*this);
        copy.global_offset += i;
        ASSERT_TRUE(copy.global_offset < copy.total_counts);
        size_t t = copy.global_offset;
        if (t >= local_counts) {
          t = t - local_counts; copy.proc = 0;
          while (copy.proc < remote_counts.size() && 
                 t >= remote_counts[copy.proc]) {
            t -= remote_counts[copy.proc];
            ++copy.proc;
          }
        }
        copy.offset = t;
        return copy;
      }

    private:
      l_iterator l_iter;
      size_t local_counts;
      std::vector<r_iterator> r_iters;
      std::vector<size_t> remote_counts;
      size_t total_counts;
      size_t proc;
      size_t offset;
      size_t global_offset;
    }; // end of edge_iterator




    /** This class represents a general edge list */
    /*  Lazy list of a hybrid collection of local and remote edges*/
    class edge_list_type {
      // Type interface for boost foreach.
    public:
      typedef edge_iterator iterator;
      typedef edge_iterator const_iterator;
      typedef edge_type value_type;

    public:
      // Construct an empty edge list
      edge_list_type() : list_size(0) { }
      // Cosntruct an edge_list with begin and end. 
      edge_list_type(edge_iterator begin, edge_iterator end) : 
        begin_ptr(begin), end_ptr(end) { 
        list_size = (size_t)(end_ptr-begin_ptr);
      }
      // Copy constructor
      edge_list_type(const edge_list_type& other) : 
        begin_ptr(other.begin_ptr), end_ptr(other.end_ptr), 
        list_size(other.list_size) { }

      inline size_t size() const { return list_size;}
      inline edge_type operator[](size_t i) const {
        ASSERT_LT(i, list_size);
        return *(begin_ptr + i);
      }
      iterator begin() const { return begin_ptr; }
      iterator end() const { return end_ptr; }
      bool empty() const { return size() == 0; }
    private:
      edge_iterator begin_ptr;
      edge_iterator end_ptr;
      size_t list_size;
    }; // end of class edge_list.


    /**
     * The vertex record stores information associated with each
     * vertex on this proc
     */
    struct vertex_record {
      /// The official owning processor for this vertex
      procid_t owner; 
      /// The local vid of this vertex on this proc
      vertex_id_type gvid;
      /// The number of in edges
      size_t num_in_edges;
      /// The nubmer of out edges
      size_t num_out_edges;
      /// The set of proc that mirror this vertex.
      std::vector<procid_t> mirrors;
      vertex_record() : 
        owner(-1), gvid(-1), num_in_edges(0), num_out_edges(0) { }
      vertex_record(const vertex_id_type& vid) : 
        owner(-1), gvid(vid), num_in_edges(0), num_out_edges(0) { }
      procid_t get_owner () const {
        return owner;
      }
      const std::vector<procid_t>& get_replicas () const {
        return mirrors;
      }
    }; // end of vertex_record


    /// The master vertex record map
    // typedef boost::unordered_map<vertex_id_type, vertex_record>  vid2record_type;
    typedef std::vector<vertex_record> lvid2record_type;
      
    // PRIVATE DATA MEMBERS ===================================================> 
    /** The rpc interface for this class */
    mutable dc_dist_object<distributed_graph> rpc;

    /** The local graph data */
    local_graph_type local_graph;
    mutex local_graph_lock;
    
    /** The map from global vertex ids to vertex records */
    lvid2record_type lvid2record;
    mutex lvid2record_lock;
    
    /** The map from local vertex ids back to global vertex ids */
    boost::unordered_map<vertex_id_type, lvid_type> vid2lvid;
        
    /** The global number of vertices and edges */
    size_t nverts, nedges;

    /** The number of vertices owned by this proc */
    size_t local_own_nverts;

    /** The global number of vertex replica */
    size_t nreplica;

    /** The beginning edge id for this machine */
    size_t begin_eid;

    /** pointer to the distributed ingress object*/
    idistributed_ingress_type* ingress_ptr; 

  public:

    // CONSTRUCTORS ==========================================================>
    distributed_graph(distributed_control& dc) : 
      rpc(dc, this), nverts(0), nedges(0), local_own_nverts(0), nreplica(0) {
      rpc.barrier();
      typedef distributed_batch_ingress<vertex_data_type, edge_data_type>
        distributed_batch_ingress_type;
      ingress_ptr = new distributed_batch_ingress_type(dc, *this);
    }



    // METHODS ===============================================================>

    /**
     * Finalize is used to complete graph ingress by resolving vertex
     * ownship and completing local data structures.
     */
    void finalize() {
      ASSERT_NE(ingress_ptr, NULL);
      ingress_ptr->finalize();
      rpc.barrier(); delete ingress_ptr; ingress_ptr = NULL;
    }
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return nverts; }

    /** \brief Get the number of edges */
    size_t num_edges() const { return nedges; }

    /** \brief Get the size of replica */
    size_t num_replica() const { return nreplica; }

    /** Return the color of a vertex */
    vertex_color_type color(vertex_id_type vid) const { 
      //TODO: IMPLEMENT
      logstream(LOG_FATAL) << "Color not implemented" << std::endl; 
      return -1;
    }

    //! Get the rerverse edge 
    edge_type reverse_edge(const edge_type& edge) const {      
      //TODO: IMPLEMENT
      logstream(LOG_FATAL) << "Reverse edge not implemented" << std::endl; 
      return edge;
    }

    //! find a particular edge
    edge_type find(vertex_id_type source,
                   vertex_id_type target) const {
      //TODO: IMPLEMENT
      logstream(LOG_FATAL) << "find edge not implemented" << std::endl; 
      return edge_type();
    }



    /** \brief Get the number of vertices local to this proc */
    size_t num_local_vertices() const { return local_graph.num_vertices(); }

    /** \brief Get the number of edges local to this proc */
    size_t num_local_edges() const { return local_graph.num_edges(); }

    /** \brief Get the number of vertices owned by this proc */
    size_t num_local_own_vertices() const { return local_own_nverts; }

    /** \brief get the local vertex id */
    lvid_type local_vid (const vertex_id_type vid) const {
      typename boost::unordered_map<vertex_id_type, lvid_type>::
        const_iterator iter = vid2lvid.find(vid);
      ASSERT_TRUE(iter != vid2lvid.end());
      return iter->second;
    } // end of local_vertex_id

    vertex_id_type global_vid (const lvid_type lvid) const { 
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid].gvid;
    } // end of global_vertex_id

    vertex_record& get_vertex_record(const vertex_id_type vid) {
      ASSERT_LT(local_vid(vid), lvid2record.size());
      return lvid2record[local_vid(vid)];
    }

    const vertex_record& get_vertex_record(const vertex_id_type vid) const {
      ASSERT_LT(local_vid(vid), lvid2record.size());
      return lvid2record[local_vid(vid)];
    }

    vertex_record& l_get_vertex_record(const lvid_type lvid) {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid];
    }

    const vertex_record& l_get_vertex_record(const vertex_id_type lvid) const {
      ASSERT_LT(lvid, lvid2record.size());
      return lvid2record[lvid];
    }

    local_graph_type& get_local_graph() {
      return local_graph;
    }

    const local_graph_type& get_local_graph() const {
      return local_graph;
    }

    edge_id_type global_eid(const leid_type eid) const {
      return (begin_eid + eid);
    } 

    leid_type local_eid(const edge_id_type eid) const {
      return (eid - begin_eid);
    } 

    leid_type local_eid(const local_edge_type& eid) const {
      return local_graph.edge_id(eid);
    }

    //! Get all the edge which edge.target() == v
    edge_list_type in_edges(const vertex_id_type vid) const {
      // Not implemented.
      logstream(LOG_WARNING) << "in_edges not implemented. " << std::endl;
      ASSERT_TRUE(false);
      return edge_list_type();
    }


    //! Get the number of edges which edge.target() == v
    size_t num_in_edges(const vertex_id_type vid) const {
      return get_vertex_record(vid).num_in_edges;
    }

    //! Get all the edges which edge.source() == v
    edge_list_type out_edges(const vertex_id_type vid) const {
      // Not implemented;
      logstream(LOG_WARNING) << "out_edges not implemented. " << std::endl;
      ASSERT_TRUE(false);
      return edge_list_type();
    }

    //! Get the number of edges which edge.source() == v
    size_t num_out_edges(const vertex_id_type vid) const {
      return get_vertex_record(vid).num_out_edges;
    }

    // Get all the local edge which edge.target() == v
    local_edge_list_type l_in_edges(const vertex_id_type vid) const {
      return local_edge_list_type(this, local_graph.in_edges(local_vid(vid)));
    }

    // Get the number of local edges which edge.target() == v
    size_t l_num_in_edges(const vertex_id_type vid) const { 
      return local_graph.num_in_edges(local_vid(vid));
    }

    // Get all the local edges which edge.source() == v
    local_edge_list_type l_out_edges(const vertex_id_type vid) const {
      return local_edge_list_type(this, local_graph.out_edges(local_vid(vid)));
    }

    // Get the number of local edges which edge.source() == v
    size_t l_num_out_edges(const vertex_id_type vid) const {
      return local_graph.num_out_edges(local_vid(vid));
    }

    /** \brief Returns a reference to the data stored on the vertex
        v. */
    VertexData& vertex_data(vertex_id_type vid) {
      return local_graph.vertex_data(local_vid(vid));
    }
    
    /** \brief Returns a constant reference to the data stored on the
        vertex v */
    const VertexData& vertex_data(vertex_id_type vid) const {
      return local_graph.vertex_data(local_vid(vid));
    }

    /** \brief Returns a reference to the data stored on the edge
        source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target){
      return local_graph.edge_data(local_vid(source), local_vid(target));
    }
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, 
                              vertex_id_type target) const {
      return local_graph.edge_data(local_vid(source), local_vid(target));
    }

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(const edge_type& edge) {
      local_edge_type l_edge(local_vid(edge.source()), local_vid(edge.target()),
                             local_eid(edge.edge_id()), local_edge_type::OUTEDGE);
      return local_graph.edge_data(l_edge);
    }
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(const edge_type& edge) const {
      local_edge_type l_edge(local_vid(edge.source()), local_vid(edge.target()),
                             local_eid(edge.edge_id()), local_edge_type::OUTEDGE);
      return local_graph.edge_data(l_edge);
    }

   
    /** 
     * \brief Creates a vertex containing the vertex data
     */
    void add_vertex(const vertex_id_type& vid, 
                    const VertexData& vdata = VertexData() ) {
      ASSERT_NE(ingress_ptr, NULL);
      ingress_ptr->add_vertex(vid, vdata);
    }

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData()) {
      ASSERT_NE(ingress_ptr, NULL);
      ingress_ptr->add_edge(source, target, edata);
    }



    void resize (size_t n) { }


  }; // End of graph

 
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

