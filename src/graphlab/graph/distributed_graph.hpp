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


#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>


#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/mpi_tools.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/graph/graph.hpp>




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
    typedef graphlab::graph<VertexData, EdgeData> local_graph_type;
    /** 
     * The type of the local vertex id.  While this is the same as the
     * vertex_id_type giving it a separate name will make method calls
     * clearer.
     */
    typedef typename local_graph_type::vertex_id_type lvid_type;

    
    /** This class represents an edge with source() and target()*/
    class edge_type {
      typedef typename local_graph_type::edge_type local_edge_type;
      local_edge_type local_edge;
      const distributed_graph* graph_ptr;
      edge_type(const distributed_graph* graph_ptr, 
                const local_edge_type local_edge) : 
        graph_ptr(graph_ptr), local_edge(local_edge) { }     
    public:
      edge_type() : graph_ptr(NULL) { }
      inline vertex_id_type source() const {
        ASSERT_FALSE(empty()); 
        // TODO: implement local to global vid conversion
      }
      inline vertex_id_type target() const { 
        ASSERT_FALSE(empty());
        // TODO: implement local to global vid conversion
      }
      bool empty() const { return graph_ptr == NULL; }
      friend class distributed_graph;
    }; // end of class edge_type.
    friend class edge_type;


    /** This class represents a lazy list of edge_type. */
    class edge_list_type {
    public:
      typedef typename edge_type::local_edge_type local_edge_type;
      typedef typename local_graph_type::edge_list_type local_edge_list_type;
      typedef typename local_edge_list_type::const_iterator 
      local_const_iterator_type;
      struct edge_functor : 
        public std::unary_function<local_edge_type, edge_type> {
        const distributed_graph* graph_ptr;
        edge_functor(const distributed_graph* graph_ptr) : 
          graph_ptr(graph_ptr) { }
        edge_type operator()(const local_edge_type& edge) const {
          return edge_type(graph_ptr, edge);
        }
      }; // end of edge_functor
      typedef boost::
      transform_iterator<edge_functor, local_const_iterator_type> iterator;

      typedef iterator const_iterator;
      typedef edge_type value_type;
    private:
      const_iterator begin_iter, end_iter;
    public:
      edge_list_type(const distributed_graph* graph_ptr = NULL,
                     const local_edge_list_type& local_edge_list =
                     local_edge_list_type()) :
        begin_iter(local_edge_list.begin(), edge_functor(graph_ptr)),
        end_iter(local_edge_list.end(), edge_functor(graph_ptr)) { }
      size_t size() const { return end_iter - begin_iter; }
      edge_type operator[](size_t i) const {
        ASSERT_LT(i, size()); return *(begin_iter + i);
      }
      const_iterator begin() const { return begin_iter; }
      const_iterator end() const { return end_iter; }
      bool empty() const { return size() == 0; }  
    }; // end of class edge_list_type


  private:
    /**
     * The vertex record stores information associated with each
     * vertex on this proc
     */
    struct vertex_record {
      /// The official owning processor for this vertex
      procid_t owner; 
      /// The local vid of this vertex on this proc
      lvid_type lvid;
      /// The number of in edges
      size_t num_in_edges;
      /// The nubmer of out edges
      size_t num_out_edges;
      /// The set of proc that mirror this vertex.
      std::vector<procid_t> mirrors;
      vertex_record() : 
        owner(-1), lvid(-1), 
        num_in_edges(0), num_out_edges(0) { }
    }; // vertex_record

    /// The master vertex record map
    typedef boost::unordered_map<vertex_id_type, vertex_record> 
    vid2record_type;
      
    // PRIVATE DATA MEMBERS ===================================================>        
    /** The rpc interface for this class */
    mutable dc_dist_object<distributed_graph> rpc;

    /** The local graph data */
    local_graph_type local_graph;
    mutex local_graph_lock;
    
    /** The map from global vertex ids to vertex records */
    vid2record_type vid2record;
    mutex vid2record_lock;
    
    /** The map from local vertex ids back to global vertex ids */
    std::vector<vertex_id_type> lvid2vid;

    /** The global number of vertices and edges */
    size_t nverts, nedges;
   
  public:

    // CONSTRUCTORS ==========================================================>
    distributed_graph(distributed_control& dc) : rpc(dc, this) {
      rpc.barrier();
    }



    // METHODS ===============================================================>

    /**
     * Finalize is used to complete graph ingress by resolving vertex
     * ownship and completing local data structures.
     */
    void finalize();
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return nverts; }

    /** \brief Get the number of edges */
    size_t num_edges() const { return nedges; }

    /** \brief Get the number of vertices local to this proc */
    size_t num_local_vertices() const { return local_graph.num_vertices(); }

    /** \brief Get the number of edges local to this proc */
    size_t num_local_edges() const { return local_graph.num_edges(); }

    /** \brief get the local vertex id */
    lvid_type local_to_global(const vertex_id_type vid) const { 
      return vrecord(vid).lvid;
    } // end of local_vertex_id

    vertex_id_type global_to_local(const lvid_type lvid) const { 
      ASSERT_LT(lvid, lvid2vid.size());
      return lvid2vid[lvid];
    } // end of global_vertex_id


    //! Get all the edge which edge.target() == v
    edge_list_type in_edges(const vertex_id_type vid) const;

    //! Get the number of edges which edge.target() == v
    size_t num_in_edges(const vertex_id_type v) const;

    //! Get all the edges which edge.source() == v
    edge_list_type out_edges(const vertex_id_type v) const;

    //! Get the number of edges which edge.source() == v
    size_t num_out_edges(const vertex_id_type v) const;

    /** \brief Returns a reference to the data stored on the vertex
        v. */
    VertexData& vertex_data(vertex_id_type v);
    
    /** \brief Returns a constant reference to the data stored on the
        vertex v */
    const VertexData& vertex_data(vertex_id_type v) const;

    /** \brief Returns a reference to the data stored on the edge
        source->target. */
    EdgeData& edge_data(vertex_id_type source, vertex_id_type target);
    
    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, 
                              vertex_id_type target) const;

    /** \brief Returns a reference to the data stored on the edge e */
    EdgeData& edge_data(const edge_type& edge);
    
    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(const edge_type& edge) const;
   
    /** 
     * \brief Creates a vertex containing the vertex data
     */
    void add_vertex(const vertex_id_type& vid, 
                    const VertexData& vdata = VertexData() );

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData());

  private:

    // HELPER ROUTINES =======================================================>    
    procid_t edge_to_proc(vertex_id_type source, vertex_id_type target) const {
      if(source > target) std::swap(source, target);
      boost::hash< std::pair<vertex_id_type, vertex_id_type> > hash_function;
      return hash_function(std::make_pair(source, target)) % rpc.numprocs();
    }

    procid_t is_local(vertex_id_type source, vertex_id_type target) const {
      return edge_to_proc(source, target) == rpc.procid();
    }

    procid_t vertex_to_init_proc(vertex_id_type vid) const { 
      return vid % rpc.numprocs();
    }    
    
    bool is_local_init(vertex_id_type vid) const {
      return vertex_to_init_proc(vid) == rpc.procid();
    }

    const vertex_record& vrecord(const vertex_id_type& vid) const {
      typedef typename vid2record_type::const_iterator iterator_type;
      iterator_type iter = vid2record.find(vid);
      ASSERT_TRUE(iter != vid2record.end());
      return iter->second;
    }


    vertex_record& vrecord(const vertex_id_type& vid) {
      typedef typename vid2record_type::iterator iterator_type;
      iterator_type iter = vid2record.find(vid);
      ASSERT_TRUE(iter != vid2record.end());
      return iter->second;
    }




    // Helper type used to synchronize the vertex data and assignments
    struct shuffle_record{
      procid_t owner;
      size_t num_in_edges, num_out_edges;
      std::vector<procid_t> mirrors;
      vertex_data_type vdata;
      shuffle_record() : 
        owner(0), num_in_edges(0), num_out_edges(0) { }
      void load(iarchive& arc) { 
        arc >> owner >> num_in_edges >> num_out_edges
            >> mirrors >> vdata;
      } // end of load
      void save(oarchive& arc) const { 
        arc << owner << num_in_edges << num_out_edges
            << mirrors << vdata;
      } // end of save     
    }; // end of vdata_shuffle_record;

    /// temporary map for vertexdata
    typedef boost::unordered_map<vertex_id_type, shuffle_record> vid2shuffle_type;
    vid2shuffle_type vid2shuffle;
    mutex vid2shuffle_lock;
    


  }; // End of graph











  //////////////////////////////////////////////////////////////////////////////
  //// Implementation 
  namespace distributed_graph_impl {
    struct preshuffle_record : public graphlab::IS_POD_TYPE {
      vertex_id_type vid, num_in_edges, num_out_edges;
      preshuffle_record() : vid(0), num_in_edges(0), num_out_edges(0) { }     
    };
  };


  template<typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::finalize()  {   
    rpc.full_barrier();
    // Check conditions on graph
    ASSERT_EQ(local_graph.num_vertices(), vid2record.size());
    // resize local vid map
    lvid2vid.resize(vid2record.size());
    
    using namespace distributed_graph_impl;
    // For all the vertices that this processor has seen determine the
    // "negotiator" and send that machine the negotiator.
    typedef std::vector< std::vector<preshuffle_record> > proc2vids_type;
    proc2vids_type proc2vids(rpc.numprocs());
    typedef typename vid2record_type::value_type vid2record_pair_type;
    foreach(vid2record_pair_type& pair, vid2record) {
      const vertex_id_type vid = pair.first;
      vertex_record& vrecord = pair.second;
      lvid2vid[vrecord.lvid] = vid;
      preshuffle_record pre_rec;
      pre_rec.vid = vid;
      pre_rec.num_in_edges = local_graph.num_in_edges(vrecord.lvid);
      pre_rec.num_out_edges = local_graph.num_out_edges(vrecord.lvid);
      proc2vids[vertex_to_init_proc(vid)].push_back(pre_rec);
    }
    // The returned local vertices are the vertices from each
    // machine for which this machine is a negotiator.
    mpi_tools::all2all(proc2vids, proc2vids);
    
   
    // Update the vid2shuffle
    for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
      foreach(const preshuffle_record& pre_rec, proc2vids[proc]) {
        shuffle_record& shuffle_rec = vid2shuffle[pre_rec.vid];
        shuffle_rec.num_in_edges += pre_rec.num_in_edges;
        shuffle_rec.num_out_edges += pre_rec.num_out_edges;
        shuffle_rec.mirrors.push_back(proc);
      }
    }

    // Construct the assignments
    std::vector<size_t> counts(rpc.numprocs());
    typedef typename vid2shuffle_type::value_type shuffle_pair_type;
    foreach(shuffle_pair_type& pair, vid2shuffle) {
      shuffle_record& record = pair.second;
      // Find the best (least loaded) processor to assign the vertex.
      std::pair<size_t, procid_t> 
        best_asg(counts[record.mirrors[0]], record.mirrors[0]);
      foreach(procid_t proc, record.mirrors)
        best_asg = std::min(best_asg, std::make_pair(counts[proc], proc));
      record.owner = best_asg.second;
      counts[record.owner]++;
    } // end of loop over 

    // Send the data to all the processors
    for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
      foreach(const preshuffle_record& pre_rec, proc2vids[proc]) {
        const std::pair<vertex_id_type, shuffle_record> 
          pair(pre_rec.vid, vid2shuffle[pre_rec.vid]);
        rpc.send_to(proc, pair);
      }
    }

    // Receive the data from all the processors.  Here we are a little
    // "clever" in that we loop over the vertices we have locally
    // managed and use them to determine how many times to recv_from which machines
    foreach(const vertex_id_type& tmpvid, lvid2vid) {     
      const procid_t proc = vertex_to_init_proc(tmpvid);
      std::pair<vertex_id_type, shuffle_record> vid_and_rec;
      rpc.recv_from(proc, vid_and_rec);
      const vertex_id_type& vid = vid_and_rec.first;
      shuffle_record& shuffle_rec = vid_and_rec.second;      
      vertex_record& vrecord = vid2record[vid];
      vrecord.mirrors.swap(shuffle_rec.mirrors);
      vrecord.owner = shuffle_rec.owner;
      local_graph.vertex_data(vrecord.lvid) = shuffle_rec.vdata;
      vrecord.num_in_edges = shuffle_rec.num_in_edges;
      vrecord.num_out_edges = shuffle_rec.num_out_edges;
    }




    // // Receive assignments from coordinators
    // mpi_tools::all2all(vdata_shuffle, vdata_shuffle);
    
    // // Incorporate all the vertex data assigned to each machine
    // for(size_t i = 0; i < vdata_shuffle.size(); ++i) {
    //   foreach(vdata_shuffle_record& shuffle_record, vdata_shuffle[i]) {
    //     ASSERT_TRUE(vid2record.find(shuffle_record.vid) != vid2record.end());
    //     vertex_record& record = vid2record[shuffle_record.vid];
    //     record.owner = shuffle_record.owner; 
    //     record.mirrors.swap(shuffle_record.mirrors);
    //     local_graph.vertex_data(record.lvid) = shuffle_record.vdata;
    //   } // end of loop over vdata
    // } // end of loop over sending machines
   
    // std::cout << "Save debugging information" << std::endl;
    // {
    //   const std::string fname = 
    //     "file_" + boost::lexical_cast<std::string>(rpc.procid());
    //   std::ofstream fout(fname.c_str());
    //   typedef typename vid2record_type::value_type pair_type;
    //   foreach(const pair_type& pair, vid2record) {      
    //     fout << pair.first << '\t' << pair.second.owner << '\t';
    //     std::vector<bool> bitset(rpc.numprocs(), false);
    //     foreach(const procid_t& proc, pair.second.mirrors)
    //       bitset[proc] = true;
    //     for(size_t i = 0; i < bitset.size(); ++i) {
    //       fout << (bitset[i]? '1' : '0') 
    //            << (i+1 < bitset.size()? '\t' : '\n');
    //     }
    //   }
    //   fout.close();
    // }   

  } // End of finalize
  
  
  template<typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  add_vertex(const vertex_id_type& vid,
             const VertexData& vdata) {
    // determine if the vertex is local
    if(is_local_init(vid)) {
      vid2shuffle_lock.lock();
      vid2shuffle[vid].vdata = vdata;
      vid2shuffle_lock.unlock();
    } else {
      rpc.remote_call(vertex_to_init_proc(vid),
                      &distributed_graph::add_vertex,
                      vid, vdata);
    }
  } // End of add vertex;
  
  
  template<typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>:: 
  add_edge(vertex_id_type source, vertex_id_type target, 
           const EdgeData& edata) {
    // determine if the edge is locally managed
    if(is_local(source,target)) {
      // get (or create) the local ids for source and target 
      vid2record_lock.lock();
      // if the record was just created (it will have an invalid
      // lvid=-1) then assign it the next lvid
      vertex_record& source_rec = vid2record[source];
      if(source_rec.lvid == vertex_id_type(-1))
        source_rec.lvid = vid2record.size() - 1;
      vertex_record& target_rec = vid2record[target];
      if(target_rec.lvid == vertex_id_type(-1))
        target_rec.lvid = vid2record.size() - 1;
      vid2record_lock.unlock();
      // Add the record to the local graph
      local_graph_lock.lock();
      const vertex_id_type max_lvid = std::max(source_rec.lvid, target_rec.lvid);
      if(max_lvid >= local_graph.num_vertices()) 
        local_graph.resize(max_lvid + 1);
      local_graph.add_edge(source_rec.lvid, target_rec.lvid, edata);
      local_graph_lock.unlock();
    } else {      
      rpc.remote_call(edge_to_proc(source, target),
                      &distributed_graph::add_edge,
                      source, target, edata);
    }
  } // End of add edge



  
  template<typename VertexData, typename EdgeData>
  typename distributed_graph<VertexData, EdgeData>::edge_list_type 
  distributed_graph<VertexData, EdgeData>:: 
  in_edges(const vertex_id_type vid) const {
    return edge_list_type(this, local_graph.in_edges(vrecord(vid).lvid));
  } // end of in_edges
  

  template<typename VertexData, typename EdgeData>
  size_t distributed_graph<VertexData, EdgeData>:: 
  num_in_edges(const vertex_id_type v) const {
    return vrecord(v).num_in_edges;
  } // end of num_in_edges
  
  template<typename VertexData, typename EdgeData>
  typename distributed_graph<VertexData, EdgeData>::edge_list_type 
  distributed_graph<VertexData, EdgeData>:: 
  out_edges(const vertex_id_type vid) const { 
    return edge_list_type(this, local_graph.out_edges(vrecord(vid).lvid));
  } // end of out_edges
  
  template<typename VertexData, typename EdgeData>
  size_t distributed_graph<VertexData, EdgeData>:: 
  num_out_edges(const vertex_id_type v) const {
    return vrecord(v).num_out_edges;
  } // end of num out edges





  template<typename VertexData, typename EdgeData>
  VertexData& distributed_graph<VertexData, EdgeData>:: 
  vertex_data(vertex_id_type vid) {
    return local_graph.vertex_data(vrecord(vid).lvid);
  } // end of vertex data

    

  template<typename VertexData, typename EdgeData>
  const VertexData& distributed_graph<VertexData, EdgeData>:: 
  vertex_data(vertex_id_type vid) const {
    return local_graph.vertex_data(vrecord(vid).lvid);
  } // end of const vertex data


  template<typename VertexData, typename EdgeData>
  EdgeData& distributed_graph<VertexData, EdgeData>:: 
  edge_data(vertex_id_type source, vertex_id_type target) {
    return local_graph.edge_data(vrecord(source).lvid, vrecord(target).lvid);
  } // end of edge data


  template<typename VertexData, typename EdgeData>
  const EdgeData& distributed_graph<VertexData, EdgeData>:: 
  edge_data(vertex_id_type source, vertex_id_type target) const {
    return local_graph.edge_data(vrecord(source).lvid, vrecord(target).lvid);
  } // end of const edge data

  template<typename VertexData, typename EdgeData>
  EdgeData& distributed_graph<VertexData, EdgeData>:: 
  edge_data(const edge_type& edge) {
    return local_graph.edge_data(edge.local_edge);
  } // end of edge data


  template<typename VertexData, typename EdgeData>
  const EdgeData& distributed_graph<VertexData, EdgeData>:: 
  edge_data(const edge_type& edge) const {
    return local_graph.edge_data(edge.local_edge);
  } // end of const edge data





} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

