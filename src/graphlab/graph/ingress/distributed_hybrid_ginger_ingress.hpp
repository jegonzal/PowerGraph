/**
 * Copyright (c) 2013 Institute of Parallel and Distributed Systems, Shanghai Jiao Tong University.
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
 *      http://ipads.se.sjtu.edu.cn/projects/powerlyra.html
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_HYBRID_GINGER_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_HYBRID_GINGER_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/graph_hash.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/logger/logger.hpp>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <graphlab/macros_def.hpp>

#undef TUNING
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges using a hybrid method.
   *        That is, for high degree edges, 
   *        for low degree edges, hashing from its target vertex.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_hybrid_ginger_ingress : 
    public distributed_ingress_base<VertexData, EdgeData> {

  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    typedef distributed_ingress_base<VertexData, EdgeData> base_type;

    typedef typename graph_type::vertex_record vertex_record;
    typedef typename graph_type::mirror_type mirror_type;


    typedef typename buffered_exchange<vertex_id_type>::buffer_type
        vertex_id_buffer_type;

    /// one-by-one ingress. e.g., SNAP
    typedef typename base_type::edge_buffer_record edge_buffer_record;
    typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
            edge_buffer_type;

    typedef typename base_type::vertex_buffer_record vertex_buffer_record;
    typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;

    /// batch ingress. e.g., R-ADJ
    struct batch_edge_buffer_record {
      std::vector<vertex_id_type> sources;
      vertex_id_type target;
      std::vector<edge_data_type> edatas;
      
      batch_edge_buffer_record(
        const std::vector<vertex_id_type>& sources = std::vector<vertex_id_type>() , 
        const vertex_id_type& target = vertex_id_type(-1), 
        const std::vector<edge_data_type>& edatas = std::vector<edge_data_type>()) :
        sources(sources), target(target), edatas(edatas) { }

      void load(iarchive& arc) { arc >> sources >> target >> edatas; }
      void save(oarchive& arc) const { arc << sources << target << edatas; }
    };
    typedef typename buffered_exchange<batch_edge_buffer_record>::buffer_type 
        batch_edge_buffer_type;
    
    /// detail vertex record for the second pass coordination. 
    typedef typename base_type::vertex_negotiator_record 
      vertex_negotiator_record;

    /// ginger structure
    /** Type of the master location hash table: [vertex-id, location-of-master] */
    typedef typename boost::unordered_map<vertex_id_type, procid_t> 
        master_hash_table_type;
    typedef typename std::pair<vertex_id_type, procid_t> 
        master_pair_type;
    typedef typename buffered_exchange<master_pair_type>::buffer_type 
        master_buffer_type;

    typedef typename std::pair<procid_t,size_t> 
        proc_edges_pair_type;
    typedef typename buffered_exchange<proc_edges_pair_type >::buffer_type
        proc_edges_buffer_type;

    /* ingress edges */
    buffered_exchange<edge_buffer_record> single_edge_exchange;
    buffered_exchange<edge_buffer_record> high_edge_exchange;
    buffered_exchange<edge_buffer_record>  low_edge_exchange;
    /* ingress vertex data */
    buffered_exchange<vertex_buffer_record> hybrid_vertex_exchange;
    buffered_exchange<vertex_buffer_record> temporary_vertex_exchange;

    /** master hash table:
     * mht is the synchronized one across the cluster,
     * mht_incremental is the new added mht which stays local since the last sync.
     */
    master_hash_table_type mht;
    master_hash_table_type mht_incremental;

    buffered_exchange<master_pair_type> master_exchange;

    /// The rpc interface for this object
    dc_dist_object<distributed_hybrid_ginger_ingress> hybrid_rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;

    /// threshold to divide high-degree and low-degree vertices
    size_t threshold;

    bool standalone;

    std::vector<edge_buffer_record> hybrid_edges;
    simple_spinlock hybrid_edges_lock;

    //these are used for edges balance
    buffered_exchange< proc_edges_pair_type >  proc_edges_exchange;
    std::vector<size_t> proc_edges_incremental;
    std::vector<size_t> proc_num_vertices;

    /// heuristic model from fennel
    /// records about the number of edges and vertices in the graph
    /// given from the commandline
    size_t nedges;
    size_t nverts;
    /// threshold for incremental mht to be synced across the cluster
    /// when the incremental mht size reaches the preset interval,
    /// we will perform a synchronization on mht across the cluster
    size_t interval;
    /// arguments for the ginger algorithm
    double alpha;
    double gamma;


  public:
    distributed_hybrid_ginger_ingress(distributed_control& dc, graph_type& graph, 
        size_t threshold = 100, size_t nedges = 0, size_t nverts = 0, 
        size_t interval = std::numeric_limits<size_t>::max()) :
        base_type(dc, graph), 
        single_edge_exchange(dc), high_edge_exchange(dc), low_edge_exchange(dc),
        hybrid_vertex_exchange(dc), temporary_vertex_exchange(dc), master_exchange(dc),
        hybrid_rpc(dc, this), graph(graph), threshold(threshold), 
        proc_edges_exchange(dc), proc_edges_incremental(dc.numprocs()), 
        proc_num_vertices(dc.numprocs()), 
        nedges(nedges), nverts(nverts), interval(interval) { 
      ASSERT_GT(nedges, 0); ASSERT_GT(nverts, 0);
      
      gamma = 1.5;
      alpha = sqrt(dc.numprocs()) * nedges / pow(nverts, gamma);

      /* fast pass for standalone case. */
      standalone = hybrid_rpc.numprocs() == 1;
      hybrid_rpc.barrier();
    } // end of constructor

    ~distributed_hybrid_ginger_ingress() { }

    /** Add an edge to the ingress object using random hashing assignment.
     *  This function acts as the first phase for SNAP graph to deliver edges
     *  via the hashing value of its target vertex.
     */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      const edge_buffer_record record(source, target, edata);
      if (standalone) {
        /* Fast pass for standalone case. */
        hybrid_edges_lock.lock();
        hybrid_edges.push_back(record);
        hybrid_edges_lock.unlock();
      } else {
        const procid_t owning_proc = 
          graph_hash::hash_vertex(target) % hybrid_rpc.numprocs();
        single_edge_exchange.send(owning_proc, record);
      }
    } // end of add edge

    /** Add edges to the ingress object using different assignment policies.
     *  This function handles the RADJ graph in different ways.
     *  For high degree edges, hashing from its source vertex;
     *  for low degree edges, using a heuristic way named ginger.
     */
    void add_edges(std::vector<vertex_id_type>& in, vertex_id_type out,
                      const std::vector<EdgeData>& edatas) {
      if (standalone) {
        /* fast pass for standalone case. */
        hybrid_edges_lock.lock();
        for(size_t i = 0; i < in.size();i++){
          const edge_buffer_record record(in[i], out, edatas[i]);
          hybrid_edges.push_back(record);
        }
        hybrid_edges_lock.unlock();
      } else {
        procid_t owning_proc = 0;
        if (in.size() >= threshold) {
          // TODO: no need send, just resend latter
          owning_proc = graph_hash::hash_vertex(out) % base_type::rpc.numprocs();
          for (size_t i = 0; i < in.size(); ++i) {
            edge_buffer_record record(in[i], out, edatas[i]);
            high_edge_exchange.send(owning_proc, record);
          }
        } else {
          owning_proc = base_type::edge_decision.edge_to_proc_ginger(
                        out, mht, mht_incremental, proc_num_vertices, 
                        alpha, gamma, in);
          for (size_t i = 0; i < in.size(); ++i) {
            edge_buffer_record record(in[i], out, edatas[i]);
            low_edge_exchange.send(owning_proc, record);
          }
          
          // update local counter
          proc_edges_incremental[owning_proc] += in.size();
        }

        // update mapping table for ginger
        size_t numprocs = base_type::rpc.numprocs();
        mht_incremental[out] = owning_proc;

        // broadcast local counters
        if (mht_incremental.size() > interval) {
          //update mht
          for (typename master_hash_table_type::iterator it = mht_incremental.begin();
                  it != mht_incremental.end(); ++it) {
            for (procid_t i = 0; i < numprocs; ++i) {
              if (i != hybrid_rpc.procid())
                master_exchange.send(i, master_pair_type(it->first, it->second));
            }
            mht[it->first] = it->second;
          }
          master_exchange.partial_flush(0);
          mht_incremental.clear();
          master_buffer_type master_buffer;
          procid_t proc;
          while(master_exchange.recv(proc, master_buffer)) {
            foreach(const master_pair_type& pair, master_buffer) {
              mht[pair.first] = pair.second;
              //proc_num_vertices[pair.second]++;
            }
          }
          
          //update proc_edges for edges balance
          for (procid_t p=0; p<hybrid_rpc.numprocs(); p++) {
            for (procid_t i=0; i<hybrid_rpc.numprocs(); i++) {
              proc_edges_exchange.send(p,std::make_pair(i,proc_edges_incremental[i]));
            }
          }
          proc_edges_exchange.partial_flush(0);
          
          proc_edges_buffer_type proc_edge_buffer;
          while (proc_edges_exchange.recv(proc, proc_edge_buffer)) {
            foreach (const proc_edges_pair_type& pair, proc_edge_buffer) {
              proc_edges_incremental[pair.first] += pair.second;
            }
          }
          proc_edges_exchange.clear();

          for (int i=0; i<hybrid_rpc.numprocs(); i++) {
            proc_num_vertices[i] += proc_edges_incremental[i]*1.0*nverts/nedges;
            proc_edges_incremental[i] = 0;
          }
        }
      }
    } // end of add edges

    /* add vdata */
    void add_vertex(vertex_id_type vid, const VertexData& vdata) { 
      const vertex_buffer_record record(vid, vdata);
      if (standalone) {
        /* fast pass for redundant finalization with no graph changes. */
        hybrid_vertex_exchange.send(0, record);
      } else {
        const procid_t owning_proc = 
          graph_hash::hash_vertex(vid) % hybrid_rpc.numprocs();
        temporary_vertex_exchange.send(owning_proc, record);
      }
    } // end of add vertex


    void assign_edges() {
      single_edge_exchange.flush();
      if (single_edge_exchange.size() > 0) {        
        typedef typename boost::unordered_map<vertex_id_type, batch_edge_buffer_record> 
                  batch_record_map_type;
        batch_record_map_type batch_map;

        edge_buffer_type edge_buffer;
        procid_t proc = -1;
        while (single_edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            batch_map[rec.target].sources.push_back(rec.source);
            batch_map[rec.target].edatas.push_back(rec.edata);
          }
        }
        // re-assigne
        for (typename batch_record_map_type::iterator it = batch_map.begin(); 
            it != batch_map.end(); ++it) {
          add_edges(it->second.sources, it->first,it->second.edatas);
        }

        single_edge_exchange.clear();
        //hybrid_rpc.full_barrier();
      }
    }

    void finalize() {
      graphlab::timer ti;
      
      size_t nprocs = hybrid_rpc.numprocs();
      procid_t l_procid = hybrid_rpc.procid();
      size_t nedges = 0;

      if (l_procid == 0) {
        memory_info::log_usage("start finalizing");
        logstream(LOG_EMPH) << "ginger finalizing ..."
                            << " #vertices=" << graph.local_graph.num_vertices()
                            << " #edges=" << graph.local_graph.num_edges()
                            << " threshold=" << threshold
                            << std::endl;
      }

      /**************************************************************************/
      /*                                                                        */
      /*               Assign edges for one-by-one format              */
      /*                                                                        */
      /**************************************************************************/
      if (!standalone) assign_edges();

      /**
       * Fast pass for redundant finalization with no graph changes. 
       */
      size_t changed_size;
      
      // flush and count changes
      if (standalone) {
        nedges = hybrid_edges.size();
        hybrid_vertex_exchange.flush();
        
        changed_size = nedges + hybrid_vertex_exchange.size();
      } else {
        high_edge_exchange.flush(); 
        low_edge_exchange.flush();
        temporary_vertex_exchange.flush();

        changed_size = high_edge_exchange.size() 
                     + low_edge_exchange.size() 
                     + hybrid_vertex_exchange.size();
      }
      
      hybrid_rpc.all_reduce(changed_size);
      if (changed_size == 0) {
        logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
        return;
      }
      
      /**************************************************************************/
      /*                                                                        */
      /*                       prepare hybrid ingress                           */
      /*                                                                        */
      /**************************************************************************/
      if (!standalone) {
        /* exchange mapping table*/
        for (typename master_hash_table_type::iterator it = mht_incremental.begin();
                 it != mht_incremental.end(); ++it) {
          for (procid_t i = 0; i < nprocs; ++i) {
            if (i != l_procid)
              master_exchange.send(i, master_pair_type(it->first, it->second));
          }
          mht[it->first] = it->second;
        }
        mht_incremental.clear();

        master_exchange.flush();
        master_buffer_type master_buffer;
        procid_t proc = -1;
        while(master_exchange.recv(proc, master_buffer)) {
          foreach(const master_pair_type& pair, master_buffer) {
            mht[pair.first] = pair.second;
          }
        }
        master_exchange.clear();

#ifdef TUNING
        if(l_procid == 0) { 
          memory_info::log_usage("exchange mapping done.");
          logstream(LOG_EMPH) << "exchange mapping: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif


        /* collect edges for batch ingress  (e.g. RADJ) */
        nedges = low_edge_exchange.size();
        hybrid_edges.reserve(nedges + high_edge_exchange.size());

        edge_buffer_type edge_buffer;
        proc = -1;
        while(low_edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            if (mht.find(rec.source) == mht.end())
              mht[rec.source] = graph_hash::hash_vertex(rec.source) % nprocs;   
            hybrid_edges.push_back(rec);
          }
        }
        low_edge_exchange.clear();
        hybrid_rpc.full_barrier(); // sync before reusing

        // re-send edges of high-degree vertices
        proc = -1;
        while(high_edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            if (mht.find(rec.source) == mht.end())
              mht[rec.source] = graph_hash::hash_vertex(rec.source) % nprocs; 
            
            const procid_t source_owner_proc = mht[rec.source];
            if (source_owner_proc == l_procid) {
              hybrid_edges.push_back(rec);
              ++nedges;
            } else {
              low_edge_exchange.send(source_owner_proc, rec);
            }
          }
        }
        high_edge_exchange.clear();

        // receive edges of high-degree vertices
        low_edge_exchange.flush();        
        proc = -1;
        while(low_edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            mht[rec.source] = l_procid;
            hybrid_edges.push_back(rec);
            ++nedges;
          }
        }
        low_edge_exchange.clear();
      }

      if(l_procid == 0) {
        memory_info::log_usage("prepare ginger finalizing done.");
        logstream(LOG_EMPH) << "prepare ginger finalizing. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }

      // connect to base finalize()
      modified_base_finalize(nedges);
      if(l_procid == 0) {
        memory_info::log_usage("base finalizing done.");
        logstream(LOG_EMPH) << "base finalizing. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }

      set_vertex_type();
      if(l_procid == 0) {
        memory_info::log_usage("set vertex type done.");
        logstream(LOG_EMPH) << "set vertex type. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }
      
      if(l_procid == 0)
        logstream(LOG_EMPH) << "ginger finalizing graph done. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
    } // end of finalize

    void set_vertex_type() {
      graphlab::timer ti;
      procid_t l_procid = hybrid_rpc.procid();
      size_t high_master = 0, high_mirror = 0, low_master = 0, low_mirror = 0;

      for (size_t lvid = 0; lvid < graph.num_local_vertices(); lvid++) {
        vertex_record& vrec = graph.lvid2record[lvid];
        if (vrec.num_in_edges >= threshold) {
          if (vrec.owner == l_procid) {
            vrec.type = graph_type::HIGH_MASTER; 
            high_master ++;
          } else {
            vrec.type = graph_type::HIGH_MIRROR;
            high_mirror ++;
          }
        } else {
          if (vrec.owner == l_procid) {
            vrec.type = graph_type::LOW_MASTER; 
            low_master ++;
          } else {
            vrec.type = graph_type::LOW_MIRROR;
            low_mirror ++;
          }
        }        
      }

//#ifdef TUNING
      // Compute the total number of high-degree and low-degree vertices
      std::vector<size_t> swap_counts(hybrid_rpc.numprocs());

      swap_counts[l_procid] = high_master;
      hybrid_rpc.all_gather(swap_counts);
      high_master = 0;
      foreach(size_t count, swap_counts) high_master += count;

      swap_counts[l_procid] = high_mirror;
      hybrid_rpc.all_gather(swap_counts);
      high_mirror = 0;
      foreach(size_t count, swap_counts) high_mirror += count;

      swap_counts[l_procid] = low_master;
      hybrid_rpc.all_gather(swap_counts);
      low_master = 0;
      foreach(size_t count, swap_counts) low_master += count;

      swap_counts[l_procid] = low_mirror;
      hybrid_rpc.all_gather(swap_counts);
      low_mirror = 0;
      foreach(size_t count, swap_counts) low_mirror += count;

      if(l_procid == 0) {
        logstream(LOG_EMPH) << "hybrid info: master [" 
                            << high_master << " " 
                            << low_master << " " 
                            << ((high_master*1.0)/(high_master+low_master)) << "]"
                            << std::endl;
        if ((high_mirror + low_mirror) > 0)
        logstream(LOG_EMPH) << "hybrid info: mirror [" 
                            << high_mirror << " " 
                            << low_mirror << " " 
                            << ((high_mirror*1.0)/(high_mirror+low_mirror)) << "]"
                            << std::endl;

        memory_info::log_usage("set vertex type done."); 
        logstream(LOG_EMPH) << "set vertex type: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
//#endif
    }


    /* do the same job as original base finalize except for
     * extracting edges from hybrid_edges instead of original edge_buffer;
     * and using mht to tracing the master location of each vertex.
     */
    void modified_base_finalize(size_t nedges) {
      graphlab::timer ti;
      procid_t l_procid = hybrid_rpc.procid();
      size_t nprocs = hybrid_rpc.numprocs();

      hybrid_rpc.full_barrier();
      
      bool first_time_finalize = false;
      /**
       * Fast pass for first time finalization. 
       */
      if (graph.is_dynamic()) {
        size_t nverts = graph.num_local_vertices();
        hybrid_rpc.all_reduce(nverts);
        first_time_finalize = (nverts == 0);
      } else {
        first_time_finalize = false;
      }

      
      typedef typename hopscotch_map<vertex_id_type, lvid_type>::value_type
          vid2lvid_pair_type;

      /**
       * \internal
       * Buffer storage for new vertices to the local graph.
       */
      typedef typename graph_type::hopscotch_map_type vid2lvid_map_type;
      vid2lvid_map_type vid2lvid_buffer;

      /**
       * \internal
       * The begining id assinged to the first new vertex.
       */
      const lvid_type lvid_start  = graph.vid2lvid.size();

      /**
       * \internal
       * Bit field incidate the vertex that is updated during the ingress. 
       */
      dense_bitset updated_lvids(graph.vid2lvid.size());


      /**************************************************************************/
      /*                                                                        */
      /*                         Construct local graph                          */
      /*                                                                        */
      /**************************************************************************/
      { 
        // Add all the edges to the local graph
        graph.local_graph.reserve_edge_space(nedges + 1);

        foreach(const edge_buffer_record& rec, hybrid_edges) {
          // Get the source_vlid;
          lvid_type source_lvid(-1);
          if(graph.vid2lvid.find(rec.source) == graph.vid2lvid.end()) {
            if (vid2lvid_buffer.find(rec.source) == vid2lvid_buffer.end()) {
              source_lvid = lvid_start + vid2lvid_buffer.size();
              vid2lvid_buffer[rec.source] = source_lvid;
            } else {
              source_lvid = vid2lvid_buffer[rec.source];
            }
          } else {
            source_lvid = graph.vid2lvid[rec.source];
            updated_lvids.set_bit(source_lvid);
          }
          // Get the target_lvid;
          lvid_type target_lvid(-1);
          if(graph.vid2lvid.find(rec.target) == graph.vid2lvid.end()) {
            if (vid2lvid_buffer.find(rec.target) == vid2lvid_buffer.end()) {                
              target_lvid = lvid_start + vid2lvid_buffer.size();
              vid2lvid_buffer[rec.target] = target_lvid;
            } else {
              target_lvid = vid2lvid_buffer[rec.target];
            }
          } else {
            target_lvid = graph.vid2lvid[rec.target];
            updated_lvids.set_bit(target_lvid);
          }
          graph.local_graph.add_edge(source_lvid, target_lvid, rec.edata);
        } // end for loop over buffers
        hybrid_edges.clear();

        ASSERT_EQ(graph.vid2lvid.size() + vid2lvid_buffer.size(), 
                  graph.local_graph.num_vertices());
#ifdef TUNING
        if(l_procid == 0)  {
          memory_info::log_usage("base::populating local graph done.");
          logstream(LOG_EMPH) << "base::populating local graph: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif
        // Finalize local graph
        graph.local_graph.finalize();
#ifdef TUNING
        logstream(LOG_INFO) << "base::local graph info: " << std::endl
                            << "\t nverts: " << graph.local_graph.num_vertices()
                            << std::endl
                            << "\t nedges: " << graph.local_graph.num_edges()
                            << std::endl;
        
        if(l_procid == 0) {
          memory_info::log_usage("base::finalizing local graph done."); 
          logstream(LOG_EMPH) << "base::finalizing local graph: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif
      }


      /**************************************************************************/
      /*                                                                        */
      /*             Receive and add vertex data to masters                     */
      /*                                                                        */
      /**************************************************************************/
      // Setup the map containing all the vertices being negotiated by this machine
      {
        if (!standalone) {
          // has flushed in finalize()
          if (temporary_vertex_exchange.size() > 0) {
            vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
            while (temporary_vertex_exchange.recv(sending_proc, vertex_buffer)) {
              foreach (const vertex_buffer_record& rec, vertex_buffer) {
                if (mht.find(rec.vid) == mht.end())
                  mht[rec.vid] = graph_hash::hash_vertex(rec.vid) % nprocs; 
                hybrid_vertex_exchange.send(mht[rec.vid], rec);
              }
            }
            temporary_vertex_exchange.clear();
          }

          // match standalone code
          hybrid_vertex_exchange.flush();
        }

        // receive any vertex data sent by other machines        
        if (hybrid_vertex_exchange.size() > 0) {
          vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
          while(hybrid_vertex_exchange.recv(sending_proc, vertex_buffer)) {
            foreach(const vertex_buffer_record& rec, vertex_buffer) {
              lvid_type lvid(-1);
              if (graph.vid2lvid.find(rec.vid) == graph.vid2lvid.end()) {
                if (vid2lvid_buffer.find(rec.vid) == vid2lvid_buffer.end()) {
                  lvid = lvid_start + vid2lvid_buffer.size();
                  vid2lvid_buffer[rec.vid] = lvid;
                } else {
                  lvid = vid2lvid_buffer[rec.vid];
                }
              } else {
                lvid = graph.vid2lvid[rec.vid];
                updated_lvids.set_bit(lvid);
              }
              if (distributed_hybrid_ginger_ingress::vertex_combine_strategy 
                  && lvid < graph.num_local_vertices()) {
                distributed_hybrid_ginger_ingress::vertex_combine_strategy(
                  graph.l_vertex(lvid).data(), rec.vdata);
              } else {
                graph.local_graph.add_vertex(lvid, rec.vdata);
              }
            }
          }
          hybrid_vertex_exchange.clear();
#ifdef TUNING
          logstream(LOG_INFO) << "base::#vert-msgs=" << hybrid_vertex_exchange.size()
                              << std::endl;
          if(l_procid == 0) {        
            memory_info::log_usage("base::adding vertex data done.");
            logstream(LOG_EMPH) << "base::adding vertex data: " 
                                << ti.current_time()
                                << " secs" 
                                << std::endl;
          }
#endif
        }
      } // end of loop to populate vrecmap


      /**************************************************************************/
      /*                                                                        */
      /*        Assign vertex data and allocate vertex (meta)data  space        */
      /*                                                                        */
      /**************************************************************************/
      {
        // determine masters for all negotiated vertices
        const size_t local_nverts = graph.vid2lvid.size() + vid2lvid_buffer.size();
        graph.lvid2record.reserve(local_nverts);
        graph.lvid2record.resize(local_nverts);
        graph.local_graph.resize(local_nverts);
        foreach(const vid2lvid_pair_type& pair, vid2lvid_buffer) {
          vertex_record& vrec = graph.lvid2record[pair.second];
          vrec.gvid = pair.first;
          if (standalone) {
            vrec.owner = 0;
          } else {
            if (mht.find(pair.first) == mht.end())
                mht[pair.first] = graph_hash::hash_vertex(pair.first) % nprocs; 
            const procid_t source_owner_proc = mht[pair.first];
            vrec.owner = source_owner_proc;
          }
        }
        ASSERT_EQ(local_nverts, graph.local_graph.num_vertices());
        ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
#ifdef TUNING
        if(l_procid == 0) {
          memory_info::log_usage("allocating lvid2record done.");
          logstream(LOG_EMPH) << "allocating lvid2record: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif
        mht.clear();
      }

      /**************************************************************************/
      /*                                                                        */
      /*                          Master handshake                              */
      /*                                                                        */
      /**************************************************************************/      
      if (!standalone) {
#ifdef _OPENMP
        buffered_exchange<vertex_id_type> vid_buffer(hybrid_rpc.dc(), omp_get_max_threads());
#else
        buffered_exchange<vertex_id_type> vid_buffer(hybrid_rpc.dc());
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
        // send not owned vids to their master
        for (lvid_type i = lvid_start; i < graph.lvid2record.size(); ++i) {
          procid_t master = graph.lvid2record[i].owner;
          if (master != l_procid)
#ifdef _OPENMP
            vid_buffer.send(master, graph.lvid2record[i].gvid, omp_get_thread_num());
#else
            vid_buffer.send(master, graph.lvid2record[i].gvid);
#endif
        }
        vid_buffer.flush();
        hybrid_rpc.barrier();

        // receive all vids owned by me
        mutex flying_vids_lock;
        boost::unordered_map<vertex_id_type, mirror_type> flying_vids;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          typename buffered_exchange<vertex_id_type>::buffer_type buffer;
          procid_t recvid = -1;
          while(vid_buffer.recv(recvid, buffer)) {
            foreach(const vertex_id_type vid, buffer) {
              if (graph.vid2lvid.find(vid) == graph.vid2lvid.end()) {
                if (vid2lvid_buffer.find(vid) == vid2lvid_buffer.end()) {
                  flying_vids_lock.lock();
                  mirror_type& mirrors = flying_vids[vid];
                  mirrors.set_bit(recvid);
                  flying_vids_lock.unlock();
                } else {
                  lvid_type lvid = vid2lvid_buffer[vid];
                  graph.lvid2record[lvid]._mirrors.set_bit(recvid);
                }
              } else {
                lvid_type lvid = graph.vid2lvid[vid];
                graph.lvid2record[lvid]._mirrors.set_bit(recvid);
                updated_lvids.set_bit(lvid);
              }
            }
          }
        }
        vid_buffer.clear();

        if (!flying_vids.empty()) {
          logstream(LOG_INFO) << "base::#flying-own-nverts="
                              << flying_vids.size() 
                              << std::endl;

          // reallocate spaces for the flying vertices. 
          size_t vsize_old = graph.lvid2record.size();
          size_t vsize_new = vsize_old + flying_vids.size();
          graph.lvid2record.resize(vsize_new);
          graph.local_graph.resize(vsize_new);
          for (typename boost::unordered_map<vertex_id_type, mirror_type>::iterator it = flying_vids.begin();
               it != flying_vids.end(); ++it) {
            lvid_type lvid = lvid_start + vid2lvid_buffer.size();
            vertex_record& vrec = graph.lvid2record[lvid];
            vertex_id_type gvid = it->first;
            vrec.owner = l_procid;
            vrec.gvid = gvid;
            vrec._mirrors = it->second;
            vid2lvid_buffer[gvid] = lvid;
          }
        }
      } // end of master handshake

#ifdef TUNING
      if(l_procid == 0) {
        memory_info::log_usage("base::master handshake done.");
        logstream(LOG_EMPH) << "base::master handshake: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      /**************************************************************************/
      /*                                                                        */
      /*                        Merge in vid2lvid_buffer                        */
      /*                                                                        */
      /**************************************************************************/
      {
        if (graph.vid2lvid.size() == 0) {
          graph.vid2lvid.swap(vid2lvid_buffer);
        } else {
          graph.vid2lvid.rehash(graph.vid2lvid.size() + vid2lvid_buffer.size());
          foreach (const typename vid2lvid_map_type::value_type& pair, vid2lvid_buffer) {
            graph.vid2lvid.insert(pair);
          }
          vid2lvid_buffer.clear();
        }
      }


      /**************************************************************************/
      /*                                                                        */
      /*              Synchronize vertex data and meta information              */
      /*                                                                        */
      /**************************************************************************/
      // TODO:  optimization for standalone
      {
        // construct the vertex set of changed vertices

        // Fast pass for first time finalize;
        vertex_set changed_vset(true);

        // Compute the vertices that needs synchronization 
        if (!first_time_finalize) {
          vertex_set changed_vset = vertex_set(false);
          changed_vset.make_explicit(graph);
          updated_lvids.resize(graph.num_local_vertices());
          for (lvid_type i = lvid_start; i <  graph.num_local_vertices(); ++i) {
            updated_lvids.set_bit(i);
          }
          changed_vset.localvset = updated_lvids; 
          buffered_exchange<vertex_id_type> vset_exchange(hybrid_rpc.dc());
          // sync vset with all mirrors
          changed_vset.synchronize_mirrors_to_master_or(graph, vset_exchange);
          changed_vset.synchronize_master_to_mirrors(graph, vset_exchange);
        }

        graphlab::graph_gather_apply<graph_type, vertex_negotiator_record> 
            vrecord_sync_gas(graph, 
                             boost::bind(&distributed_hybrid_ginger_ingress::finalize_gather, this, _1, _2), 
                             boost::bind(&distributed_hybrid_ginger_ingress::finalize_apply, this, _1, _2, _3));
        vrecord_sync_gas.exec(changed_vset);

        if(l_procid == 0) {
          memory_info::log_usage("synchronizing vertex (meta)data done.");
          logstream(LOG_EMPH) << "synchrionizing vertex (meta)data: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
      }

      base_type::exchange_global_info(standalone);

#ifdef TUNING
      if(l_procid == 0) {
        memory_info::log_usage("exchange global info done.");
        logstream(LOG_EMPH) << "exchange global info: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif
    } // end of modified base finalize

  private:
    boost::function<void(vertex_data_type&, const vertex_data_type&)> vertex_combine_strategy;

    /**
     * \brief Gather the vertex distributed meta data.
     */
    vertex_negotiator_record finalize_gather(lvid_type& lvid, graph_type& graph) {
      vertex_negotiator_record accum;
      accum.num_in_edges = graph.local_graph.num_in_edges(lvid);
      accum.num_out_edges = graph.local_graph.num_out_edges(lvid);
      if (graph.l_is_master(lvid)) {
        accum.has_data = true;
        accum.vdata = graph.l_vertex(lvid).data();
        accum.mirrors = graph.lvid2record[lvid]._mirrors;
      }
      return accum;
    }

    /**
     * \brief Update the vertex data structures with the gathered vertex metadata.  
     */
    void finalize_apply(lvid_type lvid, const vertex_negotiator_record& accum, graph_type& graph) {
      typename graph_type::vertex_record& vrec = graph.lvid2record[lvid];
      vrec.num_in_edges = accum.num_in_edges;
      vrec.num_out_edges = accum.num_out_edges;
      graph.l_vertex(lvid).data() = accum.vdata;
      vrec._mirrors = accum.mirrors;    
    }
  }; // end of distributed_hybrid_ginger_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
