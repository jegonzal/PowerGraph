/*  
 * Copyright (c) 2013 Shanghai Jiao Tong University. 
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
 *
 * 2013.11  implement hybrid random partitioning
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_HYBRID_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_HYBRID_INGRESS_HPP

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

#include <graphlab/macros_def.hpp>

#define TUNING
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges using a hybrid method.
   *        That is, for high degree edge, hashing from its source vertex;
   *        for low degree edge, hashing from its target vertex.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_hybrid_ingress : 
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

    /// The rpc interface for this object
    dc_dist_object<distributed_hybrid_ingress> hybrid_rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;
    
    /// threshold to divide high-degree and low-degree vertices
    size_t threshold;

    bool standalone;

    typedef typename base_type::edge_buffer_record edge_buffer_record;
    typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
        edge_buffer_type;

    typedef typename base_type::vertex_buffer_record vertex_buffer_record;
    typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;

    std::vector<edge_buffer_record> hybrid_edges;

    /* ingress exchange */
    buffered_exchange<edge_buffer_record> hybrid_edge_exchange;
    buffered_exchange<vertex_buffer_record> hybrid_vertex_exchange;

    /// detail vertex record for the second pass coordination. 
    typedef typename base_type::vertex_negotiator_record 
      vertex_negotiator_record;

  public:
    distributed_hybrid_ingress(distributed_control& dc, 
        graph_type& graph, size_t threshold = 100) :
        base_type(dc, graph), hybrid_rpc(dc, this), 
        graph(graph), threshold(threshold),
#ifdef _OPENMP
        hybrid_edge_exchange(dc, omp_get_max_threads()), 
        hybrid_vertex_exchange(dc, omp_get_max_threads())
#else
        hybrid_edge_exchange(dc), 
        hybrid_vertex_exchange(dc)
#endif
    {
      /* fast pass for standalone case. */
      standalone = hybrid_rpc.numprocs() == 1;
      hybrid_rpc.barrier();
    } // end of constructor

    ~distributed_hybrid_ingress() { }

    /** Add an edge to the ingress object using random hashing assignment.
     *  This function acts as the first phase for SNAP graph to deliver edges
     *  via the hashing value of its target vertex.
     */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      const edge_buffer_record record(source, target, edata);      
      const procid_t owning_proc = standalone ? 0 :
        graph_hash::hash_vertex(target) % hybrid_rpc.numprocs();
#ifdef _OPENMP
      hybrid_edge_exchange.send(owning_proc, record, omp_get_thread_num());
#else
      hybrid_edge_exchange.send(owning_proc, record);
#endif
    } // end of add edge


    /* add vdata */
    void add_vertex(vertex_id_type vid, const VertexData& vdata) { 
      const vertex_buffer_record record(vid, vdata);
      const procid_t owning_proc = standalone ? 0 :
        graph_hash::hash_vertex(vid) % hybrid_rpc.numprocs();        
#ifdef _OPENMP
      hybrid_vertex_exchange.send(owning_proc, record, omp_get_thread_num());
#else
      hybrid_vertex_exchange.send(owning_proc, record);
#endif
    } // end of add vertex


    void finalize() {
      
      graphlab::timer ti;
      
      size_t nprocs = hybrid_rpc.numprocs();
      procid_t l_procid = hybrid_rpc.procid();
      size_t nedges = 0;

      hybrid_rpc.full_barrier();

      if (l_procid == 0) {
        memory_info::log_usage("start finalizing");
        logstream(LOG_EMPH) << "hybrid finalizing ..."
                            << " #vertices=" << graph.local_graph.num_vertices()
                            << " #edges=" << graph.local_graph.num_edges()
                            << " threshold=" << threshold
                            << std::endl;
      }


      /**************************************************************************/
      /*                                                                        */
      /*                       Flush any additional data                        */
      /*                                                                        */
      /**************************************************************************/
      hybrid_edge_exchange.flush(); hybrid_vertex_exchange.flush();

      /**
       * Fast pass for redundant finalization with no graph changes. 
       */
      {
        size_t changed_size = hybrid_edge_exchange.size() + hybrid_vertex_exchange.size();
        hybrid_rpc.all_reduce(changed_size);
        if (changed_size == 0) {
          logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
          return;
        }
      }
      

      /**************************************************************************/
      /*                                                                        */
      /*                       Prepare hybrid ingress                           */
      /*                                                                        */
      /**************************************************************************/ 
      {
        edge_buffer_type edge_buffer;
        procid_t proc;
        nedges = hybrid_edge_exchange.size();

        hybrid_edges.reserve(nedges);
        if (standalone) { /* fast pass for standalone */
          proc = -1;
          while(hybrid_edge_exchange.recv(proc, edge_buffer))
            foreach(const edge_buffer_record& rec, edge_buffer)
              hybrid_edges.push_back(rec);          
          hybrid_edge_exchange.clear();
        } else {
          hopscotch_map<vertex_id_type, size_t> in_degree_set;

          proc = -1;
          while(hybrid_edge_exchange.recv(proc, edge_buffer)) {
            foreach(const edge_buffer_record& rec, edge_buffer) {
              hybrid_edges.push_back(rec);
              in_degree_set[rec.target]++;
            }
          }
          hybrid_edge_exchange.clear();
          hybrid_edge_exchange.barrier(); // barrier before reusing
#ifdef TUNING
          if(l_procid == 0) {
            logstream(LOG_INFO) << "save local edges and count in-degree: " 
                                << ti.current_time()
                                << " secs" 
                                << std::endl;
          }
#endif

          // re-send edges of high-degree vertices
          for (size_t i = 0; i < hybrid_edges.size(); i++) {
            edge_buffer_record& rec = hybrid_edges[i];
            if (in_degree_set[rec.target] > threshold) {
              const procid_t source_owner_proc = 
                graph_hash::hash_vertex(rec.source) % nprocs;
              if(source_owner_proc != l_procid){
                // re-send the edge of high-degree vertices according to source
                hybrid_edge_exchange.send(source_owner_proc, rec);
                // set re-sent edges as empty for skipping
                hybrid_edges[i] = edge_buffer_record();
                --nedges;
              }
            }
          }
#ifdef TUNING
          if(l_procid == 0) {
            logstream(LOG_INFO) << "resend edges of high-degree vertices: " 
                                << ti.current_time()
                                << " secs" 
                                << std::endl;
          }
#endif

          // receive edges of high-degree vertices
          hybrid_edge_exchange.flush();
#ifdef TUNING
          logstream(LOG_INFO) << "receive high-degree edges: "  
                              << hybrid_edge_exchange.size() << std::endl;
#endif
          proc = -1;
          while(hybrid_edge_exchange.recv(proc, edge_buffer)) {
            foreach(const edge_buffer_record& rec, edge_buffer) {
              hybrid_edges.push_back(rec);
              ++nedges;
            }
          }
          hybrid_edge_exchange.clear();
          in_degree_set.clear();
#ifdef TUNING
          if(l_procid == 0) {
            logstream(LOG_INFO) << "receive high-degree edges: " 
                                << ti.current_time()
                                << " secs" 
                                << std::endl;
          }
#endif
        }
      }

      if(l_procid == 0) {
        memory_info::log_usage("prepare hybrid finalizing done.");
        logstream(LOG_EMPH) << "prepare hybrid finalizing. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }
      
      // connect to base finalize()
      modified_base_finalize(nedges);
      
      // set vertex degree type for hybrid engine
      set_degree_type();

      if(l_procid == 0) {
        memory_info::log_usage("hybrid finalizing graph done.");
        logstream(LOG_EMPH) << "hybrid finalizing graph. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }
    } // end of finalize

    void set_degree_type() {
      graphlab::timer ti;
      procid_t l_procid = hybrid_rpc.procid();
      size_t high_master = 0, high_mirror = 0, low_master = 0, low_mirror = 0;
      
      for (size_t lvid = 0; lvid < graph.num_local_vertices(); lvid++) {
        vertex_record& vrec = graph.lvid2record[lvid];
        if (vrec.num_in_edges > threshold) {
          vrec.dtype = graph_type::HIGH; 
          if (vrec.owner == l_procid) high_master ++;
          else high_mirror ++;
        } else {
          vrec.dtype = graph_type::LOW; 
          if (vrec.owner == l_procid) low_master ++;
          else low_mirror ++;
        }        
      }

#ifdef TUNING
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
                            << (float(high_master)/(high_master+low_master)) << "]"
                            << std::endl;
        if ((high_mirror + low_mirror) > 0)
        logstream(LOG_EMPH) << "hybrid info: mirror [" 
                            << high_mirror << " " 
                            << low_mirror << " " 
                            << (float(high_mirror)/(high_mirror+low_mirror)) << "]"
                            << std::endl;

        memory_info::log_usage("set vertex type done."); 
        logstream(LOG_EMPH) << "set vertex type: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif
    }


    /* 
     * do the same job as original base finalize except for
     * extracting edges from hybrid_edges instead of original edge_buffer
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
      { // Add all the edges to the local graph
        graph.local_graph.reserve_edge_space(nedges + 1);

        foreach(const edge_buffer_record& rec, hybrid_edges) {
          // skip re-sent edges
          if (rec.source == vertex_id_type(-1)) continue;

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
          logstream(LOG_INFO) << "populating local graph: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif
        // Finalize local graph
        graph.local_graph.finalize();
#ifdef TUNING
        logstream(LOG_INFO) << "local graph info: " << std::endl
                            << "\t nverts: " << graph.local_graph.num_vertices()
                            << std::endl
                            << "\t nedges: " << graph.local_graph.num_edges()
                            << std::endl;
        
        if(l_procid == 0) {
          logstream(LOG_INFO) << "finalizing local graph: " 
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
              if (distributed_hybrid_ingress::vertex_combine_strategy 
                && lvid < graph.num_local_vertices()) {
                distributed_hybrid_ingress::vertex_combine_strategy(
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
            logstream(LOG_INFO) << "adding vertex data: " 
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
          if (standalone)
            vrec.owner = 0;
          else
            vrec.owner = graph_hash::hash_vertex(pair.first) % nprocs;
        }
        ASSERT_EQ(local_nverts, graph.local_graph.num_vertices());
        ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
#ifdef TUNING
        if(l_procid == 0) {
          logstream(LOG_INFO) << "allocating lvid2record: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif
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
          logstream(LOG_INFO) << "#flying-own-nverts="
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
        logstream(LOG_INFO) << "master handshake: " 
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
                             boost::bind(&distributed_hybrid_ingress::finalize_gather, this, _1, _2), 
                             boost::bind(&distributed_hybrid_ingress::finalize_apply, this, _1, _2, _3));
        vrecord_sync_gas.exec(changed_vset);

#ifdef TUNING
        if(l_procid == 0) {
          logstream(LOG_INFO) << "synchrionizing vertex (meta)data: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
#endif
      }

      base_type::exchange_global_info(standalone);
#ifdef TUNING
      if(l_procid == 0) {
        logstream(LOG_INFO) << "exchange global info: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif

      if(l_procid == 0) {
        memory_info::log_usage("base finalizing done.");
        logstream(LOG_EMPH) << "base finalizing. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }
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
  }; // end of distributed_hybrid_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
