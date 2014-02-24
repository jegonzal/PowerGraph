/**  
 * Copyright (c) 2009 Carnegie Mellon University.  All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not
 *  use this file except in compliance with the License.
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

#ifndef GRAPHLAB_DISTRIBUTED_INGRESS_BASE_HPP
#define GRAPHLAB_DISTRIBUTED_INGRESS_BASE_HPP

#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/graph_hash.hpp>
#include <graphlab/graph/ingress/ingress_edge_decision.hpp>
#include <graphlab/graph/graph_gather_apply.hpp>
#include <graphlab/util/memory_info.hpp>
#include <graphlab/util/hopscotch_map.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \brief Implementation of the basic ingress functionality.
   */
  template <typename VertexType, typename EdgeType> 
  class distributed_graph;

  template<typename VertexData, typename EdgeData>
  class distributed_ingress_base {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;
    
    typedef typename graph_type::vertex_record vertex_record;
    typedef typename graph_type::mirror_type mirror_type;

   
    /// The rpc interface for this object
    dc_dist_object<distributed_ingress_base> rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;

    /// Temporary buffers used to store vertex data on ingress
    struct vertex_buffer_record {
      vertex_id_type vid;
      vertex_data_type vdata;
      vertex_buffer_record(vertex_id_type vid = -1,
                           vertex_data_type vdata = vertex_data_type()) :
        vid(vid), vdata(vdata) { }
      void load(iarchive& arc) { arc >> vid >> vdata; }
      void save(oarchive& arc) const { arc << vid << vdata; }
    }; 
    buffered_exchange<vertex_buffer_record> vertex_exchange;

    /// Temporar buffers used to store edge data on ingress
    struct edge_buffer_record {
      vertex_id_type source, target;
      edge_data_type edata;
      edge_buffer_record(const vertex_id_type& source = vertex_id_type(-1), 
                         const vertex_id_type& target = vertex_id_type(-1), 
                         const edge_data_type& edata = edge_data_type()) :
        source(source), target(target), edata(edata) { }
      void load(iarchive& arc) { arc >> source >> target >> edata; }
      void save(oarchive& arc) const { arc << source << target << edata; }
    };
    buffered_exchange<edge_buffer_record> edge_exchange;

    /// Detail vertex record for the second pass coordination. 
    struct vertex_negotiator_record {
      mirror_type mirrors;
      vertex_id_type num_in_edges, num_out_edges;
      bool has_data;
      vertex_data_type vdata;
      vertex_negotiator_record() : num_in_edges(0), num_out_edges(0), has_data(false) { }

      void load(iarchive& arc) { 
        arc >> num_in_edges >> num_out_edges >> mirrors >> has_data >> vdata;
      }
      void save(oarchive& arc) const { 
        arc << num_in_edges << num_out_edges << mirrors << has_data << vdata;
      }

      vertex_negotiator_record operator+=(const vertex_negotiator_record& v2) {
        num_in_edges += v2.num_in_edges;
        num_out_edges += v2.num_out_edges;
        mirrors |= v2.mirrors;
        if (v2.has_data) {
          vdata = v2.vdata;
        }
        return *this;
      }
    };

    /// Ingress decision object for computing the edge destination. 
    ingress_edge_decision<VertexData, EdgeData> edge_decision;

  public:
    distributed_ingress_base(distributed_control& dc, graph_type& graph) :
      rpc(dc, this), graph(graph), vertex_exchange(dc), edge_exchange(dc),
      edge_decision(dc) {
      rpc.barrier();
    } // end of constructor

    virtual ~distributed_ingress_base() { }

    /** \brief Add an edge to the ingress object. */
    virtual void add_edge(vertex_id_type source, vertex_id_type target,
                          const EdgeData& edata) {
      const procid_t owning_proc = 
        edge_decision.edge_to_proc_random(source, target, rpc.numprocs());
      const edge_buffer_record record(source, target, edata);
      edge_exchange.send(owning_proc, record);
    } // end of add edge


    /** \brief Add an vertex to the ingress object. */
    virtual void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      const procid_t owning_proc = graph_hash::hash_vertex(vid) % rpc.numprocs();
      const vertex_buffer_record record(vid, vdata);
      vertex_exchange.send(owning_proc, record);
    } // end of add vertex


    void set_duplicate_vertex_strategy(
        boost::function<void(vertex_data_type&,
                             const vertex_data_type&)> combine_strategy) {
      vertex_combine_strategy = combine_strategy;
    }



    /** \brief Finalize completes the local graph data structure 
     * and the vertex record information. 
     *
     * \internal
     * The finalization goes through 5 steps:
     *
     * 1. Construct local graph using the received edges, during which
     * the vid2lvid map is built.
     *
     * 2. Construct lvid2record map (of empty entries) using the received vertices. 
     *
     * 3. Complete lvid2record map by exchanging the vertex_info. 
     *
     * 4. Exchange the negotiation records, including singletons. (Local graph 
     * handling singletons). 
     *
     * 5. Exchange global graph statistics.
     */
    virtual void finalize() {

      rpc.full_barrier();

      bool first_time_finalize = false;
      /**
       * Fast pass for first time finalization. 
       */
      if (graph.is_dynamic()) {
        size_t nverts = graph.num_local_vertices();
        rpc.all_reduce(nverts);
        first_time_finalize = (nverts == 0);
      } else {
        first_time_finalize = false;
      }


      if (rpc.procid() == 0) {
        logstream(LOG_EMPH) << "Finalizing Graph..." << std::endl;
      }

      typedef typename hopscotch_map<vertex_id_type, lvid_type>::value_type
        vid2lvid_pair_type;

      typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
        edge_buffer_type;

      typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;

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
      /*                       Flush any additional data                        */
      /*                                                                        */
      /**************************************************************************/
      edge_exchange.flush(); vertex_exchange.flush();     

      /**
       * Fast pass for redundant finalization with no graph changes. 
       */
      {
        size_t changed_size = edge_exchange.size() + vertex_exchange.size();
        rpc.all_reduce(changed_size);
        if (changed_size == 0) {
          logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
          return;
        }
      }

      if(rpc.procid() == 0)       
        memory_info::log_usage("Post Flush");

     
      /**************************************************************************/
      /*                                                                        */
      /*                         Construct local graph                          */
      /*                                                                        */
      /**************************************************************************/
      { // Add all the edges to the local graph
        logstream(LOG_INFO) << "Graph Finalize: constructing local graph" << std::endl;
        const size_t nedges = edge_exchange.size()+1;
        graph.local_graph.reserve_edge_space(nedges + 1);      
        edge_buffer_type edge_buffer;
        procid_t proc;
        while(edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
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
            // std::cout << "add edge " << rec.source << "\t" << rec.target << std::endl;
          } // end of loop over add edges
        } // end for loop over buffers
        edge_exchange.clear();

        ASSERT_EQ(graph.vid2lvid.size()  + vid2lvid_buffer.size(), graph.local_graph.num_vertices());
        if(rpc.procid() == 0)  {
          memory_info::log_usage("Finished populating local graph.");
        }

        // Finalize local graph
        logstream(LOG_INFO) << "Graph Finalize: finalizing local graph." 
                            << std::endl;
        graph.local_graph.finalize();
        logstream(LOG_INFO) << "Local graph info: " << std::endl
                            << "\t nverts: " << graph.local_graph.num_vertices()
                            << std::endl
                            << "\t nedges: " << graph.local_graph.num_edges()
                            << std::endl;
        
        if(rpc.procid() == 0) {
          memory_info::log_usage("Finished finalizing local graph."); 
          // debug
          // std::cout << graph.local_graph << std::endl;
        }
      }

      /**************************************************************************/
      /*                                                                        */
      /*             Receive and add vertex data to masters                     */
      /*                                                                        */
      /**************************************************************************/
      // Setup the map containing all the vertices being negotiated by this machine
      { // Receive any vertex data sent by other machines
        vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
        while(vertex_exchange.recv(sending_proc, vertex_buffer)) {
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
            if (vertex_combine_strategy && lvid < graph.num_local_vertices()) {
              vertex_combine_strategy(graph.l_vertex(lvid).data(), rec.vdata);
            } else {
              graph.local_graph.add_vertex(lvid, rec.vdata);
            }
          }
        }
        vertex_exchange.clear();
        if(rpc.procid() == 0)         
          memory_info::log_usage("Finished adding vertex data");
      } // end of loop to populate vrecmap



      /**************************************************************************/
      /*                                                                        */
      /*        assign vertex data and allocate vertex (meta)data  space        */
      /*                                                                        */
      /**************************************************************************/
      { // Determine masters for all negotiated vertices
        const size_t local_nverts = graph.vid2lvid.size() + vid2lvid_buffer.size();
        graph.lvid2record.reserve(local_nverts);
        graph.lvid2record.resize(local_nverts);
        graph.local_graph.resize(local_nverts);
        foreach(const vid2lvid_pair_type& pair, vid2lvid_buffer) {
            vertex_record& vrec = graph.lvid2record[pair.second];
            vrec.gvid = pair.first;
            vrec.owner = graph_hash::hash_vertex(pair.first) % rpc.numprocs();
        }
        ASSERT_EQ(local_nverts, graph.local_graph.num_vertices());
        ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
        if(rpc.procid() == 0)       
          memory_info::log_usage("Finihsed allocating lvid2record");
      }

      /**************************************************************************/
      /*                                                                        */
      /*                          Master handshake                              */
      /*                                                                        */
      /**************************************************************************/
      {
#ifdef _OPENMP
        buffered_exchange<vertex_id_type> vid_buffer(rpc.dc(), omp_get_max_threads());
#else
        buffered_exchange<vertex_id_type> vid_buffer(rpc.dc());
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
        // send not owned vids to their master
        for (lvid_type i = lvid_start; i < graph.lvid2record.size(); ++i) {
          procid_t master = graph.lvid2record[i].owner;
          if (master != rpc.procid())
#ifdef _OPENMP
            vid_buffer.send(master, graph.lvid2record[i].gvid, omp_get_thread_num());
#else
            vid_buffer.send(master, graph.lvid2record[i].gvid);
#endif
        }
        vid_buffer.flush();
        rpc.barrier();

        // receive all vids owned by me
        mutex flying_vids_lock;
        boost::unordered_map<vertex_id_type, mirror_type> flying_vids;
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          typename buffered_exchange<vertex_id_type>::buffer_type buffer;
          procid_t recvid;
          while(vid_buffer.recv(recvid, buffer)) {
            foreach(const vertex_id_type vid, buffer) {
              if (graph.vid2lvid.find(vid) == graph.vid2lvid.end()) {
                if (vid2lvid_buffer.find(vid) == vid2lvid_buffer.end()) {
                  flying_vids_lock.lock();
                  mirror_type& mirrors = flying_vids[vid];
                  flying_vids_lock.unlock();
                  mirrors.set_bit(recvid);
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
        // reallocate spaces for the flying vertices. 
        size_t vsize_old = graph.lvid2record.size();
        size_t vsize_new = vsize_old + flying_vids.size();
        graph.lvid2record.resize(vsize_new);
        graph.local_graph.resize(vsize_new);
        for (typename boost::unordered_map<vertex_id_type, mirror_type>::iterator it = flying_vids.begin();
             it != flying_vids.end(); ++it) {
          lvid_type lvid = lvid_start + vid2lvid_buffer.size();
          vertex_id_type gvid = it->first; 
          graph.lvid2record[lvid].owner = rpc.procid();
          graph.lvid2record[lvid].gvid = gvid;
          graph.lvid2record[lvid]._mirrors= it->second;
          vid2lvid_buffer[gvid] = lvid;
          // std::cout << "proc " << rpc.procid() << " recevies flying vertex " << gvid << std::endl;
        }
      } // end of master handshake

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
          // vid2lvid_buffer.swap(vid2lvid_map_type(-1));
        }
      }


      /**************************************************************************/
      /*                                                                        */
      /*              synchronize vertex data and meta information              */
      /*                                                                        */
      /**************************************************************************/
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
          buffered_exchange<vertex_id_type> vset_exchange(rpc.dc());
          // sync vset with all mirrors
          changed_vset.synchronize_mirrors_to_master_or(graph, vset_exchange);
          changed_vset.synchronize_master_to_mirrors(graph, vset_exchange);
        }

        graphlab::graph_gather_apply<graph_type, vertex_negotiator_record> 
            vrecord_sync_gas(graph, 
                             boost::bind(&distributed_ingress_base::finalize_gather, this, _1, _2), 
                             boost::bind(&distributed_ingress_base::finalize_apply, this, _1, _2, _3));
        vrecord_sync_gas.exec(changed_vset);

        if(rpc.procid() == 0)       
          memory_info::log_usage("Finished synchronizing vertex (meta)data");
      }

      exchange_global_info();
    } // end of finalize


    /* Exchange graph statistics among all nodes and compute
     * global statistics for the distributed graph. */
    void exchange_global_info () {
      // Count the number of vertices owned locally
      graph.local_own_nverts = 0;
      foreach(const vertex_record& record, graph.lvid2record)
        if(record.owner == rpc.procid()) ++graph.local_own_nverts;

      // Finalize global graph statistics. 
      logstream(LOG_INFO)
        << "Graph Finalize: exchange global statistics " << std::endl;

      // Compute edge counts
      std::vector<size_t> swap_counts(rpc.numprocs());
      swap_counts[rpc.procid()] = graph.num_local_edges();
      rpc.all_gather(swap_counts);
      graph.nedges = 0;
      foreach(size_t count, swap_counts) graph.nedges += count;


      // compute vertex count
      swap_counts[rpc.procid()] = graph.num_local_own_vertices();
      rpc.all_gather(swap_counts);
      graph.nverts = 0;
      foreach(size_t count, swap_counts) graph.nverts += count;

      // compute replicas
      swap_counts[rpc.procid()] = graph.num_local_vertices();
      rpc.all_gather(swap_counts);
      graph.nreplicas = 0;
      foreach(size_t count, swap_counts) graph.nreplicas += count;


      if (rpc.procid() == 0) {
        logstream(LOG_EMPH) << "Graph info: "  
                            << "\n\t nverts: " << graph.num_vertices()
                            << "\n\t nedges: " << graph.num_edges()
                            << "\n\t nreplicas: " << graph.nreplicas
                            << "\n\t replication factor: " << (double)graph.nreplicas/graph.num_vertices()
                            << std::endl;
      }
    }


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
     * \brief Update the vertex datastructures with the gathered vertex metadata.  
     */
    void finalize_apply(lvid_type lvid, const vertex_negotiator_record& accum, graph_type& graph) {
        typename graph_type::vertex_record& vrec = graph.lvid2record[lvid];
        vrec.num_in_edges = accum.num_in_edges;
        vrec.num_out_edges = accum.num_out_edges;
        graph.l_vertex(lvid).data() = accum.vdata;
        vrec._mirrors = accum.mirrors;
    }
  }; // end of distributed_ingress_base
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
