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

#ifndef GRAPHLAB_DISTRIBUTED_INGRESS_BASE_HPP
#define GRAPHLAB_DISTRIBUTED_INGRESS_BASE_HPP

#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/ingress_edge_decision.hpp>
#include <graphlab/graph/ingress/vertex_channel.hpp>
#include <graphlab/parallel/mutex.hpp>
#include <graphlab/util/memory_info.hpp>
#include <graphlab/util/cuckoo_map_pow2.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \brief Implementation of the basic ingress functionality.
   */
  template<typename VertexData, typename EdgeData>
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
    typedef vertex_channel<VertexData, EdgeData> vertex_channel_type;

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

    typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
        edge_buffer_type;
    typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;

    typedef boost::unordered_map<vertex_id_type, size_t> vchannel_map_type;

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


    /**
     * Defines the strategy to use when duplicate vertices are inserted.
     * The default behavior is that an arbitrary vertex data is picked.
     * This allows you to define a combining strategy.
     */
    void set_duplicate_vertex_strategy(boost::function<void(vertex_data_type&,
                                                        const vertex_data_type&)>
                                       combine_strategy) {
      vertex_combine_strategy = combine_strategy;
    }

    /** \brief Finalize completes the local graph data structure 
     * and the vertex record information. 
     *
     * \internal
     * The finalization goes through 5 steps:
     *
     * 1. Receive edges and construct local graph, during which
     * the part of the vid2lvid map and vchannel_map is built.
     *
     * 2. Receive vertices and fillin the data entry in vchannel_map.
     *
     * 3. Gather "singleton" vertices, and finalize the structure of vid2lvid, vchannel, lvid2record.
     *
     * 4. Exchange the vertex information by gather/scatter of vchannel_map.
     *
     * 5. Exchange global graph statistics.
     */
    virtual void finalize() {
      timer mytimer;

      rpc.full_barrier();
      if (rpc.procid() == 0) {
        logstream(LOG_EMPH) << "Finalizing Graph..." << std::endl;
      }

      // Flush any additional data
      edge_exchange.flush(); vertex_exchange.flush();     
      if(rpc.procid() == 0)       
        memory_info::log_usage("Post Flush " + 
                               boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");


      ////////////////////////////// Phase 1 ////////////////////////////
      logstream(LOG_INFO) << "Graph Finalize: constructing local graph" << std::endl;
      { // Add all the edges to the local graph
        const size_t nedges = edge_exchange.size()+1;
        graph.local_graph.reserve_edge_space(nedges + 1);      
        edge_buffer_type edge_buffer;
        procid_t proc;
        while(edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            // Get the source_vid, target_vid;
            lvid_type source_lvid = get_or_add_vid2lvid(rec.source);
            lvid_type target_lvid = get_or_add_vid2lvid(rec.target);

            // Add edge to the local graph
            graph.local_graph.add_edge(source_lvid, target_lvid, rec.edata);

            // Add vertices to vertex channel for synchronization 
            add_vertex_channel(rec.source);
            add_vertex_channel(rec.target);

          } // end of loop over add edges
        } // end for loop over buffers
        edge_exchange.clear();
      }
      if(rpc.procid() == 0) 
        memory_info::log_usage("Finish populating local graph " + 
                               boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");
      
      ASSERT_EQ(graph.vid2lvid.size(), graph.local_graph.num_vertices());
      
      // Finalize local graph
      logstream(LOG_INFO) << "Graph Finalize: finalizing local graph." 
                          << std::endl;
      graph.local_graph.finalize();
      logstream(LOG_INFO) << "Local graph info: " << std::endl
                          << "\t nverts: " << graph.local_graph.num_vertices()
                          << std::endl
                          << "\t nedges: " << graph.local_graph.num_edges()
                          << std::endl;
      if(rpc.procid() == 0)       
        memory_info::log_usage("Finish finalizing local graph " + 
                               boost::lexical_cast<std::string>(mytimer.current_time()) + " secs"); 

      ////////////////////////////// Phase 2 ////////////////////////////
      // Receive vertex data, setup the map containing all the updated vertices  
      { 
        vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
        while(vertex_exchange.recv(sending_proc, vertex_buffer)) {
          foreach(const vertex_buffer_record& rec, vertex_buffer) {
            get_or_add_vid2lvid(rec.vid);
            add_vertex_channel(rec.vid);
            vertex_channel_type& vchannel = vchannel_list[vchannel_map[rec.vid]];
            if (vertex_combine_strategy &&  vchannel.data_is_set) {
              vertex_combine_strategy(vchannel.vdata, rec.vdata);
            } else {
              vchannel.vdata = rec.vdata;
              vchannel.data_is_set = true;
            }
          }
        }
        vertex_exchange.clear();
      } // end of loop to populate vrecmap

      rpc.barrier();
      if(rpc.procid() == 0)         
        memory_info::log_usage("Finish receving vertex data " + 
                               boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");

      /////////////////////////// Phase 3 /////////////////////////
      { // synchornize all updated vertices between master and mirrors
        // mirror vertex notify the owner proc
        // collect vids to owners 
        std::vector<std::vector<vertex_id_type> > vid_exchange(rpc.numprocs());
        std::vector<mutex> vid_exchange_locks(rpc.numprocs());
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (ssize_t i = 0; i < vchannel_list.size(); ++i) {
          vertex_channel_type& vc = vchannel_list[i]; 
          procid_t owner = vc.owner;
          vid_exchange_locks[owner].lock();
          vid_exchange[vc.owner].push_back(vc.vid);
          vid_exchange_locks[owner].unlock();
        }

        rpc.all_to_all(vid_exchange);

        // check any new vids
        mutex mtx;
        std::set<vertex_id_type> new_vids;
        for (size_t i = 0; i < rpc.numprocs(); ++i) {
          if (i == rpc.procid()) 
            continue;
          std::vector<vertex_id_type>& vid_vec = vid_exchange[i];
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for (ssize_t j = 0; j < vid_vec.size(); ++j) {
            if (vchannel_map.find(vid_vec[j]) == vchannel_map.end()) {
              mtx.lock();
              new_vids.insert(vid_vec[j]);
              mtx.unlock();
            }
          }
          std::vector<vertex_id_type>().swap(vid_vec);
        }
        foreach(vertex_id_type vid, new_vids) {
          add_vertex_channel(vid);
          get_or_add_vid2lvid(vid);
        }
        
        rpc.barrier();
        if(rpc.procid() == 0)         
          memory_info::log_usage("Finish master vid exchange " + 
                                 boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");

#ifdef _OPENMP
#pragma omp parallel for
#endif
        // owner vertex gather from mirrors 
        for (ssize_t i = 0; i < vchannel_list.size(); ++i) {
          vertex_channel_type& vc = vchannel_list[i];
          vertex_id_type vid = vc.vid;
          lvid_type lvid = graph.vid2lvid[vid];
          vc.num_in_edges = graph.local_graph.num_in_edges(lvid);
          vc.num_out_edges = graph.local_graph.num_out_edges(lvid);
          if (vc.owner != rpc.procid())
            rpc.remote_call(vc.owner,
                            &distributed_ingress_base::vchannel_gather, vc);
        }
        rpc.full_barrier();
        if (rpc.procid() == 0)
          memory_info::log_usage("Finish vertex gather " + 
                                 boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");

#ifdef _OPENMP
#pragma omp parallel for
#endif
       for (ssize_t i = 0; i < vchannel_list.size(); ++i) {
          vertex_channel_type& vc = vchannel_list[i];
          if (vc.owner == rpc.procid()) {
            //foreach mirror
            typename mirror_type::iterator it = vc.mirrors.begin();
            while (it != vc.mirrors.end()) {
              rpc.remote_call(*it,
                              &distributed_ingress_base::vchannel_scatter, vc);
              ++it;
            }
          }
        }

        if (rpc.procid() == 0)
          memory_info::log_usage("Finish vertex scatter " + 
                                 boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");

        // reserve and resize vertices datastructre in the graph 
        size_t vsize = graph.vid2lvid.size();
        if (vsize > graph.lvid2record.size()) {
          graph.lvid2record.reserve(vsize);
          graph.lvid2record.resize(vsize);
        }
        if (vsize > graph.local_graph.num_vertices()) {
          graph.local_graph.reserve(vsize);
          graph.local_graph.resize(vsize);
        }

        rpc.full_barrier();

        // update local graph vertex record with updated vertices info
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (ssize_t i = 0; i < vchannel_list.size(); ++i) {
          const vertex_channel_type& vc = vchannel_list[i];
          vertex_id_type gvid = vc.vid; 
          ASSERT_TRUE(graph.vid2lvid.find(gvid) != graph.vid2lvid.end());
          lvid_type lvid = graph.vid2lvid[gvid];
          ASSERT_LT(lvid, graph.lvid2record.size());
          vertex_record& vrecord = graph.lvid2record[lvid];
          vrecord.gvid = gvid;      
          vrecord._mirrors = vc.mirrors;
          vrecord.owner = vc.owner;
          vrecord.num_in_edges = vc.num_in_edges;
          vrecord.num_out_edges = vc.num_out_edges;
          if (vc.data_is_set) {
            graph.local_graph.add_vertex(lvid, vc.vdata);
          }
        }
      } // end of phase 3
      rpc.full_barrier();
      if(rpc.procid() == 0) 
        memory_info::log_usage("Finish updating vertices " + 
                               boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");

      ASSERT_EQ(graph.vid2lvid.size(), graph.local_graph.num_vertices());
      ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
 
      exchange_global_info();
      clear();
      rpc.full_barrier();

      if(rpc.procid() == 0)  {
        memory_info::log_usage("Complete finalizing graph " + 
                               boost::lexical_cast<std::string>(mytimer.current_time()) + " secs");
        logstream(LOG_EMPH) << "Graph info: "  
                            << "\n\t nverts: " << graph.num_vertices()
                            << "\n\t nedges: " << graph.num_edges()
                            << "\n\t nreplicas: " << graph.nreplicas
                            << "\n\t replication factor: " << (double)graph.nreplicas/graph.num_vertices()
                            << std::endl;
      }
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

      // compute begin edge id
      graph.begin_eid = 0;
      for(size_t i = 0; i < rpc.procid(); ++i) 
        graph.begin_eid += swap_counts[i];

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
    }


    ///////// Helper Functions /////////////
  private:
    lvid_type get_or_add_vid2lvid(vertex_id_type vid) {
      lvid_type lvid(-1);
      if(graph.vid2lvid.find(vid) == graph.vid2lvid.end()) {
        lvid = graph.vid2lvid.size();
        graph.vid2lvid[vid] = lvid;
      } else {
        lvid = graph.vid2lvid[vid];
      }
      return lvid;
    }

    void add_vertex_channel(vertex_id_type vid) {
      if (vchannel_map.find(vid) == vchannel_map.end()) {
        vchannel_map[vid] = vchannel_list.size();
        vertex_channel_type vc;
        vc.owner = graph_hash::hash_vertex(vid) % rpc.numprocs();
        vc.vid = vid;
        if (vc.owner != rpc.procid()) {
          vc.mirrors.set_bit(rpc.procid());
        }
        vchannel_list.push_back(vc);
      }
    }

    void vchannel_gather(const vertex_channel_type& acc) {
      ASSERT_TRUE(vchannel_map.find(acc.vid) != vchannel_map.end());
      vertex_channel_type& vc = vchannel_list[vchannel_map[acc.vid]];
      ASSERT_TRUE(vc.owner == rpc.procid());
      vc.mtx.lock();
      vc += acc;
      vc.mtx.unlock();
    }

    void vchannel_scatter(const vertex_channel_type& acc) {
      ASSERT_TRUE(vchannel_map.find(acc.vid) != vchannel_map.end());
      vertex_channel_type& vc = vchannel_list[vchannel_map[acc.vid]];
      ASSERT_TRUE(vc.owner != rpc.procid());
      vc = acc;
    }

    void clear() { 
      vchannel_map.clear();
      std::vector<vertex_channel_type>().swap(vchannel_list);
      vertex_exchange.clear();
      edge_exchange.clear();
    }

    ///////// Helper Datastructure //////////
  protected:
    /// The rpc interface for this object
    dc_dist_object<distributed_ingress_base> rpc;

    /// The underlying distributed graph object that is being loaded
    graph_type& graph;

    buffered_exchange<vertex_buffer_record> vertex_exchange;

    buffered_exchange<edge_buffer_record> edge_exchange;

    /// Ingress decision object for computing the edge destination. 
    ingress_edge_decision<VertexData, EdgeData> edge_decision;

    vchannel_map_type vchannel_map;
    std::vector<vertex_channel_type> vchannel_list;

    boost::function<void(vertex_data_type&,
                         const vertex_data_type&)> vertex_combine_strategy;

  }; // end of distributed_ingress_base
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
