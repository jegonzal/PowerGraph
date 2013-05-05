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

    boost::function<void(vertex_data_type&,
                         const vertex_data_type&)> vertex_combine_strategy;
    
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

   
    /// Struct of minimal vertex information for the first pass coordination. 
    struct vertex_info : public graphlab::IS_POD_TYPE {
      vertex_id_type vid, num_in_edges, num_out_edges;
      vertex_info(vertex_id_type vid = 0, vertex_id_type num_in_edges = 0,
                  vertex_id_type num_out_edges = 0) : 
        vid(vid), num_in_edges(num_in_edges), num_out_edges(num_out_edges) { }     
    }; // end of vertex_info


    /// Detail vertex record for the second pass coordination. 
    struct vertex_negotiator_record {
      mirror_type mirrors;
      vertex_data_type vdata;
      vertex_id_type num_in_edges, num_out_edges;
      procid_t owner;
      bool data_is_set;
      vertex_negotiator_record() : 
        vdata(vertex_data_type()), num_in_edges(0), num_out_edges(0), 
        owner(-1), data_is_set(false) { }
      void load(iarchive& arc) { 
        arc >> num_in_edges >> num_out_edges >> owner >> mirrors >> vdata;
      }
      void save(oarchive& arc) const { 
        arc << num_in_edges << num_out_edges << owner << mirrors << vdata;
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
     * 1. Construct local graph using the received edges, during which
     * the vid2lvid map is built.
     *
     * 2. Construct lvid2record map (of empty entries) using the received vertices. 
     *
     * 3. Complete lvid2record map by exchanging the vertex_info. 
     *
     * 4. Exchange the negotiation records, including singletons. (Local graph 
     * 
     * handling singletons). 
     * 5. Exchange global graph statistics.
     */
    virtual void finalize() {
      rpc.full_barrier();
      if (rpc.procid() == 0) {
        logstream(LOG_EMPH) << "Finalizing Graph..." << std::endl;
      }
      // typedef typename boost::unordered_map<vertex_id_type, lvid_type>::value_type 
      //   vid2lvid_pair_type;
      typedef typename cuckoo_map_pow2<vertex_id_type, lvid_type, 3, 
                                       uint32_t>::value_type
        vid2lvid_pair_type;
      typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
        edge_buffer_type;
      typedef boost::unordered_map<vertex_id_type, vertex_negotiator_record>
        vrec_map_type;
      typedef typename vrec_map_type::value_type vrec_pair_type;
      typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;
      typedef typename buffered_exchange<vertex_info>::buffer_type 
        vinfo_buffer_type;

      
      // Flush any additional data
      edge_exchange.flush(); vertex_exchange.flush();     
      if(rpc.procid() == 0)       
        memory_info::log_usage("Post Flush");

      logstream(LOG_INFO) << "Graph Finalize: constructing local graph" << std::endl;
      { // Add all the edges to the local graph
        const size_t nedges = edge_exchange.size()+1;
        graph.local_graph.reserve_edge_space(nedges + 1);      
        edge_buffer_type edge_buffer;
        procid_t proc;
        while(edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            // Get the source_vlid;
            lvid_type source_lvid(-1);
            if(graph.vid2lvid.find(rec.source) == graph.vid2lvid.end()) {
              source_lvid = graph.vid2lvid.size();
              graph.vid2lvid[rec.source] = source_lvid;
              // graph.local_graph.resize(source_lvid + 1);
            } else source_lvid = graph.vid2lvid[rec.source];
            // Get the target_lvid;
            lvid_type target_lvid(-1);
            if(graph.vid2lvid.find(rec.target) == graph.vid2lvid.end()) {
              target_lvid = graph.vid2lvid.size();
              graph.vid2lvid[rec.target] = target_lvid;
              // graph.local_graph.resize(target_lvid + 1);
            } else target_lvid = graph.vid2lvid[rec.target];
            graph.local_graph.add_edge(source_lvid, target_lvid, rec.edata);
          } // end of loop over add edges
        } // end for loop over buffers
        edge_exchange.clear();
      }

      if(rpc.procid() == 0) 
        memory_info::log_usage("Finished populating local graph.");
      
      ASSERT_EQ(graph.vid2lvid.size(), graph.local_graph.num_vertices());
      logstream(LOG_INFO) << "Vid2lvid size: " << graph.vid2lvid.size() << "\t" << "Max lvid in edge buffer: " << graph.local_graph.maxlvid() << std::endl;
      
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
        memory_info::log_usage("Finished finalizing local graph."); 
       
      // Setup the map containing all the vertices being negotiated by
      // this machine
      vrec_map_type vrec_map;
      { // Receive any vertex data sent by other machines
        vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
        while(vertex_exchange.recv(sending_proc, vertex_buffer)) {
          foreach(const vertex_buffer_record& rec, vertex_buffer) {
            vertex_negotiator_record& negotiator_rec = vrec_map[rec.vid];
            if (vertex_combine_strategy &&  negotiator_rec.data_is_set) {
              vertex_combine_strategy(negotiator_rec.vdata, rec.vdata);
            } else {
              negotiator_rec.vdata = rec.vdata;
              negotiator_rec.data_is_set = true;
            }
          }
        }
        vertex_exchange.clear();
      } // end of loop to populate vrecmap

      if(rpc.procid() == 0)         
        memory_info::log_usage("Emptied vertex data exchange");



      { // Compute the mirroring information for all vertices
        // negotiated by this machine
        buffered_exchange<vertex_info> vinfo_exchange(rpc.dc());
        vinfo_buffer_type recv_buffer; procid_t sending_proc = -1;
        size_t iter = 0; size_t last_iter = graph.vid2lvid.size() - 1;

        // special case: flush there is no vertex info to send.
        if (graph.vid2lvid.size() == 0) {
          vinfo_exchange.flush();
          // recv any buffers 
          while(vinfo_exchange.recv(sending_proc, recv_buffer)) {
            foreach(vertex_info vinfo, recv_buffer) {
              vertex_negotiator_record& rec = vrec_map[vinfo.vid];
              rec.num_in_edges += vinfo.num_in_edges;
              rec.num_out_edges += vinfo.num_out_edges;
              rec.mirrors.set_bit(sending_proc);
            } // end of for loop over all vertex info
          } // end of recv while loop
        } else {
          // usual case : send vertex info and receive when necessary.
          foreach(const vid2lvid_pair_type& pair, graph.vid2lvid) {
            // Send a vertex
            const vertex_id_type vid = pair.first;
            const vertex_id_type lvid = pair.second;
            const procid_t negotiator = graph_hash::hash_vertex(vid) % rpc.numprocs();
            const vertex_info vinfo(vid, graph.local_graph.num_in_edges(lvid),
                                    graph.local_graph.num_out_edges(lvid));
            vinfo_exchange.send(negotiator, vinfo);
            if (iter == last_iter)  vinfo_exchange.flush();
            // recv any buffers if necessary
            while(vinfo_exchange.recv(sending_proc, recv_buffer)) {
              foreach(vertex_info vinfo, recv_buffer) {
                vertex_negotiator_record& rec = vrec_map[vinfo.vid];
                rec.num_in_edges += vinfo.num_in_edges;
                rec.num_out_edges += vinfo.num_out_edges;
                rec.mirrors.set_bit(sending_proc);
              } // end of for loop over all vertex info
            } // end of recv while loop
            ++iter;
          } // end of loop over edge info        
        } // end of compute mirror information
      }

      if(rpc.procid() == 0) 
        memory_info::log_usage("Exchanged basic vertex info");

      { // Determine masters for all negotiated vertices
        logstream(LOG_INFO) 
          << "Graph Finalize: Constructing and sending vertex assignments" 
          << std::endl;

        std::vector<size_t> counts(rpc.numprocs()); 
        size_t num_singletons = 0;
        // For simplicity simply assign this machine as the master
        foreach(vrec_pair_type& pair, vrec_map) {
          vertex_negotiator_record& rec = pair.second;
          procid_t master = rpc.procid();
          rec.owner = master; 
          rec.mirrors.clear_bit(master);
          // singleton that does not have in/out edges on this machine 
          if (graph.vid2lvid.find(pair.first) == graph.vid2lvid.end())
            ++num_singletons;
        }

        if(rpc.procid() == 0) 
          memory_info::log_usage("Finished computing masters");

        { // Initialize vertex records
          graph.lvid2record.reserve(graph.vid2lvid.size() + num_singletons);
          graph.lvid2record.resize(graph.vid2lvid.size() + num_singletons);
          foreach(const vid2lvid_pair_type& pair, graph.vid2lvid) 
            graph.lvid2record[pair.second].gvid = pair.first;      
          // Check conditions on graph
          // ASSERT_EQ(graph.local_graph.num_vertices(), graph.lvid2record.size());
          graph.local_graph.reserve(graph.local_graph.num_vertices() + num_singletons);
        }
        if(rpc.procid() == 0)       
          memory_info::log_usage("Finished lvid2record");


        // Exchange the negotiation records
        typedef std::pair<vertex_id_type, vertex_negotiator_record> 
          exchange_pair_type;
        typedef buffered_exchange<exchange_pair_type> 
          negotiator_exchange_type;
        negotiator_exchange_type negotiator_exchange(rpc.dc(), 1, 1000);
        typename negotiator_exchange_type::buffer_type recv_buffer;
        procid_t sending_proc(-1);
        size_t iter = 0; size_t last_iter = vrec_map.size() - 1;

        // special case : flush when there is no negotiate record to send.
        if (vrec_map.size() == 0) {
          negotiator_exchange.flush();
          // receive records.
          while(negotiator_exchange.recv(sending_proc, recv_buffer)) {
              foreach(const exchange_pair_type& pair, recv_buffer) {
                const vertex_id_type& vid = pair.first;
                const vertex_negotiator_record& negotiator_rec = pair.second;
                // Determine the local vid 
                lvid_type lvid(-1);
                if(graph.vid2lvid.find(vid) == graph.vid2lvid.end()) {
                  lvid = graph.vid2lvid.size();
                  graph.vid2lvid[vid] = lvid;
                  graph.local_graph.add_vertex(lvid, negotiator_rec.vdata);
                  ASSERT_LT(lvid, graph.lvid2record.size());
                  graph.lvid2record[lvid].gvid = vid;
                } else {
                  lvid = graph.vid2lvid[vid];
                  ASSERT_LT(lvid, graph.local_graph.num_vertices());
                  graph.local_graph.vertex_data(lvid) = negotiator_rec.vdata;
                }
                ASSERT_LT(lvid, graph.lvid2record.size());
                vertex_record& local_record = graph.lvid2record[lvid];
                local_record.owner = negotiator_rec.owner;
                ASSERT_EQ(local_record.num_in_edges, 0); 
                local_record.num_in_edges = negotiator_rec.num_in_edges;
                ASSERT_EQ(local_record.num_out_edges, 0);
                local_record.num_out_edges = negotiator_rec.num_out_edges;
                local_record._mirrors = negotiator_rec.mirrors;
              }
            } // end of while loop over negotiator_exchange.recv
        } else {
          // usual case : send and receive negotiate record.
          foreach(vrec_pair_type& pair, vrec_map) {
            const vertex_id_type& vid = pair.first;
            const vertex_negotiator_record& negotiator_rec = pair.second;
            const exchange_pair_type exchange_pair(vid, negotiator_rec);
            // send to the owner
            negotiator_exchange.send(negotiator_rec.owner, exchange_pair);
            // send to all the mirrors
            foreach(size_t mirror, negotiator_rec.mirrors) {
              negotiator_exchange.send(mirror, exchange_pair);
            }
            if (iter == last_iter) negotiator_exchange.flush();
            // Recevie any records
            while(negotiator_exchange.recv(sending_proc, recv_buffer)) {
              foreach(const exchange_pair_type& pair, recv_buffer) {
                const vertex_id_type& vid = pair.first;
                const vertex_negotiator_record& negotiator_rec = pair.second;
                // Determine the local vid 
                lvid_type lvid(-1);
                if(graph.vid2lvid.find(vid) == graph.vid2lvid.end()) {
                  lvid = graph.vid2lvid.size();
                  graph.vid2lvid[vid] = lvid;
                  graph.local_graph.add_vertex(lvid, negotiator_rec.vdata);
                  ASSERT_LT(lvid, graph.lvid2record.size());
                  graph.lvid2record[lvid].gvid = vid;
                } else {
                  lvid = graph.vid2lvid[vid];
                  ASSERT_LT(lvid, graph.local_graph.num_vertices());
                  graph.local_graph.vertex_data(lvid) = negotiator_rec.vdata;
                }
                ASSERT_LT(lvid, graph.lvid2record.size());
                vertex_record& local_record = graph.lvid2record[lvid];
                local_record.owner = negotiator_rec.owner;
#ifndef USE_DYNAMIC_LOCAL_GRAPH
                ASSERT_EQ(local_record.num_in_edges, 0); 
                ASSERT_EQ(local_record.num_out_edges, 0);
#endif
                local_record.num_in_edges = negotiator_rec.num_in_edges;
                local_record.num_out_edges = negotiator_rec.num_out_edges;
                local_record._mirrors = negotiator_rec.mirrors;
              }
            } // end of while loop over negotiator_exchange.recv
            ++iter;
          } // end of for loop over local vertex records        
        }
      } // end of master exchange

      if(rpc.procid() == 0) 
        memory_info::log_usage("Finished sending updating mirrors");


      ASSERT_EQ(graph.vid2lvid.size(), graph.local_graph.num_vertices());
      ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
 
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


      if (rpc.procid() == 0) {
        logstream(LOG_EMPH) << "Graph info: "  
                            << "\n\t nverts: " << graph.num_vertices()
                            << "\n\t nedges: " << graph.num_edges()
                            << "\n\t nreplicas: " << graph.nreplicas
                            << "\n\t replication factor: " << (double)graph.nreplicas/graph.num_vertices()
                            << std::endl;
      }
    }
  }; // end of distributed_ingress_base

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
