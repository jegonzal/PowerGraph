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
 *  software distributed under the License is distributed on an "AS *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
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
#include <graphlab/logger/logger.hpp>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <graphlab/macros_def.hpp>

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


    typedef typename hopscotch_map<vertex_id_type, lvid_type>::value_type
        vid2lvid_pair_type;    
    
    typedef typename buffered_exchange<vertex_id_type>::buffer_type
        vertex_id_buffer_type;

    typedef typename graph_type::zone_type zone_type;

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

    typedef typename std::pair<procid_t,int> mht_entry_type;

    /** Type of the master location hash table:
     * a map from vertex id to the location of its master. */
    typedef typename boost::unordered_map<vertex_id_type, procid_t> master_hash_table_type;
    typedef typename std::pair<vertex_id_type, procid_t> master_pair_type;
    typedef typename buffered_exchange<master_pair_type>::buffer_type master_buffer_type;

    typedef typename std::pair<procid_t,size_t> proc_edges_pair_type;
    typedef typename buffered_exchange<proc_edges_pair_type >::buffer_type
        proc_edges_buffer_type;

    /** Number of edges and vertices in graph. */
    buffered_exchange<edge_buffer_record> hybrid_edge_exchange;
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

    buffered_exchange<edge_buffer_record> high_edge_exchange;
    buffered_exchange<edge_buffer_record>  low_edge_exchange;

    //these are used for edges balance
    
    buffered_exchange< proc_edges_pair_type >  proc_edges_exchange;
    std::vector<size_t> proc_edges_incremental;

    std::vector<size_t> proc_num_vertices;

    /// threshold to divide high-degree and low-degree vertices
    size_t threshold;
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
        hybrid_edge_exchange(dc), hybrid_vertex_exchange(dc),temporary_vertex_exchange(dc),
        master_exchange(dc),
        hybrid_rpc(dc, this),graph(graph),high_edge_exchange(dc),low_edge_exchange(dc),proc_edges_exchange(dc),
        proc_edges_incremental(dc.numprocs()),proc_num_vertices(dc.numprocs()), threshold(threshold), 
        nedges(nedges), nverts(nverts), interval(interval) { 
      gamma = 1.5;
      alpha = sqrt(dc.numprocs()) * nedges / pow(nverts, gamma);
      if(nedges<=0)
        nedges=1;
    } // end of constructor

    ~distributed_hybrid_ginger_ingress() { }

    /** Add an edge to the ingress object using random hashing assignment.
     *  This function acts as the first phase for SNAP graph to deliver edges
     *  via the hashing value of its target vertex.
     */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      const procid_t owning_proc = 
        graph_hash::hash_vertex(target) % hybrid_rpc.numprocs();
      const edge_buffer_record record(source, target, edata);
      hybrid_edge_exchange.send(owning_proc, record);
    } // end of add edge

    /** Add edges to the ingress object using different assignment policies.
     *  This function handles the RADJ graph in different ways.
     *  For high degree edges, hashing from its source vertex;
     *  for low degree edges, using a heuristic way named ginger.
     */
    void add_edges(std::vector<vertex_id_type>& in, vertex_id_type vid,
                  const std::vector<EdgeData>& edatas) {
      procid_t owning_proc;
      if (in.size()>=threshold) {
        owning_proc = graph_hash::hash_vertex(vid) % base_type::rpc.numprocs();
      } else {
        owning_proc = 
            base_type::edge_decision.edge_to_proc_ginger(vid, mht, 
            mht_incremental, proc_num_vertices, alpha, gamma, in);
      }

      typedef typename base_type::edge_buffer_record edge_buffer_record;
      if(in.size()>=threshold){
        for (size_t i = 0; i < in.size(); ++i) {
          edge_buffer_record record(in[i], vid, edatas[i]);
          high_edge_exchange.send(owning_proc, record);
        }
      }
      else{
        for (size_t i = 0; i < in.size(); ++i) {
          edge_buffer_record record(in[i], vid, edatas[i]);
          low_edge_exchange.send(owning_proc, record);
        }

        proc_edges_incremental[owning_proc]+=in.size(); //edge balance
      }
      size_t numprocs = base_type::rpc.numprocs();
      mht_incremental[vid] = owning_proc;
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
        for(procid_t p=0;p<hybrid_rpc.numprocs();p++){
          for(procid_t i=0;i<hybrid_rpc.numprocs();i++){
            proc_edges_exchange.send(p,std::make_pair(i,proc_edges_incremental[i]));
          }
        }
        proc_edges_exchange.partial_flush(0);
        
        proc_edges_buffer_type proc_edge_buffer;
        while(proc_edges_exchange.recv(proc, proc_edge_buffer)) {
          foreach(const proc_edges_pair_type& pair, proc_edge_buffer) {
            proc_edges_incremental[pair.first]+=pair.second;
          }
        }
        proc_edges_exchange.clear();

        for(int i=0;i<hybrid_rpc.numprocs();i++){
          proc_num_vertices[i]+=proc_edges_incremental[i]*1.0*nverts/nedges;
          proc_edges_incremental[i]=0;
        }
      }
    } // end of add edges

    /* add vdata */
    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      const procid_t owning_proc = 
        graph_hash::hash_vertex(vid) % hybrid_rpc.numprocs();
      const vertex_buffer_record record(vid, vdata);
      temporary_vertex_exchange.send(owning_proc, record);
    } // end of add vertex


    void snap_ingress(){
      typedef typename boost::unordered_map<vertex_id_type,
      batch_edge_buffer_record> batch_record_map_type;
      batch_record_map_type batch_map;

      hybrid_edge_exchange.flush();

      edge_buffer_type edge_buffer;
      procid_t proc = -1;
      while(hybrid_edge_exchange.recv(proc, edge_buffer)) {
        foreach(const edge_buffer_record& rec, edge_buffer) {
          batch_map[rec.target].sources.push_back(rec.source);
          batch_map[rec.target].edatas.push_back(rec.edata);
        }
      }
      hybrid_edge_exchange.clear();

      hybrid_rpc.full_barrier();

      for (typename batch_record_map_type::iterator it = batch_map.begin();it != batch_map.end(); ++it) {
        add_edges(it->second.sources, it->first,it->second.edatas);
      }
    }

    void finalize() {
      graphlab::timer ti;
      
      snap_ingress();

      size_t nprocs = hybrid_rpc.numprocs();
      procid_t l_procid = hybrid_rpc.procid();
      std::vector<edge_buffer_record> hybrid_edges;
      size_t nedges;


      {  // send mht to all other nodes
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
        procid_t proc;
        while(master_exchange.recv(proc, master_buffer)) {
          foreach(const master_pair_type& pair, master_buffer) {
            mht[pair.first] = pair.second;
          }
        }
        master_exchange.clear();
  
        

        std::vector<edge_buffer_record> hybrid_edges;
        hopscotch_map<vertex_id_type, size_t> in_degree_set;
        
        
        if (l_procid == 0) {
          memory_info::log_usage("start finalizing");
          logstream(LOG_EMPH) << "hybrid ginger finalizing Graph ..."
                              << " #vertices=" << graph.local_graph.num_vertices()
                              << " #edges=" << graph.local_graph.num_edges()
                              << " threshold=" << threshold
                              << std::endl;
        }
      }
  




      /**************************************************************************/
      /*                                                                        */
      /*                         batch ingress                         */
      /*                                                                        */
      /**************************************************************************/
      {
        high_edge_exchange.flush();
        low_edge_exchange.flush();

        nedges=low_edge_exchange.size();
        hybrid_edges.reserve(nedges+high_edge_exchange.size());

        edge_buffer_type edge_buffer;
        procid_t proc = -1;
        while(low_edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            if (mht.find(rec.source) == mht.end())
                mht[rec.source] = graph_hash::hash_vertex(rec.source) % hybrid_rpc.numprocs();   
            hybrid_edges.push_back(rec);
          }
        }
        low_edge_exchange.clear();
  
        hybrid_rpc.full_barrier();
  
        while(high_edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            if (mht.find(rec.source) == mht.end())
                mht[rec.source] = graph_hash::hash_vertex(rec.source) % hybrid_rpc.numprocs(); 
            const procid_t source_owner_proc = mht[rec.source];
  
            if(source_owner_proc == l_procid){
              hybrid_edges.push_back(rec);
              nedges++;
            } else {
              low_edge_exchange.send(source_owner_proc, rec);
            }
          }
        }
        high_edge_exchange.clear();
  
        low_edge_exchange.flush();
  
        if(l_procid == 0)
          logstream(LOG_INFO) << "receive high-degree mirrors: "  
                              << low_edge_exchange.size() << std::endl;
        {
          edge_buffer_type edge_buffer;
          procid_t proc = -1;
          while(low_edge_exchange.recv(proc, edge_buffer)) {
            foreach(const edge_buffer_record& rec, edge_buffer) {
              mht[rec.source] = l_procid;
              hybrid_edges.push_back(rec);
              nedges++;
            }
          }
        }
        low_edge_exchange.clear();
      }

      // connect to base finalize()
      // pass nedges and hybrid_edges to the base finalize
      modified_base_finalize(nedges, hybrid_edges);

      set_vertex_type();
      
      if(l_procid == 0)
        logstream(LOG_EMPH) << "orignal hybrid finalizing graph done. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
  } // end of finalize

    void set_vertex_type() {
      graphlab::timer ti;
      procid_t l_procid = hybrid_rpc.procid();

      for (size_t lvid = 0; lvid < graph.num_local_vertices(); lvid++) {
        vertex_record& vrec = graph.lvid2record[lvid];
        if (vrec.owner == l_procid) {
          if (vrec.num_in_edges >= threshold) 
            vrec.type = graph_type::HIGH_MASTER;
          else vrec.type = graph_type::LOW_MASTER;
        } else {
          if (vrec.num_in_edges >= threshold) 
            vrec.type = graph_type::HIGH_MIRROR;
          else vrec.type = graph_type::LOW_MIRROR;
        }
      }

      if(l_procid == 0) {
        memory_info::log_usage("set vertex type done."); 
        logstream(LOG_EMPH) << "set vertex type: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
    }

    /* do the same job as original base finalize except for
     * extracting edges from hybrid_edges instead of original edge_buffer;
     * and using mht to tracing the master location of each vertex.
     */
    void modified_base_finalize(size_t nedges,
        std::vector<edge_buffer_record>& hybrid_edges) {

      graphlab::timer ti;
      procid_t l_procid = hybrid_rpc.procid();

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
      /*                       flush any additional data                        */
      /*                                                                        */
      /**************************************************************************/
      temporary_vertex_exchange.flush();
      

      /**************************************************************************/
      /*                                                                        */
      /*                         Construct local graph                          */
      /*                                                                        */
      /**************************************************************************/
      { 
        // Add all the edges to the local graph
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
        if(l_procid == 0)  {
          memory_info::log_usage("populating local graph done.");
          logstream(LOG_EMPH) << "populating local graph: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }

        // Finalize local graph
        graph.local_graph.finalize();
        logstream(LOG_INFO) << "local graph info: " << std::endl
                            << "\t nverts: " << graph.local_graph.num_vertices()
                            << std::endl
                            << "\t nedges: " << graph.local_graph.num_edges()
                            << std::endl;
        
        if(l_procid == 0) {
          memory_info::log_usage("finalizing local graph done."); 
          logstream(LOG_EMPH) << "finalizing local graph: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }
      }


      /**************************************************************************/
      /*                                                                        */
      /*             receive and add vertex data to masters                     */
      /*                                                                        */
      /**************************************************************************/
      // Setup the map containing all the vertices being negotiated by this machine
      {
        if (temporary_vertex_exchange.size() > 0) {
          vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
          while(temporary_vertex_exchange.recv(sending_proc, vertex_buffer)) {
            foreach(const vertex_buffer_record& rec, vertex_buffer) {
              if (mht.find(rec.vid) == mht.end())
                mht[rec.vid] = graph_hash::hash_vertex(rec.vid) % hybrid_rpc.numprocs(); 
              hybrid_vertex_exchange.send(mht[rec.vid], rec);
            }
          }
          temporary_vertex_exchange.clear();
        }
        hybrid_vertex_exchange.flush();

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

          logstream(LOG_INFO) << " #vert-msgs=" << hybrid_vertex_exchange.size()
                              << std::endl;
          if(l_procid == 0) {        
            memory_info::log_usage("adding vertex data done.");
            logstream(LOG_EMPH) << "adding vertex data: " 
                                << ti.current_time()
                                << " secs" 
                                << std::endl;
          }
        }
      } // end of loop to populate vrecmap


      /**************************************************************************/
      /*                                                                        */
      /*        assign vertex data and allocate vertex (meta)data  space        */
      /*                                                                        */
      /**************************************************************************/

#ifdef INGRESS_DEBUG
      std::vector<size_t> degree;
      degree.resize(graph_type::NUM_ZONE_TYPES, 0);
#endif
      {
        // determine masters for all negotiated vertices
        const size_t local_nverts = graph.vid2lvid.size() + vid2lvid_buffer.size();
        graph.lvid2record.reserve(local_nverts);
        graph.lvid2record.resize(local_nverts);
        graph.local_graph.resize(local_nverts);
        foreach(const vid2lvid_pair_type& pair, vid2lvid_buffer) {
          vertex_record& vrec = graph.lvid2record[pair.second];
          vrec.gvid = pair.first;
          if (mht.find(pair.first) == mht.end())
              mht[pair.first] = graph_hash::hash_vertex(pair.first) % hybrid_rpc.numprocs(); 
          const procid_t source_owner_proc = mht[pair.first];
          vrec.owner = source_owner_proc;
        }
        ASSERT_EQ(local_nverts, graph.local_graph.num_vertices());
        ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
        if(l_procid == 0) {
          memory_info::log_usage("allocating lvid2record done.");
          logstream(LOG_EMPH) << "allocating lvid2record: " 
                              << ti.current_time()
                              << " secs" 
                              << std::endl;
        }

        mht.clear();
      }

      /**************************************************************************/
      /*                                                                        */
      /*                          master handshake                              */
      /*                                                                        */
      /**************************************************************************/      
      {
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

      if(l_procid == 0) {
        memory_info::log_usage("master handshake done.");
        logstream(LOG_EMPH) << "master handshake: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }


      /**************************************************************************/
      /*                                                                        */
      /*                        merge in vid2lvid_buffer                        */
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

      base_type::exchange_global_info();
      if(l_procid == 0) {
        memory_info::log_usage("exchange global info done.");
        logstream(LOG_EMPH) << "exchange global info: " 
                            << ti.current_time()
                            << " secs" 
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
  }; // end of distributed_hybrid_ginger_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
