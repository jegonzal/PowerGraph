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
 * 2014.04  implement bipartite-aware partitioning with affinity
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_BIPARTITE_AFFINITY_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BIPARTITE_AFFINITY_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>

#include <graphlab/macros_def.hpp>
#include <map>
#include <set>
#include <vector>
#include <graphlab/parallel/pthread_tools.hpp>

#define TUNING
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges with data affinity for bipartite graph.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_bipartite_affinity_ingress : 
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData edge_data_type;


    typedef distributed_ingress_base<VertexData, EdgeData> base_type;

    typedef typename graph_type::vertex_record vertex_record;

    typedef typename base_type::edge_buffer_record edge_buffer_record;
    typedef typename base_type::vertex_buffer_record vertex_buffer_record;
    
    /// The rpc interface for this object
    dc_dist_object<distributed_bipartite_affinity_ingress> bipartite_rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;

    simple_spinlock bipartite_vertex_lock;
    std::vector<vertex_buffer_record> bipartite_vertexs;
    simple_spinlock bipartite_edge_lock;
    std::vector<edge_buffer_record> bipartite_edges;

    bool favorite_source;

    typedef typename boost::unordered_map<vertex_id_type, procid_t> 
        master_hash_table_type;
    typedef typename std::pair<vertex_id_type, procid_t> 
        master_pair_type;
    typedef typename buffered_exchange<master_pair_type>::buffer_type 
        master_buffer_type;

    master_hash_table_type mht;
    buffered_exchange<master_pair_type> mht_exchange;

  public:
    distributed_bipartite_affinity_ingress(distributed_control& dc, graph_type& graph, const std::string& favorite):
    base_type(dc, graph), bipartite_rpc(dc, this), graph(graph), mht_exchange(dc) {
      favorite_source = favorite == "source" ? true : false;
    } // end of constructor

    ~distributed_bipartite_affinity_ingress() { }

    /** Add an edge to the ingress object using random assignment. */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      const edge_buffer_record record(source, target, edata);
      bipartite_edge_lock.lock();
      bipartite_edges.push_back(record);
      bipartite_edge_lock.unlock();
    } // end of add edge
    
    void add_vertex(vertex_id_type vid, const VertexData& vdata) { 
      const vertex_buffer_record record(vid, vdata);
      bipartite_vertex_lock.lock();
      bipartite_vertexs.push_back(record);
      mht[vid] = bipartite_rpc.procid();
      bipartite_vertex_lock.unlock();
    } // end of add vertex

    void finalize() {
      graphlab::timer ti;
      
      size_t nprocs = bipartite_rpc.numprocs();
      procid_t l_procid = bipartite_rpc.procid();

      bipartite_rpc.full_barrier();

      if (l_procid == 0) {
        memory_info::log_usage("start finalizing");
        logstream(LOG_EMPH) << "bipartite w/ affinity finalizing ..."
                            << " #verts=" << graph.local_graph.num_vertices()
                            << " #edges=" << graph.local_graph.num_edges()
                            << " favorite=" << (favorite_source ? "source" : "target")
                            << std::endl;
      }


      /**
       * Fast pass for redundant finalization with no graph changes. 
       */
      {
        size_t changed_size = bipartite_vertexs.size() + bipartite_edges.size();
        bipartite_rpc.all_reduce(changed_size);
        if (changed_size == 0) {
          logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
          return;
        }
      }

      /* directly add vertices loaded from local file to the local graph */
      const size_t local_nverts = bipartite_vertexs.size();
      graph.lvid2record.resize(local_nverts);
      graph.local_graph.resize(local_nverts);
      lvid_type lvid = 0;
      foreach(const vertex_buffer_record& rec, bipartite_vertexs){
        graph.vid2lvid[rec.vid] = lvid;
        graph.local_graph.add_vertex(lvid,rec.vdata);
        vertex_record& vrec = graph.lvid2record[lvid];
        vrec.gvid = rec.vid;
        vrec.owner = bipartite_rpc.procid();
        lvid++;
      }

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "add " << local_nverts << " vertex: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      /* exchange mapping table using mht_exchange */
      for (typename master_hash_table_type::iterator it = mht.begin(); 
          it != mht.end(); ++it) {
        for (procid_t i = 0; i < nprocs; ++i) {
          if (i != l_procid)
            mht_exchange.send(i, master_pair_type(it->first, it->second));
        }
      }

      mht_exchange.flush();
      master_buffer_type master_buffer;
      procid_t proc = -1;
      while(mht_exchange.recv(proc, master_buffer)) {
        foreach(const master_pair_type& pair, master_buffer) {
          mht[pair.first] = pair.second;
        }
      }
      mht_exchange.clear();

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "exchange mapping: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      /* distribute edges */ 
      foreach(const edge_buffer_record& rec, bipartite_edges){
        vertex_id_type favorite = favorite_source ? rec.source : rec.target;
        if(mht.find(favorite) == mht.end())
          mht[favorite] = graph_hash::hash_vertex(favorite) % nprocs;
        const procid_t owning_proc = mht[favorite];
        // save to the buffer of edge_exchange in ingress_base
        base_type::edge_exchange.send(owning_proc, rec);
      }
      bipartite_vertexs.clear();
      bipartite_edges.clear();
      mht.clear();

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "distribute edges: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif

      // call base finalize()
      base_type::finalize();
      if(l_procid == 0) {
        memory_info::log_usage("bipartite w/ affinity finalizing graph done.");
        logstream(LOG_EMPH) << "bipartite w/ affinity finalizing graph. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }
    } // end of finalize

  }; // end of distributed_bipartite_affinity_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
