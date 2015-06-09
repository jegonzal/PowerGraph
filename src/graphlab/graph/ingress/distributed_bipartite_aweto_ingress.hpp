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
 * 2014.04  implement bipartite-aware partitioning with heuristic (aweto)
 *
 */


#ifndef GRAPHLAB_DISTRIBUTED_BIPARTITE_AWETO_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BIPARTITE_AWETO_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/graph/ingress/sharding_constraint.hpp>
#include <graphlab/graph/ingress/ingress_edge_decision.hpp>


#include <graphlab/macros_def.hpp>
#include <map>
#include <set>
#include <vector>

#define TUNING
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object benefit for bipartite graph.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_bipartite_aweto_ingress : 
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
    typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
        edge_buffer_type;

    typedef typename base_type::vertex_buffer_record vertex_buffer_record;
    typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;


    /// The rpc interface for this object
    dc_dist_object<distributed_bipartite_aweto_ingress> bipartite_rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;

    std::vector<edge_buffer_record> bipartite_edges;

    bool favorite_source;

    /* ingress exchange */
    buffered_exchange<vertex_buffer_record> bipartite_vertex_exchange;
    buffered_exchange<edge_buffer_record> bipartite_edge_exchange;

  public:
    distributed_bipartite_aweto_ingress(distributed_control& dc, graph_type& graph, const std::string& favorite) :
    base_type(dc, graph), bipartite_rpc(dc, this), graph(graph),
#ifdef _OPENMP
    bipartite_vertex_exchange(dc, omp_get_max_threads()), 
    bipartite_edge_exchange(dc, omp_get_max_threads())
#else
    bipartite_vertex_exchange(dc),bipartite_edge_exchange(dc)
#endif
    {
      favorite_source = (favorite == "source") ? true : false;      
    } // end of constructor

    ~distributed_bipartite_aweto_ingress() { 
      
    }

    /** accumulate edges temporal rally point using random of "favorite" assignment. */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      vertex_id_type favorite = favorite_source ? source : target;
      const procid_t owning_proc = 
        graph_hash::hash_vertex(favorite) % bipartite_rpc.numprocs();
      const edge_buffer_record record(source, target, edata);
#ifdef _OPENMP
      bipartite_edge_exchange.send(owning_proc, record, omp_get_thread_num());
#else
      bipartite_edge_exchange.send(owning_proc, record);
#endif
    } // end of add edge

    /** accumulate edges temporal rally point using random of "favorite" assignment. */
    void add_vertex(vertex_id_type vid, const VertexData& vdata) { 
      const procid_t owning_proc =
        graph_hash::hash_vertex(vid) % bipartite_rpc.numprocs();
      const vertex_buffer_record record(vid, vdata);
#ifdef _OPENMP
      bipartite_vertex_exchange.send(owning_proc, record, omp_get_thread_num());
#else
      bipartite_vertex_exchange.send(owning_proc, record);
#endif
    } // end of add vertex

    void finalize() {
      graphlab::timer ti;

      size_t nprocs = bipartite_rpc.numprocs();
      procid_t l_procid = bipartite_rpc.procid();


      bipartite_rpc.full_barrier();

      if (l_procid == 0) {
        memory_info::log_usage("start finalizing");
        logstream(LOG_EMPH) << "bipartite aweto finalizing ..."
                            << " #verts=" << graph.local_graph.num_vertices()
                            << " #edges=" << graph.local_graph.num_edges()
                            << " favorite=" << (favorite_source ? "source" : "target")
                            << std::endl;
      }

      /**************************************************************************/
      /*                                                                        */
      /*                       Flush any additional data                        */
      /*                                                                        */
      /**************************************************************************/
      bipartite_edge_exchange.flush(); bipartite_vertex_exchange.flush();

      /**
       * Fast pass for redundant finalization with no graph changes. 
       */
      {
        size_t changed_size = bipartite_edge_exchange.size() + bipartite_vertex_exchange.size();
        bipartite_rpc.all_reduce(changed_size);
        if (changed_size == 0) {
          logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
          return;
        }
      }

      /**************************************************************************/
      /*                                                                        */
      /*     calculate the distribution of favorite vertex's neighbors    */
      /*                                                                        */
      /**************************************************************************/
      boost::unordered_map<vertex_id_type,std::vector<int> > count_map;
      edge_buffer_type edge_buffer;
      procid_t proc(-1);
      while(bipartite_edge_exchange.recv(proc, edge_buffer)) {
        foreach(const edge_buffer_record& rec, edge_buffer) {
          vertex_id_type favorite, second;
          if (favorite_source)  { favorite = rec.source; second = rec.target; }
          else { favorite = rec.target; second = rec.source; }
          
          if(count_map.find(favorite) == count_map.end())
            count_map[favorite].resize(nprocs);

          const procid_t owner_proc = graph_hash::hash_vertex(second) % nprocs;
          count_map[favorite][owner_proc] += 1;

          bipartite_edges.push_back(rec);
        }
      }
      bipartite_edge_exchange.clear();


      /**************************************************************************/
      /*                                                                        */
      /*      record heuristic location of favorite vertex      */
      /*                                                                        */
      /**************************************************************************/
      buffered_exchange<vertex_id_type> vid_buffer(bipartite_rpc.dc());
      std::set<vertex_id_type> own_vid_set;
      
      // record current nedges distributed from this machine.
      std::vector<int> proc_num_edges(nprocs);      
      boost::unordered_map<vertex_id_type, procid_t> mht;

      for(typename boost::unordered_map<vertex_id_type, std::vector<int> >::iterator it = count_map.begin();
          it != count_map.end(); ++it) {
        procid_t best_proc = l_procid;
        // heuristic score
        double best_score = (it->second)[best_proc] 
                          - sqrt(1.0*proc_num_edges[best_proc]);

        for(size_t i = 0; i < bipartite_rpc.numprocs(); i++) {
          double score = (it->second)[i] 
                               - sqrt(1.0*proc_num_edges[i]);
          if(score > best_score) {
            best_proc  = i;
            best_score = score;
          }
        }

        // update nedges
        for(size_t i = 0; i < bipartite_rpc.numprocs(); i++)
          proc_num_edges[best_proc] += (it->second)[i];

        mht[it->first] = best_proc;
        vid_buffer.send(best_proc, it->first);
      }

      // find all favorite vertices this machine own
      vid_buffer.flush();
      {
        typename buffered_exchange<vertex_id_type>::buffer_type buffer;
        procid_t recvid(-1);
        while(vid_buffer.recv(recvid, buffer)) {
          foreach(const vertex_id_type vid, buffer)
            own_vid_set.insert(vid);
        }
      }
      vid_buffer.clear();

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "hold " << own_vid_set.size() << " masters: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      /**************************************************************************/
      /*                                                                        */
      /*                          exchange edges                              */
      /*                                                                        */
      /**************************************************************************/
      for (size_t i = 0; i < bipartite_edges.size(); i++) {
        edge_buffer_record& rec = bipartite_edges[i];
        vertex_id_type favorite = favorite_source ? rec.source : rec.target;
        procid_t owner_proc = mht[favorite];        
        // save to the buffer of edge_exchange in ingress_base
        base_type::edge_exchange.send(owner_proc,rec);   
      }

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "exchange edges: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      /**************************************************************************/
      /*                                                                        */
      /*                     add vertices  to local graph                 */
      /*                                                                        */
      /**************************************************************************/
      graph.lvid2record.resize(own_vid_set.size());
      graph.local_graph.resize(own_vid_set.size());
      lvid_type lvid = 0;
      foreach(const vertex_id_type& vid, own_vid_set){
        graph.vid2lvid[vid] = lvid;
        vertex_record& vrec = graph.lvid2record[lvid];
        vrec.gvid = vid;
        vrec.owner = bipartite_rpc.procid();
        lvid++;
      }
      own_vid_set.clear();

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "add vertices to local graph: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      /**************************************************************************/
      /*                                                                        */
      /*                   re-send favorite vertex data                    */
      /*                                                                        */
      /**************************************************************************/
      {
        vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
        while(bipartite_vertex_exchange.recv(sending_proc, vertex_buffer)) {
          foreach(const vertex_buffer_record& rec, vertex_buffer) {
            if(mht.find(rec.vid) != mht.end()) {
              base_type::vertex_exchange.send(mht[rec.vid], rec);
            } else {
              base_type::vertex_exchange.send(l_procid, rec);
            }
          }
        }
        bipartite_vertex_exchange.clear();
        mht.clear();
      }

#ifdef TUNING
      if(l_procid == 0) { 
        logstream(LOG_INFO) << "exchange vertex data: " 
                            << ti.current_time()
                            << " secs" 
                            << std::endl;
      }
#endif


      // call base finalize()
      base_type::finalize();
      if(l_procid == 0) {
        memory_info::log_usage("bipartite aweto finalizing graph done.");
        logstream(LOG_EMPH) << "bipartite aweto finalizing graph. (" 
                            << ti.current_time() 
                            << " secs)" 
                            << std::endl;
      }
    } // end of finalize

  }; // end of distributed_bipartite_aweto_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
