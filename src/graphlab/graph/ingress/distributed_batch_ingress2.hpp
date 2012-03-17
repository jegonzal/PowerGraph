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

#ifndef GRAPHLAB_DISTRIBUTED_BATCH_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BATCH_INGRESS_HPP

#include <boost/unordered_set.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/idistributed_ingress.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
    class distributed_graph;

  template<typename VertexData, typename EdgeData>
  class distributed_batch_ingress : 
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData edge_data_type;
    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;
    /// Type for vertex colors 
    typedef graphlab::vertex_color_type vertex_color_type;
    typedef typename graph_type::lvid_type  lvid_type;
    typedef typename graph_type::vertex_record vertex_record;
    typedef typename graph_type::mirror_type mirror_type;

    dc_dist_object<distributed_batch_ingress> rpc;
    typedef distributed_ingress_base<VertexData, EdgeData> base_type;

    mutex local_graph_lock;
    mutex lvid2record_lock;

    typedef fixed_dense_bitset<graph_type::MAX_MACHINES> bin_counts_type;
    /** The map from vertex id to pairs of <pid, local_degree_of_v> */
    typedef typename boost::unordered_map<vertex_id_type, bin_counts_type> 
    dht_degree_table_type;

    /** distributed hash table stored on local machine */ 
    std::vector<bin_counts_type > dht_degree_table;

    // must be called with a readlock acquired on dht_degree_table_lock
    size_t vid_to_dht_entry_with_readlock(vertex_id_type vid) {
      size_t idx = (vid - rpc.procid()) / rpc.numprocs();
      if (dht_degree_table.size() <= idx) {
        dht_degree_table_lock.unlock();
        dht_degree_table_lock.writelock();
        size_t newsize = std::max(dht_degree_table.size() * 2, idx + 1);
        dht_degree_table.resize(newsize);
        dht_degree_table_lock.unlock();
        dht_degree_table_lock.readlock();
      }
      return idx;
    }
    rwlock dht_degree_table_lock;

    /** Local minibatch buffer */
    size_t num_edges; // number of edges in the current buffer
    size_t bufsize; // capacity of the local buffer 
    std::vector<std::pair<vertex_id_type, vertex_id_type> > edgesend;
    std::vector<EdgeData> edatasend;
    mutex edgesend_lock;
    std::vector<boost::unordered_set<vertex_id_type> > query_set;
    /** The map from proc_id to num_edges on that proc */
    std::vector<size_t> proc_num_edges;

    DECLARE_TRACER(batch_ingress_add_edge);
    DECLARE_TRACER(batch_ingress_add_edges);
    DECLARE_TRACER(batch_ingress_compute_assignments);
    DECLARE_TRACER(batch_ingress_request_degree_table);
    DECLARE_TRACER(batch_ingress_get_degree_table);
    DECLARE_TRACER(batch_ingress_update_degree_table);

    bool usehash;
    bool userecent; 

  public:
    distributed_batch_ingress(distributed_control& dc, graph_type& graph, 
        size_t bufsize = 50000, bool usehash = false, bool userecent = false) :
      base_type(dc, graph), rpc(dc, this), 
      num_edges(0), bufsize(bufsize), query_set(dc.numprocs()),
      proc_num_edges(dc.numprocs()), usehash(usehash), userecent(userecent) { 
       rpc.barrier(); 

      INITIALIZE_TRACER(batch_ingress_add_edge, "Time spent in add edge");
      INITIALIZE_TRACER(batch_ingress_add_edges, "Time spent in add block edges" );
      INITIALIZE_TRACER(batch_ingress_compute_assignments, "Time spent in compute assignment");
      INITIALIZE_TRACER(batch_ingress_request_degree_table, "Time spent in requesting assignment");
      INITIALIZE_TRACER(batch_ingress_get_degree_table, "Time spent in retrieve degree table");
      INITIALIZE_TRACER(batch_ingress_update_degree_table, "Time spent in update degree table");
     }

// override base::add_edge
    void add_edge(vertex_id_type source, vertex_id_type target, const EdgeData& edata) {
      BEGIN_TRACEPOINT(batch_ingress_add_edge);
      edgesend_lock.lock();
      ASSERT_LT(edgesend.size(), bufsize);
      edgesend.push_back(std::make_pair(source, target)); 
      edatasend.push_back(edata);        
      query_set[base_type::vertex_to_proc(source)].insert(source);
      query_set[base_type::vertex_to_proc(target)].insert(target);
      ++num_edges;
      edgesend_lock.unlock();
      END_TRACEPOINT(batch_ingress_add_edge);
      if (is_full()) flush();
    } // end of add_edge

    // overide base finalize; 
    void finalize() { 
      flush(); 
      base_type::finalize();
    } // end of finalize


  private:

    // HELPER ROUTINES =======================================================>    
    // Add a block of edges to storage.
    void add_edges(const std::vector<vertex_id_type>& source_arr, 
        const std::vector<vertex_id_type>& target_arr, 
        const std::vector<EdgeData>& edata_arr) {
      BEGIN_TRACEPOINT(batch_ingress_add_edges);
      ASSERT_TRUE((source_arr.size() == target_arr.size())
          && (source_arr.size() == edata_arr.size())); 
      if (source_arr.size() == 0) return;

      std::vector<lvid_type> local_source_arr; 
      local_source_arr.reserve(source_arr.size());
      std::vector<lvid_type> local_target_arr;
      local_target_arr.reserve(source_arr.size());
      /** The map from vertex_id to its degree on this proc.*/
      std::vector<std::vector<vertex_id_type> > local_degree_count(rpc.numprocs());

      lvid_type max_lvid = 0;

      lvid2record_lock.lock();
      for (size_t i = 0; i < source_arr.size(); ++i) {
        vertex_id_type source = source_arr[i];
        vertex_id_type target = target_arr[i]; 
        lvid_type lvid_source;
        lvid_type lvid_target;
        typedef typename boost::unordered_map<vertex_id_type, lvid_type>::iterator 
          vid2lvid_iter;
        vid2lvid_iter iter;

          iter = base_type::graph.vid2lvid.find(source);
          if (iter == base_type::graph.vid2lvid.end()) {
            lvid_source = base_type::graph.vid2lvid.size();
            base_type::graph.vid2lvid.insert(std::make_pair(source, lvid_source));
            base_type::graph.lvid2record.push_back(vertex_record(source));
          } else {
            lvid_source = iter->second;
          }

          iter = base_type::graph.vid2lvid.find(target);
          if (iter == base_type::graph.vid2lvid.end()) {
            lvid_target = base_type::graph.vid2lvid.size();
            base_type::graph.vid2lvid.insert(std::make_pair(target , lvid_target));
            base_type::graph.lvid2record.push_back(vertex_record(target));
          } else {
            lvid_target = iter->second;
          }

        local_source_arr.push_back(lvid_source);
        local_target_arr.push_back(lvid_target);
        max_lvid = std::max(std::max(lvid_source, lvid_target), 
            max_lvid);

        local_degree_count[base_type::vertex_to_proc(source)].push_back(source);
        local_degree_count[base_type::vertex_to_proc(target)].push_back(target);
      }
      lvid2record_lock.unlock();
      // Send out local_degree count;
      for (size_t i = 0; i < rpc.numprocs(); ++i) {
        if (i != rpc.procid()) {
          rpc.remote_call(i, 
                          &distributed_batch_ingress::block_add_degree_counts, 
                          rpc.procid(),
                          local_degree_count[i]);
        } else {
          block_add_degree_counts(rpc.procid(), local_degree_count[i]);
        }
        local_degree_count[i].clear();
      }

      // Add edges to local graph
      local_graph_lock.lock();
      if (max_lvid > 0 && max_lvid >= base_type::graph.local_graph.num_vertices()) {
        base_type::graph.local_graph.resize(max_lvid + 1);
      }
      base_type::graph.local_graph.add_edges(local_source_arr, 
                                             local_target_arr, edata_arr);
      local_graph_lock.unlock();
 
      END_TRACEPOINT(batch_ingress_add_edges);
    } // end of add edges

    void block_add_degree_counts (procid_t pid, std::vector<vertex_id_type>& whohas) {
      BEGIN_TRACEPOINT(batch_ingress_update_degree_table);
      dht_degree_table_lock.readlock();
      foreach (vertex_id_type& vid, whohas) {
        size_t idx = vid_to_dht_entry_with_readlock(vid);
        dht_degree_table[idx].set_bit_unsync(pid);
      }
      dht_degree_table_lock.unlock();
      END_TRACEPOINT(batch_ingress_update_degree_table);
    }


    dht_degree_table_type 
    block_get_degree_table(const boost::unordered_set<vertex_id_type>& vid_query) {
      BEGIN_TRACEPOINT(batch_ingress_get_degree_table);
      dht_degree_table_type answer;
      dht_degree_table_lock.readlock();
      foreach (vertex_id_type qvid, vid_query) {
        answer[qvid] = dht_degree_table[vid_to_dht_entry_with_readlock(qvid)]; 
      }
      dht_degree_table_lock.unlock();
      END_TRACEPOINT(batch_ingress_get_degree_table);
      return answer;
    }  // end of block get degree table


   void assign_edges(std::vector<std::vector<vertex_id_type> >& proc_src,
                     std::vector<std::vector<vertex_id_type> >& proc_dst,
                     std::vector<std::vector<EdgeData> >& proc_edata) {
     ASSERT_EQ(num_edges, edgesend.size());
     if (num_edges == 0) return;

     edgesend_lock.lock();
     // Get the degree table.
     BEGIN_TRACEPOINT(batch_ingress_request_degree_table);
     std::vector<dht_degree_table_type> degree_table(rpc.numprocs());
     
     for (size_t i = 0; i < rpc.numprocs(); ++i) {
       if (i == rpc.procid()) {
         degree_table[i] = block_get_degree_table(query_set[i]);
       } else {
         degree_table[i] = 
           rpc.remote_request(i, 
               &distributed_batch_ingress::block_get_degree_table,
               query_set[i]);
       }
       query_set[i].clear();
     }
     END_TRACEPOINT(batch_ingress_request_degree_table);

     for (size_t i = 0; i < num_edges; ++i) {
       std::pair<vertex_id_type, vertex_id_type>& e = 
         edgesend[i];

       BEGIN_TRACEPOINT(batch_ingress_compute_assignments);
       size_t src_proc = base_type::vertex_to_proc(e.first);
       size_t dst_proc = base_type::vertex_to_proc(e.second);
       bin_counts_type& src_degree = degree_table[src_proc][e.first];
       bin_counts_type& dst_degree = degree_table[dst_proc][e.second];
       procid_t proc = base_type::edge_decision.edge_to_proc_greedy(e.first, e.second, 
           src_degree, dst_degree, proc_num_edges, usehash, userecent);
       END_TRACEPOINT(batch_ingress_compute_assignments);

       ASSERT_LT(proc, proc_src.size());
       proc_src[proc].push_back(e.first);
       proc_dst[proc].push_back(e.second);
       proc_edata[proc].push_back(edatasend[i]);
     }
     edgesend.clear();
     edatasend.clear();
     edgesend_lock.unlock();
   } // end add edge

    // Flush all edges in the buffer.
    void flush() {
      std::vector< std::vector<vertex_id_type> > proc_src(rpc.numprocs());
      std::vector< std::vector<vertex_id_type> > proc_dst(rpc.numprocs());
      std::vector< std::vector<EdgeData> > proc_edata(rpc.numprocs());
      assign_edges(proc_src, proc_dst, proc_edata);
      for (size_t i = 0; i < proc_src.size(); ++i) {
        if (proc_src[i].size() == 0) 
          continue;
        if (i == rpc.procid()) {
          add_edges(proc_src[i], proc_dst[i], proc_edata[i]);
          num_edges -= proc_src[i].size();
        } else {
          rpc.remote_call(i, &distributed_batch_ingress::add_edges,
              proc_src[i], proc_dst[i], proc_edata[i]);
          num_edges -= proc_src[i].size();
        } // end if
      } // end for
    } // end flush

    size_t size() { return num_edges; }
    bool is_full() { return size() >= bufsize; }

  }; // end of distributed_batch_ingress

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
