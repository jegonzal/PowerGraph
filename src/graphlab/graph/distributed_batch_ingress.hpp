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
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
    class distributed_graph;

  template<typename VertexData, typename EdgeData>
  class distributed_batch_ingress : 
    public idistributed_ingress<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;
    /// Type for vertex colors 
    typedef graphlab::vertex_color_type vertex_color_type;

    typedef typename graph_type::lvid_type  lvid_type;
    typedef typename graph_type::vertex_record vertex_record;

    typedef typename graph_type::mirror_type mirror_type;

    dc_dist_object<distributed_batch_ingress> rpc;
    graph_type& graph;
    mutex local_graph_lock;
    mutex lvid2record_lock;

    /// Temporar buffers used to store vertex data on ingress
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


    struct shuffle_record : public graphlab::IS_POD_TYPE {
      vertex_id_type vid, num_in_edges, num_out_edges;
      shuffle_record(vertex_id_type vid = 0, vertex_id_type num_in_edges = 0,
                     vertex_id_type num_out_edges = 0) : 
        vid(vid), num_in_edges(num_in_edges), num_out_edges(num_out_edges) { }     
    }; // end of shuffle_record

    // Helper type used to synchronize the vertex data and assignments
    struct vertex_negotiator_record {
      vertex_id_type vid;
      procid_t owner;
      size_t num_in_edges, num_out_edges;
      mirror_type mirrors;
      vertex_data_type vdata;
      vertex_negotiator_record() : 
        vid(-1), owner(-1), num_in_edges(0), num_out_edges(0) { }
      void load(iarchive& arc) { 
        arc >> vid >> owner >> num_in_edges >> num_out_edges
            >> mirrors >> vdata;
      } // end of load
      void save(oarchive& arc) const { 
        arc << vid << owner << num_in_edges << num_out_edges
            << mirrors << vdata;
      } // end of save     
    }; // end of vertex_negotiator_record 


    typedef typename boost::unordered_map<vertex_id_type, size_t>  vid2degree_type;


    /// temporary map for vertexdata
    typedef boost::unordered_map<vertex_id_type, vertex_negotiator_record> vrec_map_type;
    vrec_map_type vrec_map;

    /** The map from vertex id to pairs of <pid, local_degree_of_v> */
    typedef typename boost::unordered_map<vertex_id_type, std::vector<size_t> > 
    dht_degree_table_type;
    
    std::vector<std::vector<size_t> > dht_degree_table;
    rwlock dht_degree_table_lock;
    
    // must be called with a readlock acquired on dht_degree_table_lock
    size_t vid_to_dht_entry_with_readlock(vertex_id_type vid) {
      //ASSERT_EQ((vid - rpc.procid()) % rpc.numprocs(), 0);
      size_t idx = (vid - rpc.procid()) / rpc.numprocs();
      if (dht_degree_table.size() <= idx) {
        dht_degree_table_lock.unlock();
        dht_degree_table_lock.writelock();
        size_t newsize = std::max(dht_degree_table.size() * 2, idx + 1);
        dht_degree_table.resize(newsize, std::vector<size_t>(rpc.numprocs(), 0));
        dht_degree_table_lock.unlock();
        dht_degree_table_lock.readlock();
      }
      return idx;
    }
    

    // Local minibatch buffer 
    size_t num_edges;
    size_t limit;
    std::vector<std::pair<vertex_id_type, vertex_id_type> > edgesend;
    mutex edgesend_lock;
    std::vector<EdgeData> edatasend;
    std::vector<boost::unordered_set<vertex_id_type> > query_set;


    /** The map from proc_id to num_edges on that proc */
    std::vector<size_t> proc_num_edges;

    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);
    DECLARE_TRACER(batch_ingress_add_edge);
    DECLARE_TRACER(batch_ingress_add_edges);
    DECLARE_TRACER(batch_ingress_compute_assignments);
    DECLARE_TRACER(batch_ingress_request_degree_table);
    DECLARE_TRACER(batch_ingress_get_degree_table);
    DECLARE_TRACER(batch_ingress_update_degree_table);
    DECLARE_TRACER(batch_ingress_finalize);

    enum {
      EVENT_EDGE_SEEN_NONE_UNIQUE = 0,
      EVENT_EDGE_SEEN_NONE_TIE = 1,
      EVENT_EDGE_SEEN_ONE_UNIQUE = 2,
      EVENT_EDGE_SEEN_ONE_TIE = 3,
      EVENT_EDGE_SEEN_BOTH_UNIQUE =4,
      EVENT_EDGE_SEEN_BOTH_TIE = 5
    };

  public:
    distributed_batch_ingress(distributed_control& dc, graph_type& graph, size_t bufsize = 50000) :
      rpc(dc, this), graph(graph), vertex_exchange(dc),
      num_edges(0), limit(bufsize), 
      query_set(dc.numprocs()),
      proc_num_edges(dc.numprocs()) { 
       rpc.barrier(); 

#ifdef USE_EVENT_LOG
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 500, 
                               dist_event_log::RATE_BAR);
#else
      PERMANENT_INITIALIZE_DIST_EVENT_LOG(eventlog, dc, std::cout, 500, 
                                dist_event_log::LOG_FILE);
#endif

      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, EVENT_EDGE_SEEN_NONE_UNIQUE, "Zero end (unique)");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, EVENT_EDGE_SEEN_NONE_TIE, "Zero end (tie)");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, EVENT_EDGE_SEEN_ONE_UNIQUE, "One end (unique)");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, EVENT_EDGE_SEEN_ONE_TIE, "One end (tie)");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, EVENT_EDGE_SEEN_BOTH_UNIQUE, "Both ends (unique)");
      PERMANENT_ADD_DIST_EVENT_TYPE(eventlog, EVENT_EDGE_SEEN_BOTH_TIE, "Both ends (tie)");


      INITIALIZE_TRACER(batch_ingress_add_edge, "Time spent in add edge");
      INITIALIZE_TRACER(batch_ingress_add_edges, "Time spent in add block edges" );
      INITIALIZE_TRACER(batch_ingress_compute_assignments, "Time spent in compute assignment");
      INITIALIZE_TRACER(batch_ingress_request_degree_table, "Time spent in requesting assignment");
      INITIALIZE_TRACER(batch_ingress_get_degree_table, "Time spent in retrieve degree table");
      INITIALIZE_TRACER(batch_ingress_update_degree_table, "Time spent in update degree table");
      INITIALIZE_TRACER(batch_ingress_finalize, "Time spent in finalize");
     }
void add_edge(vertex_id_type source, vertex_id_type target, const EdgeData& edata) {
      BEGIN_TRACEPOINT(batch_ingress_add_edge);
      edgesend_lock.lock();
      ASSERT_LT(edgesend.size(), limit);
      edgesend.push_back(std::make_pair(source, target)); 
      edatasend.push_back(edata);        
      query_set[vertex_to_proc(source)].insert(source);
      query_set[vertex_to_proc(target)].insert(target);
      ++num_edges;
      edgesend_lock.unlock();
      END_TRACEPOINT(batch_ingress_add_edge);
      if (is_full()) flush();
    } // end of add_edge



    // This is a local only method
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
      std::vector<vid2degree_type> local_degree_count(rpc.numprocs());
        

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

          iter = graph.vid2lvid.find(source);
          if (iter == graph.vid2lvid.end()) {
            lvid_source = graph.vid2lvid.size();
            graph.vid2lvid.insert(std::make_pair(source, lvid_source));
            graph.lvid2record.push_back(vertex_record(source));
          } else {
            lvid_source = iter->second;
          }

          iter = graph.vid2lvid.find(target);
          if (iter == graph.vid2lvid.end()) {
            lvid_target = graph.vid2lvid.size();
            graph.vid2lvid.insert(std::make_pair(target , lvid_target));
            graph.lvid2record.push_back(vertex_record(target));
          } else {
            lvid_target = iter->second;
          }

        local_source_arr.push_back(lvid_source);
        local_target_arr.push_back(lvid_target);
        max_lvid = std::max(std::max(lvid_source, lvid_target), 
            max_lvid);

        ++local_degree_count[vertex_to_proc(source)][source];
        ++local_degree_count[vertex_to_proc(target)][target];
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
      if (max_lvid > 0 && max_lvid >= graph.local_graph.num_vertices()) {
        graph.local_graph.resize(max_lvid + 1);
      }
      graph.local_graph.add_block_edges(local_source_arr, local_target_arr, edata_arr);
      local_graph_lock.unlock();
 
      END_TRACEPOINT(batch_ingress_add_edges);
    } // end of add edges
    

    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      procid_t owning_proc = vertex_to_proc(vid);
      const vertex_buffer_record record(vid, vdata);
      vertex_exchange.send(owning_proc, record);
    } // end of add vertex

    void finalize() { 
      flush(); vertex_exchange.flush();
      rpc.full_barrier();

      BEGIN_TRACEPOINT(batch_ingress_finalize);

      // Check conditions on graph
      if (graph.local_graph.num_vertices() != graph.lvid2record.size()) {
        logstream(LOG_WARNING) << "Finalize check failed. "
                               << "local_graph size: " 
                               << graph.local_graph.num_vertices() 
                               << " not equal to lvid2record size: " 
                               << graph.lvid2record.size()
                               << std::endl;
      }
      ASSERT_EQ(graph.local_graph.num_vertices(), graph.lvid2record.size());

      logstream(LOG_INFO) << "Local graph finalizing: " << std::endl;
      // Finalize local graph
      graph.local_graph.finalize();

      logstream(LOG_INFO) << "Local graph info: " << std::endl
                          << "\t nverts: " << graph.local_graph.num_vertices()
                          << std::endl
                          << "\t nedges: " << graph.local_graph.num_edges()
                          << std::endl;

      // Begin the shuffle phase For all the vertices that this
      // processor has seen determine the "negotiator" and send the
      // negotiator the edge information for that vertex.
      typedef std::vector< std::vector<shuffle_record> > proc2vids_type;
      typedef typename boost::unordered_map<vertex_id_type, lvid_type>::value_type  
        vid2lvid_pair_type;

      proc2vids_type proc2vids(rpc.numprocs());
      foreach(const vid2lvid_pair_type& pair, graph.vid2lvid) {
        const vertex_id_type vid = pair.first;
        const vertex_id_type lvid = pair.second;
        const procid_t negotiator = vertex_to_proc(vid);
        const shuffle_record rec(vid, graph.local_graph.num_in_edges(lvid),
                                 graph.local_graph.num_out_edges(lvid));
        proc2vids[negotiator].push_back(rec);
      }

      // The returned local vertices are the vertices from each
      // machine for which this machine is a negotiator.
      logstream(LOG_INFO) 
        << "Finalize: start exchange shuffle records" << std::endl;
      mpi_tools::all2all(proc2vids, proc2vids);
      logstream(LOG_INFO) 
        << "Finalize: finish exchange shuffle records" << std::endl;

      // Receive any vertex data sent by other machines
      typedef boost::unordered_map<vertex_id_type, vertex_negotiator_record>
        vrec_map_type;
      vrec_map_type vrec_map;
      {
        typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
          vertex_buffer_type;
        vertex_buffer_type vertex_buffer;
        procid_t proc;
        while(vertex_exchange.recv(proc, vertex_buffer)) {
          foreach(const vertex_buffer_record& rec, vertex_buffer) {
            vertex_negotiator_record& negotiator_rec = vrec_map[rec.vid];
            negotiator_rec.vdata = rec.vdata;
          }
        }
      } // end of loop to populate vrecmap

   
      // Update the mirror information for all vertices negotiated by
      // this machine
      logstream(LOG_INFO) 
        << "Finalize: accumulating mirror set for each vertex" << std::endl;
      for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
        foreach(const shuffle_record& shuffle_rec, proc2vids[proc]) {
          vertex_negotiator_record& negotiator_rec = vrec_map[shuffle_rec.vid];
          negotiator_rec.num_in_edges += shuffle_rec.num_in_edges;
          negotiator_rec.num_out_edges += shuffle_rec.num_out_edges;
          negotiator_rec.mirrors.set_bit(proc);
        }
      }


      // Construct the vertex owner assignments and send assignment
      // along with vdata to all the mirrors for each vertex
      logstream(LOG_INFO) << "Constructing and sending vertex assignments" 
                          << std::endl;
      std::vector<size_t> counts(rpc.numprocs());      
      typedef typename vrec_map_type::value_type vrec_pair_type;
      buffered_exchange<vertex_negotiator_record> negotiator_exchange(rpc.dc());
      // Loop over all vertices and the vertex buffer
      foreach(vrec_pair_type& pair, vrec_map) {
        const vertex_id_type vid = pair.first;
        vertex_negotiator_record& negotiator_rec = pair.second;
        negotiator_rec.vid = vid; // update the vid if it has not been set

       // Find the best (least loaded) processor to assign the vertex.
        uint32_t first_mirror = 0; 
        ASSERT_TRUE(negotiator_rec.mirrors.first_bit(first_mirror));
        std::pair<size_t, uint32_t> 
           best_asg(counts[first_mirror], first_mirror);
        foreach(uint32_t proc, negotiator_rec.mirrors) {
            best_asg = std::min(best_asg, std::make_pair(counts[proc], proc));
        }

        negotiator_rec.owner = best_asg.second;
        counts[negotiator_rec.owner]++;
        // Notify all machines of the new assignment
        foreach(uint32_t proc, negotiator_rec.mirrors) {
            negotiator_exchange.send(proc, negotiator_rec);
        }
      } // end of loop over vertex records

      negotiator_exchange.flush();
      logstream(LOG_INFO) << "Recieving vertex data." << std::endl;
      {
        typedef typename buffered_exchange<vertex_negotiator_record>::buffer_type 
          buffer_type;
        buffer_type negotiator_buffer;
        procid_t proc;
        while(negotiator_exchange.recv(proc, negotiator_buffer)) {
          foreach(const vertex_negotiator_record& negotiator_rec, negotiator_buffer) {
            ASSERT_TRUE(graph.vid2lvid.find(negotiator_rec.vid) != 
                        graph.vid2lvid.end());
            const lvid_type lvid = graph.vid2lvid[negotiator_rec.vid];
            ASSERT_LT(lvid, graph.local_graph.num_vertices());
            graph.local_graph.vertex_data(lvid) = negotiator_rec.vdata;
            ASSERT_LT(lvid, graph.lvid2record.size());
            vertex_record& local_record = graph.lvid2record[lvid];
            local_record.owner = negotiator_rec.owner;
            ASSERT_EQ(local_record.num_in_edges, 0); // this should have not been set
            local_record.num_in_edges = negotiator_rec.num_in_edges;
            ASSERT_EQ(local_record.num_out_edges, 0); // this should have not been set
            local_record.num_out_edges = negotiator_rec.num_out_edges;

            ASSERT_TRUE(negotiator_rec.mirrors.begin() != negotiator_rec.mirrors.end());
            local_record._mirrors = negotiator_rec.mirrors;
            local_record._mirrors.clear_bit(negotiator_rec.owner);
          }
        }
      }

      rpc.full_barrier();

      // Count the number of vertices owned locally
      graph.local_own_nverts = 0;
      foreach(const vertex_record& record, graph.lvid2record)
        if(record.owner == rpc.procid()) ++graph.local_own_nverts;

      logstream(LOG_DEBUG) 
        << rpc.procid() << ": local owned vertices: " << graph.local_own_nverts
        << std::endl;

      // Finalize global graph statistics. 
      logstream(LOG_DEBUG)
        << "Finalize: exchange global statistics " << std::endl;

      // Compute edge counts
      std::vector<size_t> swap_counts(rpc.numprocs(), graph.num_local_edges());
      mpi_tools::all2all(swap_counts, swap_counts);
      graph.nedges = 0;
      foreach(size_t count, swap_counts) graph.nedges += count;

      // compute begin edge id
      graph.begin_eid = 0;
      for(size_t i = 0; i < rpc.procid(); ++i) graph.begin_eid += swap_counts[i];

      // Computer vertex count
      swap_counts.assign(rpc.numprocs(), graph.num_local_own_vertices());
      mpi_tools::all2all(swap_counts, swap_counts);
      graph.nverts = 0;
      foreach(size_t count, swap_counts) graph.nverts += count;

      // Computer replicas
      swap_counts.assign(rpc.numprocs(), graph.num_local_vertices());
      mpi_tools::all2all(swap_counts, swap_counts);
      graph.nreplicas = 0;
      foreach(size_t count, swap_counts) graph.nreplicas += count;

      END_TRACEPOINT(batch_ingress_finalize);
    } // end of finalize


  private:

    // HELPER ROUTINES =======================================================>    
    procid_t vertex_to_proc(vertex_id_type vid) const { 
      return vid % rpc.numprocs();
    }    
    
    bool is_local(vertex_id_type vid) const {
      return vertex_to_proc(vid) == rpc.procid();
    }


    void block_add_degree_counts (procid_t pid, vid2degree_type& degree) {
      BEGIN_TRACEPOINT(batch_ingress_update_degree_table);
      typedef typename vid2degree_type::value_type value_pair_type;
      dht_degree_table_lock.readlock();
      foreach (value_pair_type& pair, degree) {
        add_degree_counts(pair.first, pid, pair.second);
      }
      dht_degree_table_lock.unlock();
      END_TRACEPOINT(batch_ingress_update_degree_table);
    }

    // Thread unsafe, used as a subroutine of block add degree counts.
    void add_degree_counts(const vertex_id_type& vid, procid_t pid, 
                           size_t count) {
      size_t idx = vid_to_dht_entry_with_readlock(vid);
      __sync_add_and_fetch(&(dht_degree_table[idx][pid]), count);
    } // end of add degree counts


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

   procid_t edge_to_proc(vertex_id_type src, vertex_id_type dst,
                            std::vector<dht_degree_table_type>& degree_table) {
     size_t src_proc = vertex_to_proc(src);
     size_t dst_proc = vertex_to_proc(dst);
     std::vector<size_t>& src_degree = degree_table[src_proc][src];
     std::vector<size_t>& dst_degree = degree_table[dst_proc][dst];

     procid_t best_proc = -1; 
     double maxscore = 0.0;
     double epsilon = 0.01;
     std::vector<double> proc_score(rpc.numprocs()); 
     
     size_t minedges = *std::min_element(proc_num_edges.begin(), proc_num_edges.end());
     size_t maxedges = *std::max_element(proc_num_edges.begin(), proc_num_edges.end());
     for (size_t i = 0; i < rpc.numprocs(); ++i) {
       size_t sd = src_degree[i]; 
       size_t td = dst_degree[i];

       double bal = (maxedges - proc_num_edges[i])/(epsilon + maxedges - minedges);
       proc_score[i] = bal + ((sd > 0) + (td > 0));
     }
     maxscore = *std::max_element(proc_score.begin(), proc_score.end());

     std::vector<procid_t> top_procs; 
     for (size_t i = 0; i < rpc.numprocs(); ++i)
       if (std::fabs(proc_score[i] - maxscore) < 1e-5)
         top_procs.push_back(i);

      if (top_procs.size() > 1) {
        if (maxscore >= 2) {
          PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_BOTH_TIE, 1)
        } else if (maxscore >= 1) {
          PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_ONE_TIE, 1)
        } else {
          PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_NONE_TIE, 1); 
        }
      } else {
        if (maxscore >= 2) {
          PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_BOTH_UNIQUE, 1);
        } else if (maxscore >= 1) {
          PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_ONE_UNIQUE, 1)
        } else {
          PERMANENT_ACCUMULATE_DIST_EVENT(eventlog, EVENT_EDGE_SEEN_NONE_UNIQUE, 1); 
        }
      }

     // Hash the edge to one of the best procs.
     typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
      boost::hash< edge_pair_type >  hash_function;
      const edge_pair_type edge_pair(std::min(src, dst), 
                                     std::max(src, dst));
      best_proc = top_procs[hash_function(edge_pair) % top_procs.size()];

     ASSERT_LT(best_proc, rpc.numprocs());
     ++src_degree[best_proc];
     ++dst_degree[best_proc];
     ++proc_num_edges[best_proc];
     return best_proc;
   }

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
       procid_t proc = edge_to_proc(e.first, e.second, degree_table);
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
    bool is_full() { return size() >= limit; }

  }; // end of distributed_batch_ingress

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
