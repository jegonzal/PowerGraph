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

#ifndef GRAPHLAB_DISTRIBUTED_BATCH_INGRESS2_HPP
#define GRAPHLAB_DISTRIBUTED_BATCH_INGRESS2_HPP


#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
    class distributed_graph;

  template<typename VertexData, typename EdgeData>
  class distributed_batch_ingress2 : 
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

    dc_dist_object<distributed_batch_ingress2> rpc;
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
      std::vector<procid_t> mirrors;
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

    // Local minibatch buffer 
    size_t num_edges;
    size_t limit;
    std::vector< std::vector<vertex_id_type> > proc_src;
    std::vector< std::vector<vertex_id_type> > proc_dst;
    std::vector< std::vector<EdgeData> > proc_edata;
    std::vector<std::pair<vertex_id_type, vertex_id_type> > edgesend;
    std::vector<EdgeData> edatasend;
    std::vector<std::set<vertex_id_type> > query_set;

    /** The map from proc_id to num_edges on that proc */
    std::vector<size_t> proc_num_edges;


    /** The map from vertex id to pairs of <pid, local_degree_of_v> */
    typedef typename std::vector< std::vector<size_t> > degree_hash_table_type;
    degree_hash_table_type dht;
    mutex dht_degree_table_lock;


    DECLARE_TRACER(batch_ingress_add_edge);
    DECLARE_TRACER(batch_ingress_flush);
    DECLARE_TRACER(batch_ingress_add_edges);
    DECLARE_TRACER(batch_ingress_compute_assignments);
    DECLARE_TRACER(batch_ingress_update_degree_table);
    DECLARE_TRACER(batch_ingress_finalize);


  public:
    distributed_batch_ingress2(distributed_control& dc, graph_type& graph) :
      rpc(dc, this), graph(graph), vertex_exchange(dc),
      num_edges(0), limit(4000), 
      proc_src(dc.numprocs()), proc_dst(dc.numprocs()),
      proc_edata(dc.numprocs()), query_set(dc.numprocs()),
     proc_num_edges(dc.numprocs()), dht(1572869) { 
       rpc.barrier(); 

       for (size_t i = 0; i < dht.size(); ++i) {
         dht[i].resize(rpc.numprocs(), 0);
       }

      INITIALIZE_TRACER(batch_ingress_add_edge, "Time spent in add edge");
      INITIALIZE_TRACER(batch_ingress_flush, "Time spent in flush");
      INITIALIZE_TRACER(batch_ingress_add_edges, "Time spent in add block edges" );
      INITIALIZE_TRACER(batch_ingress_compute_assignments, "Time spent in compute assignment");
      INITIALIZE_TRACER(batch_ingress_update_degree_table, "Time spent in update degree table");
      INITIALIZE_TRACER(batch_ingress_finalize, "Time spent in finalize");
     }

    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      BEGIN_TRACEPOINT(batch_ingress_add_edge);
      ASSERT_LT(edgesend.size(), limit);
      edgesend.push_back(std::make_pair(source, target)); 
      edatasend.push_back(edata);        
      query_set[vertex_to_proc(source)].insert(source);
      query_set[vertex_to_proc(target)].insert(target);
      ++num_edges;
      END_TRACEPOINT(batch_ingress_add_edge);

      BEGIN_TRACEPOINT(batch_ingress_flush);
      if (is_full()) flush();
      END_TRACEPOINT(batch_ingress_flush);
    } // end of add_edge



    // This is a local only method
    // Copy the rpc sent edge into a local recv buffer.
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
      boost::unordered_map<vertex_id_type, size_t> local_degree_count;
     
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

        ++local_degree_count[source];
        ++local_degree_count[target];
      }
      lvid2record_lock.unlock();

      // Add edges to local graph
      local_graph_lock.lock();
      if (max_lvid > 0 && max_lvid >= graph.local_graph.num_vertices()) {
        graph.local_graph.resize(max_lvid + 1);
      }
      graph.local_graph.add_block_edges(local_source_arr, local_target_arr, edata_arr);
      local_graph_lock.unlock();

      
      // Send out local_degree count;
      for (size_t i = 0; i < rpc.numprocs(); ++i) {
        if (i != rpc.procid()) {
          rpc.remote_call(i, 
                          &distributed_batch_ingress2::block_add_degree_counts, 
                          rpc.procid(),
                          local_degree_count);
        } else {
          block_add_degree_counts(rpc.procid(), local_degree_count);
        }
      }
      END_TRACEPOINT(batch_ingress_add_edges);
    } // end of add edges
    

    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      procid_t owning_proc = vertex_to_proc(vid);
      const vertex_buffer_record record(vid, vdata);
      vertex_exchange.send(owning_proc, record);
    } // end of add vertex

    void finalize() { 
      flush();
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
          negotiator_rec.mirrors.push_back(proc);
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
        std::pair<size_t, procid_t> 
          best_asg(counts[negotiator_rec.mirrors[0]], negotiator_rec.mirrors[0]);
        foreach(procid_t proc, negotiator_rec.mirrors)
          best_asg = std::min(best_asg, std::make_pair(counts[proc], proc));
        negotiator_rec.owner = best_asg.second;
        counts[negotiator_rec.owner]++;
        // Notify all machines of the new assignment
        foreach(procid_t dest, negotiator_rec.mirrors) 
          negotiator_exchange.send(dest, negotiator_rec);
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
            ASSERT_GT(negotiator_rec.mirrors.size(), 0);
            local_record._mirrors.reserve(negotiator_rec.mirrors.size()-1);
            ASSERT_EQ(local_record._mirrors.size(), 0);
            // copy the mirrors but drop the owner
            for(size_t i = 0; i < negotiator_rec.mirrors.size(); ++i) {
              if(negotiator_rec.mirrors[i] != negotiator_rec.owner) 
                local_record._mirrors.push_back(negotiator_rec.mirrors[i]);
            }
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

    size_t vertex_to_dht_hash(vertex_id_type vid) const {
      return vid % dht.size();
    }


    void block_add_degree_counts (procid_t pid, vid2degree_type& degree) {
      BEGIN_TRACEPOINT(batch_ingress_update_degree_table);
      typedef typename vid2degree_type::value_type value_pair_type;
      dht_degree_table_lock.lock();
      foreach (value_pair_type& pair, degree) {
        size_t hashidx = vertex_to_dht_hash (pair.first);
        dht[hashidx][pid] += pair.second;
      }
      dht_degree_table_lock.unlock();
      END_TRACEPOINT(batch_ingress_update_degree_table);
    }

    procid_t edge_to_proc(vertex_id_type src, vertex_id_type dst) {
     std::vector<size_t>& src_degree = dht[vertex_to_dht_hash(src)];
     std::vector<size_t>& dst_degree = dht[vertex_to_dht_hash(dst)];

     procid_t best_proc = -1; 
     double maxscore = 0.0;
     double epsilon = 1e-5;
     std::vector<double> proc_score(rpc.numprocs()); 
     for (size_t i = 0; i < rpc.numprocs(); ++i) {
       size_t sd = src_degree[i];
       size_t td = src_degree[i];

       size_t minedges = *std::min_element(proc_num_edges.begin(), proc_num_edges.end());
       size_t maxedges = *std::max_element(proc_num_edges.begin(), proc_num_edges.end());
       double bal = (maxedges - proc_num_edges[i])/(epsilon + maxedges - minedges);
       proc_score[i] = bal;
       if (!(sd || td)) { // proc hasn't seen either src or dst
         proc_score[i] += 0; 
       } else if (!(sd && td)) { // proc has seen one but not the other
         proc_score[i] += 1; 
       } else {
         proc_score[i] += 2;
       }
       if (proc_score[i] > maxscore) {
         maxscore = proc_score[i];
       }
     }

     std::vector<procid_t> top_procs; 
     for (size_t i = 0; i < rpc.numprocs(); ++i)
       if (std::fabs(proc_score[i] - maxscore) < epsilon)
         top_procs.push_back(i);

     // Hash the edge to one of the best procs.
     best_proc = top_procs[std::max(src, dst) % top_procs.size()];
     ASSERT_LT(best_proc, rpc.numprocs());

     ++src_degree[best_proc];
     ++dst_degree[best_proc];
     ++proc_num_edges[best_proc];
     return best_proc;
   }

   void assign_edges() {
     ASSERT_EQ(num_edges, edgesend.size());
     if (num_edges == 0) return;

     for (size_t i = 0; i < num_edges; ++i) {
       std::pair<vertex_id_type, vertex_id_type>& e = 
         edgesend[i];

       BEGIN_TRACEPOINT(batch_ingress_compute_assignments);
       procid_t proc = edge_to_proc(e.first, e.second);
       END_TRACEPOINT(batch_ingress_compute_assignments);

       ASSERT_LT(proc, proc_src.size());
       proc_src[proc].push_back(e.first);
       proc_dst[proc].push_back(e.second);
       proc_edata[proc].push_back(edatasend[i]);
     }
     edgesend.clear();
     edatasend.clear();
   } // end add edge

    // Flush all edges in the buffer.
    void flush() {
      assign_edges();
      for (size_t i = 0; i < proc_src.size(); ++i) {
        if (proc_src[i].size() == 0) 
          continue;
        if (i == rpc.procid()) {
          add_edges(proc_src[i], proc_dst[i], proc_edata[i]);
          clear(i);
        } else {
          rpc.remote_call(i, &distributed_batch_ingress2::add_edges,
              proc_src[i], proc_dst[i], proc_edata[i]);
          clear(i);
        } // end if
      } // end for
    } // end flush

    size_t size() { return num_edges; }
    bool is_full() { return size() >= limit; }

    // Clear the nth slot
    void clear(size_t n) {
      ASSERT_LT(n, proc_src.size());
      num_edges -= proc_src[n].size();
      proc_src[n].clear();
      proc_dst[n].clear();
      proc_edata[n].clear();
    }
    void clear_all() {
      for (size_t i = 0; i < proc_src.size(); ++i)
        clear(i);
      ASSERT_EQ(num_edges, 0);
    }
  }; // end of distributed_batch_ingress

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
