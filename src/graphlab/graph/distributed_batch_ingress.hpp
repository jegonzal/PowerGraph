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


#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_graph.hpp>


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


    
    // Helper type used to synchronize the vertex data and assignments
    struct shuffle_record {
      procid_t owner;
      size_t num_in_edges, num_out_edges;
      std::vector<procid_t> mirrors;
      vertex_data_type vdata;
      shuffle_record() : 
        owner(0), num_in_edges(0), num_out_edges(0) { }
      void load(iarchive& arc) { 
        arc >> owner >> num_in_edges >> num_out_edges
            >> mirrors >> vdata;
      } // end of load
      void save(oarchive& arc) const { 
        arc << owner << num_in_edges << num_out_edges
            << mirrors << vdata;
      } // end of save     
    }; // end of shuffle_record

    struct preshuffle_record : public graphlab::IS_POD_TYPE {
      vertex_id_type vid, num_in_edges, num_out_edges;
      preshuffle_record() : vid(0), num_in_edges(0), num_out_edges(0) { }     
    }; // end of preshuffle_record


    dc_dist_object<distributed_batch_ingress> rpc;
    graph_type& graph;

    /** The map from vertex_id to its degree on this proc.*/
    typedef typename boost::unordered_map<vertex_id_type, size_t>  vid2degree_type;
    std::vector<vid2degree_type> local_degree_count;


    /// temporary map for vertexdata
    typedef boost::unordered_map<vertex_id_type, shuffle_record> vid2shuffle_type;
    vid2shuffle_type vid2shuffle;
    mutex vid2shuffle_lock;


    /** The map from proc_id to num_edges on that proc */
    std::vector<size_t> proc_num_edges;

    /** The map from proc_id to num_local_vertices on that proc */
    std::vector<size_t> proc_num_vertices;

    /** The map from proc_id to num_local_own_vertices on that proc */
    std::vector<size_t> proc_num_own_vertices;

    /** The map from vertex id to pairs of <pid, local_degree_of_v> */
    typedef typename boost::unordered_map<vertex_id_type, std::vector<size_t> > 
    dht_degree_table_type;
    dht_degree_table_type dht_degree_table;
    mutex dht_degree_table_lock;


    size_t num_edges;
    size_t limit;
    size_t max_degree;
    std::vector< std::vector<vertex_id_type> > proc_src;
    std::vector< std::vector<vertex_id_type> > proc_dst;
    std::vector< std::vector<EdgeData> > proc_edata;
    
    std::vector<std::pair<vertex_id_type, vertex_id_type> > edgebuf;
    std::vector<EdgeData> databuf;
    std::vector<std::set<vertex_id_type> > query_set;

  

  public:

    distributed_batch_ingress(distributed_control& dc, graph_type& graph) :
      rpc(dc, this), graph(graph),
      local_degree_count(dc.numprocs()),
      num_edges(0), limit(500), max_degree(200), 
      proc_src(dc.numprocs()), proc_dst(dc.numprocs()),
      proc_edata(dc.numprocs()), query_set(dc.numprocs()) { rpc.barrier(); }

    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      ASSERT_LT(edgebuf.size(), limit);
      edgebuf.push_back(std::make_pair(source, target)); 
      databuf.push_back(edata);        
      query_set[vertex_to_init_proc(source)].insert(source);
      query_set[vertex_to_init_proc(target)].insert(target);
      ++num_edges;
      if (is_full()) flush();
    } // end of add_edge


    void add_edges(const std::vector<vertex_id_type>& source_arr, 
                   const std::vector<vertex_id_type>& target_arr, 
                   const std::vector<EdgeData>& edata_arr) {

      // This is a local only method
      ASSERT_TRUE((source_arr.size() == target_arr.size())
                  && (source_arr.size() == edata_arr.size())); 
      if (source_arr.size() == 0) return;

      std::vector<lvid_type> local_source_arr; 
      local_source_arr.reserve(source_arr.size());
      std::vector<lvid_type> local_target_arr;
      local_target_arr.reserve(target_arr.size());

      graph.lvid2record_lock.lock();
      lvid_type max_lvid = 0;
      for (size_t i = 0; i < source_arr.size(); ++i) {
        vertex_id_type source = source_arr[i];
        vertex_id_type target = target_arr[i];

        if (graph.vid2lvid.find(source) == graph.vid2lvid.end()) {
          lvid_type lvid = graph.vid2lvid.size();
          graph.vid2lvid.insert(std::make_pair(source, lvid));
          graph.lvid2record.push_back(vertex_record(source));
        }
        ++local_degree_count[vertex_to_init_proc(source)][source];

        if (graph.vid2lvid.find(target) == graph.vid2lvid.end()) {
          lvid_type lvid = graph.vid2lvid.size();
          graph.vid2lvid.insert(std::make_pair(target, lvid));
          graph.lvid2record.push_back(vertex_record(target));
        }
        ++local_degree_count[vertex_to_init_proc(source)][source];
        local_source_arr.push_back(graph.vid2lvid[source]);
        local_target_arr.push_back(graph.vid2lvid[target]);
        max_lvid = std::max(std::max(graph.vid2lvid[source], graph.vid2lvid[target]), 
                            max_lvid);
      }

      // send out local_degree count;
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
      graph.lvid2record_lock.unlock(); 

      graph.local_graph_lock.lock();
      if (max_lvid > 0 && max_lvid >= graph.local_graph.num_vertices()) {
        graph.local_graph.resize(max_lvid + 1);
      }
      graph.local_graph.add_block_edges(local_source_arr, local_target_arr, edata_arr);
      graph.local_graph_lock.unlock();
    } // end of add edges
    

    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      // determine if the vertex is local
      if(is_local_init(vid)) {
        vid2shuffle_lock.lock();
        vid2shuffle[vid].vdata = vdata;
        vid2shuffle_lock.unlock();
      } else {
        rpc.remote_call(vertex_to_init_proc(vid),
                        &distributed_batch_ingress::add_vertex,
                        vid, vdata);
      }
    } // end of add vertex

    void finalize() { 
      flush();
      rpc.full_barrier();

      //clear dht_degree_table;
      typedef typename dht_degree_table_type::value_type dtable_entry_type;
      foreach(dtable_entry_type& ety, dht_degree_table) {
        std::vector<size_t>().swap(ety.second);
      }
      dht_degree_table_type().swap(dht_degree_table);


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

      // Finalize local graph
      graph.local_graph.finalize();

      // // resize local vid map
      // lvid2vid.resize(vid2record.size());
  
      // For all the vertices that this processor has seen determine the
      // "negotiator" and send that machine the negotiator.
      typedef std::vector< std::vector<preshuffle_record> > proc2vids_type;
      proc2vids_type proc2vids(rpc.numprocs());

      for (size_t lvid = 0; lvid < graph.local_graph.num_vertices(); ++lvid) {
        const vertex_id_type vid = graph.lvid2record[lvid].gvid;
        preshuffle_record pre_rec;
        pre_rec.vid = vid;
        pre_rec.num_in_edges = graph.local_graph.num_in_edges(lvid);
        pre_rec.num_out_edges = graph.local_graph.num_out_edges(lvid);
        proc2vids[vertex_to_init_proc(vid)].push_back(pre_rec);
      }

      // The returned local vertices are the vertices from each
      // machine for which this machine is a negotiator.
      logstream(LOG_DEBUG) 
        << "Finalize: start exchange preshuffle records" << std::endl;

      mpi_tools::all2all(proc2vids, proc2vids);

      logstream(LOG_DEBUG) 
        << "Finalize: finish exchange preshuffle records" << std::endl;


      // Estimate the size of proc2vid
      size_t proc2vid_size = 0;
   
      // Update the vid2shuffle
      logstream(LOG_DEBUG) 
        << "Finalize: update vid 2 shuffle records" << std::endl;
      for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
        foreach(const preshuffle_record& pre_rec, proc2vids[proc]) {
          shuffle_record& shuffle_rec = vid2shuffle[pre_rec.vid];
          shuffle_rec.num_in_edges += pre_rec.num_in_edges;
          shuffle_rec.num_out_edges += pre_rec.num_out_edges;
          shuffle_rec.mirrors.push_back(proc);
        }
        proc2vid_size += proc2vids[proc].capacity() * sizeof(preshuffle_record);
      }

      // Construct the assignments
      logstream(LOG_DEBUG) 
        << "Finalize: constructing assignments" << std::endl;

      std::vector<size_t> counts(rpc.numprocs());
      typedef typename vid2shuffle_type::value_type shuffle_pair_type;
      foreach(shuffle_pair_type& pair, vid2shuffle) {
        shuffle_record& record = pair.second;
        // Find the best (least loaded) processor to assign the vertex.
        std::pair<size_t, procid_t> 
          best_asg(counts[record.mirrors[0]], record.mirrors[0]);
        foreach(procid_t proc, record.mirrors)
          best_asg = std::min(best_asg, std::make_pair(counts[proc], proc));
        record.owner = best_asg.second;
        counts[record.owner]++;
      } // end of loop over 

      // Send the data to all the processors
      logstream(LOG_DEBUG) 
        << "Finalize: send out assignments" << std::endl;

      for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
        typedef std::pair<vertex_id_type, shuffle_record> vid_shuffle_type;
        std::vector<vid_shuffle_type> vertex_assign;
        foreach(const preshuffle_record& pre_rec, proc2vids[proc]) {
          const std::pair<vertex_id_type, shuffle_record> 
            pair(pre_rec.vid, vid2shuffle[pre_rec.vid]);
          vertex_assign.push_back(pair);
        }
        rpc.send_to_nonblocking(proc, vertex_assign);
      }

      // Receive the data from all the processors.  Here we are a little
      // "clever" in that we loop over the vertices we have locally
      // managed and use them to determine how many times to recv_from
      // which machines
      logstream(LOG_DEBUG) 
        << "Finalize: receiving assignments" << std::endl;
      for (size_t i = 0; i < rpc.numprocs();  ++i) {
        std::vector<std::pair<vertex_id_type, shuffle_record> > vertex_assign;
        rpc.recv_from_nonblocking(i, vertex_assign);
        typedef std::pair<vertex_id_type, shuffle_record> vid_shuffle_type;
        foreach (vid_shuffle_type vid_and_rec, vertex_assign) {
          const vertex_id_type& vid = vid_and_rec.first;
          shuffle_record& shuffle_rec = vid_and_rec.second;      


          lvid_type lvid = graph.vid2lvid[vid];
          vertex_record& vrecord = graph.lvid2record[lvid];
          vrecord.mirrors.swap(shuffle_rec.mirrors);
          vrecord.owner = shuffle_rec.owner;

          graph.local_graph.vertex_data(lvid) = shuffle_rec.vdata;
          vrecord.num_in_edges = shuffle_rec.num_in_edges;
          vrecord.num_out_edges = shuffle_rec.num_out_edges;
          if (vrecord.owner == rpc.procid()) 
            ++graph.local_own_nverts;
        }
      }

      rpc.barrier();

      // Finalize global graph statistics. 
      logstream(LOG_DEBUG)
        << "Finalize: exchange global statistics " << std::endl;

      proc_num_edges.assign(rpc.numprocs(), graph.num_local_edges());
      mpi_tools::all2all(proc_num_edges, proc_num_edges);
      graph.begin_eid = 0;
      for (procid_t i = 0; i < rpc.procid(); ++i) {
        graph.begin_eid += proc_num_edges[i];
      }
      graph.nedges = graph.begin_eid;
      for (procid_t i = rpc.procid(); i < rpc.numprocs(); ++i) {
        graph.nedges += proc_num_edges[i];
      }

      proc_num_vertices.assign(rpc.numprocs(), graph.num_local_vertices());
      mpi_tools::all2all(proc_num_vertices, proc_num_vertices);
      for (procid_t i = 0; i < rpc.numprocs(); ++i) {
        graph.nreplica += proc_num_vertices[i];
      }

      proc_num_own_vertices.assign(rpc.numprocs(), graph.num_local_own_vertices());
      mpi_tools::all2all(proc_num_own_vertices, proc_num_own_vertices);
      for (procid_t i = 0; i < rpc.numprocs(); ++i) {
        graph.nverts += proc_num_own_vertices[i];
      }


      // // Receive assignments from coordinators
      // mpi_tools::all2all(vdata_shuffle, vdata_shuffle);
    
      // // Incorporate all the vertex data assigned to each machine
      // for(size_t i = 0; i < vdata_shuffle.size(); ++i) {
      //   foreach(vdata_shuffle_record& shuffle_record, vdata_shuffle[i]) {
      //     ASSERT_TRUE(vid2record.find(shuffle_record.vid) != vid2record.end());
      //     vertex_record& record = vid2record[shuffle_record.vid];
      //     record.owner = shuffle_record.owner; 
      //     record.mirrors.swap(shuffle_record.mirrors);
      //     local_graph.vertex_data(record.lvid) = shuffle_record.vdata;
      //   } // end of loop over vdata
      // } // end of loop over sending machines
   
      // std::cout << "Save debugging information" << std::endl;
      // {
      //   const std::string fname = 
      //     "file_" + boost::lexical_cast<std::string>(rpc.procid());
      //   std::ofstream fout(fname.c_str());
      //   typedef typename vid2record_type::value_type pair_type;
      //   foreach(const pair_type& pair, vid2record) {      
      //     fout << pair.first << '\t' << pair.second.owner << '\t';
      //     std::vector<bool> bitset(rpc.numprocs(), false);
      //     foreach(const procid_t& proc, pair.second.mirrors)
      //       bitset[proc] = true;
      //     for(size_t i = 0; i < bitset.size(); ++i) {
      //       fout << (bitset[i]? '1' : '0') 
      //            << (i+1 < bitset.size()? '\t' : '\n');
      //     }
      //   }
      //   fout.close();
      // }   
    }

  private:

    // HELPER ROUTINES =======================================================>    
    procid_t edge_to_proc(vertex_id_type source, vertex_id_type target) const {
      if(source > target) std::swap(source, target);
      boost::hash< std::pair<vertex_id_type, vertex_id_type> > hash_function;
      return hash_function(std::make_pair(source, target)) % rpc.numprocs();
    }

    bool is_local(vertex_id_type source, vertex_id_type target) const {
      return edge_to_proc(source, target) == rpc.procid();
    }

    procid_t vertex_to_init_proc(vertex_id_type vid) const { 
      return vid % rpc.numprocs();
    }    
    
    bool is_local_init(vertex_id_type vid) const {
      return vertex_to_init_proc(vid) == rpc.procid();
    }


    void block_add_degree_counts (procid_t pid, vid2degree_type degree) {
      typedef typename vid2degree_type::value_type value_pair_type;
      dht_degree_table_lock.lock();
      foreach (value_pair_type& pair, degree) {
        add_degree_counts(pair.first, pid, pair.second);
      }
      dht_degree_table_lock.unlock();
    }

    // Thread unsafe, used as a subroutine of block add degree counts.
    void add_degree_counts(const vertex_id_type& vid, procid_t pid, 
                           size_t count) {
      typedef typename dht_degree_table_type::iterator iterator_type;
      iterator_type iter = dht_degree_table.find(vid);
      if (iter == dht_degree_table.end()) {
        std::vector<size_t>& dtable = dht_degree_table[vid];
        dtable.resize(rpc.numprocs(), 0);
        dtable[pid] += count;
      } else {
        (iter->second)[pid] += count;
      }
    } // end of add degree counts


    dht_degree_table_type 
    block_get_degree_table(const std::set<vertex_id_type>& vid_query) {
      dht_degree_table_type answer;
      typedef typename dht_degree_table_type::iterator iterator_type;
      dht_degree_table_lock.lock();
      foreach (vertex_id_type qvid, vid_query) {
        iterator_type iter = dht_degree_table.find(qvid);
        if (iter == dht_degree_table.end()) {
          answer[qvid] = std::vector<size_t>(rpc.numprocs(), 0);
        } else {
          answer[qvid] = iter->second;
        }
      }
      dht_degree_table_lock.unlock();
      return answer;
    }  // end of block get degree table

   procid_t edge_to_proc(vertex_id_type src, vertex_id_type dst,
                            std::vector<dht_degree_table_type>& degree_table) {
        /// TODO: What is going on here?  Do you really mean to return
        /// immediately
        // Sorry, this is to compare naive partition with greedy partition.
        // return graph_ptr->edge_to_proc(src, dst);
         
        procid_t best_proc = -1; 
        size_t src_proc = vertex_to_init_proc(src);
        size_t dst_proc = vertex_to_init_proc(dst);
        std::vector<size_t>& src_degree = degree_table[src_proc][src];
        std::vector<size_t>& dst_degree = degree_table[dst_proc][dst];

        size_t best_src_proc = -1;
        size_t max_src_degree = 0;
        size_t best_dst_proc = -1;
        size_t max_dst_degree = 0;
          
        for (size_t i = 0; i < rpc.numprocs(); ++i) {
          if (src_degree[i] <= max_degree && src_degree[i] >= max_src_degree)
            { best_src_proc = i; max_src_degree = src_degree[i]; }
          if (dst_degree[i] <= max_degree && dst_degree[i] >= max_dst_degree)
            { best_dst_proc = i; max_dst_degree = dst_degree[i]; }
        }

        // no machine has ever seen this vertex 
        if (max_src_degree == 0) {
          best_src_proc = rand() % rpc.numprocs();
        }
        if (max_dst_degree == 0) {
          best_dst_proc = rand() % rpc.numprocs();
        }

        // std::cout << "best_src_proc: " << best_src_proc
        //   << "\n max_src_degree: " << max_src_degree
        //   << "\n best_dst_proc: " << best_dst_proc
        //   << "\n max_dst_degree: " << max_dst_degree
        //   << std::endl;

        // All procs are full, increase the limit and random assign.
        if (best_src_proc == size_t(-1) && best_dst_proc == size_t(-1)) {
          std::cout << "Double degree limit to " << max_degree << std::endl;
          max_degree*= 2;
          best_proc = rand() % rpc.numprocs();
        } else {
          if (best_src_proc == size_t(-1)) {
            best_proc = best_dst_proc;
          } else if (best_dst_proc == size_t(-1)) {
            best_proc = best_src_proc;
          } else {         
            best_proc = max_src_degree > max_dst_degree ? best_src_proc : best_dst_proc;
          }
        }
        ASSERT_LT(best_proc, rpc.numprocs());
        ++src_degree[best_proc];
        ++dst_degree[best_proc];
        return best_proc;
      }

      void assign_edges() {
        // Get the degree table.
        std::vector<dht_degree_table_type> degree_table(rpc.numprocs());
        for (size_t i = 0; i < rpc.numprocs(); ++i) {
          degree_table[i] = 
            rpc.remote_request(i, 
                              &distributed_batch_ingress::block_get_degree_table,
                              query_set[i]);
          query_set[i].clear();
        }
        
        for (size_t i = 0; i < num_edges; ++i) {
          std::pair<vertex_id_type, vertex_id_type>& e = 
            edgebuf[i];
          procid_t proc = edge_to_proc(e.first, e.second, degree_table);
          ASSERT_LT(proc, proc_src.size());
          proc_src[proc].push_back(e.first);
          proc_dst[proc].push_back(e.second);
          proc_edata[proc].push_back(databuf[i]);
        }
        edgebuf.clear();
        databuf.clear();
      } // end add edge



    // Flush all edges in the buffer.
    void flush() {
      // std::cout << "Flushing edge buffer..." << std::endl;
      assign_edges();
      for (size_t i = 0; i < proc_src.size(); ++i) {
        if (proc_src[i].size() == 0) 
          continue;
        if (i == rpc.procid()) {
          add_edges(proc_src[i], proc_dst[i], proc_edata[i]);
          clear(i);
        } else {
          rpc.remote_call(i, &distributed_batch_ingress::add_edges,
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
