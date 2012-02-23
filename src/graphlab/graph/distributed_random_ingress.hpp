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

#ifndef GRAPHLAB_DISTRIBUTED_RANDOM_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_RANDOM_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_graph.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  template<typename VertexData, typename EdgeData>
  class distributed_random_ingress : 
    public idistributed_ingress<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;

    /// Vertex record
    typedef typename graph_type::lvid_type  lvid_type;
    typedef typename graph_type::vertex_record vertex_record;



    /// The rpc interface for this object
    dc_dist_object<distributed_random_ingress> rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;


    /// Temporar buffers used to store vertex data on ingress
    struct vertex_buffer_record {
      procid_t owner;
      size_t num_in_edges, num_out_edges;
      std::vector<procid_t> mirrors;
      vertex_data_type vdata;
      vertex_buffer_record() : owner(-1), num_in_edges(0), num_out_edges(0) { }
      void load(iarchive& arc) { 
        arc >> owner >> num_in_edges >> num_out_edges
            >> mirrors >> vdata;
      } // end of load
      void save(oarchive& arc) const { 
        arc << owner << num_in_edges << num_out_edges
            << mirrors << vdata;
      } // end of save     
    }; 
    std::vector< boost::unordered_map<vertex_id_type, vertex_buffer_record> > 
    vertex_buffers;
    std::vector< mutex > vertex_buffer_locks;


    /// Temporar buffers used to store edge data on ingress
    struct edge_buffer_record {
      vertex_id_type source, target;
      edge_data_type edata;
      edge_buffer_record(const vertex_id_type& source = vertex_id_type(-1), 
                         const vertex_id_type& target = vertex_id_type(-1), 
                         const edge_data_type& edata = edge_data_type()) :
        source(source), target(target), edata(edata) { }
    };
    std::vector< std::vector<edge_buffer_record> > edge_buffers;
    std::vector< mutex > edge_buffer_locks;

   
    struct shuffle_record : public graphlab::IS_POD_TYPE {
      vertex_id_type vid, num_in_edges, num_out_edges;
      shuffle_record(vertex_id_type vid = 0, vertex_id_type num_in_edges = 0,
                     vertex_id_type num_out_edges = 0) : 
        vid(vid), num_in_edges(num_in_edges), num_out_edges(num_out_edges) { }     
    }; // end of shuffle_record


  public:

    distributed_random_ingress(distributed_control& dc, graph_type& graph) :
      rpc(dc, this), graph(graph),
      vertex_buffers(rpc.numprocs()), vertex_buffer_locks(rpc.numprocs()),
      edge_buffers(rpc.numprocs()), edge_buffer_locks(rpc.numprocs())
    { rpc.barrier(); }

    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      const procid_t owning_proc = edge_to_proc(source, target);
      if(owning_proc == rpc.procid()) { 
        const size_t buffer_ind = random::rand() % edge_buffers.size();
        edge_buffer_locks[buffer_ind].lock();
        edge_buffers[buffer_ind].
          push_back(edge_buffer_record(source, target, edata));
        edge_buffer_locks[buffer_ind].unlock();        
      } else {
        rpc.remote_call(owning_proc, &distributed_random_ingress::add_edge,
                        source, target, edata);
      }
    } // end of add_edge


    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      const procid_t owning_proc = vertex_to_proc(vid);
      if(owning_proc == rpc.procid()) {
        const size_t buffer_ind = vertex_to_buffer(vid);
        vertex_buffer_locks[buffer_ind].lock();
        vertex_buffers[buffer_ind][vid].vdata = vdata;
        vertex_buffer_locks[buffer_ind].unlock();        
      } else {
        rpc.remote_call(owning_proc, &distributed_random_ingress::add_vertex,
                        vid, vdata);
      }
    } // end of add vertex


    void update_vrecord(const vertex_id_type vid, 
                        const vertex_buffer_record& buffer_rec) {
      const lvid_type lvid = graph.vid2lvid[vid];
      vertex_record& local_record = graph.lvid2record[lvid];
      local_record.owner = buffer_rec.owner;
      ASSERT_EQ(local_record.num_in_edges, 0); // this should have not been set
      local_record.num_in_edges = buffer_rec.num_in_edges;
      ASSERT_EQ(local_record.num_out_edges, 0); // this should have not been set
      local_record.num_out_edges = buffer_rec.num_out_edges;
      ASSERT_EQ(local_record.mirrors.size(), 0);
      local_record.mirrors = buffer_rec.mirrors;
    } // end of update_vrecord

    void finalize() { 
      rpc.full_barrier();
      // add all the edges to the local graph --------------------------------
      for(size_t i = 0; i < edge_buffers.size(); ++i) {
        foreach(const edge_buffer_record& rec, edge_buffers[i]) {
          // Get the source_vlid;
          lvid_type source_lvid(-1);
          if(graph.vid2lvid.find(rec.source) == graph.vid2lvid.end()) {
            source_lvid = graph.vid2lvid.size();
            graph.vid2lvid[rec.source] = source_lvid;
            graph.local_graph.resize(source_lvid + 1);
          } else source_lvid = graph.vid2lvid[rec.source];
          // Get the target_lvid;
          lvid_type target_lvid(-1);
          if(graph.vid2lvid.find(rec.target) == graph.vid2lvid.end()) {
            target_lvid = graph.vid2lvid.size();
            graph.vid2lvid[rec.target] = target_lvid;
            graph.local_graph.resize(target_lvid + 1);
          } else target_lvid = graph.vid2lvid[rec.target];
          // Add the edge data to the graph
          graph.local_graph.add_edge(source_lvid, target_lvid, rec.edata);          
        } // end of loop over add edges
        // clear the buffer
        std::vector< edge_buffer_record >().swap(edge_buffers[i]);
      } // end for loop over buffers
      logstream(LOG_INFO) << "Finalizing local graph" << std::endl;

      // Finalize local graph
      graph.local_graph.finalize();

      // Initialize vertex records
      graph.lvid2record.resize(graph.vid2lvid.size());
      typedef typename boost::unordered_map<vertex_id_type, lvid_type>::value_type 
        vid2lvid_pair_type;
      foreach(const vid2lvid_pair_type& pair, graph.vid2lvid) 
        graph.lvid2record[pair.second].gvid = pair.first;      
      // Check conditions on graph
      ASSERT_EQ(graph.local_graph.num_vertices(), graph.lvid2record.size());   
   
      // Begin the shuffle phase For all the vertices that this
      // processor has seen determine the "negotiator" and send the
      // negotiator the edge information for that vertex.
      typedef std::vector< std::vector<shuffle_record> > proc2vids_type;
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

   
      // Update the vid2shuffle
      logstream(LOG_INFO) 
        << "Finalize: update vid 2 shuffle records" << std::endl;
      size_t proc2vid_size = 0;
      for(procid_t proc = 0; proc < rpc.numprocs(); ++proc) {
        foreach(const shuffle_record& rec, proc2vids[proc]) {
          const size_t buffer_index = vertex_to_buffer(rec.vid);
          vertex_buffer_record& vbuffer = vertex_buffers[buffer_index][rec.vid];
          vbuffer.num_in_edges += rec.num_in_edges;
          vbuffer.num_out_edges += rec.num_out_edges;
          vbuffer.mirrors.push_back(proc);
        }
        proc2vid_size += proc2vids[proc].capacity() * sizeof(shuffle_record);
      }
      logstream(LOG_INFO) << "Shuffle record size (bytes): " << proc2vid_size
                          << std::endl;

      // Construct the vertex owner assignments
      logstream(LOG_INFO) << "Finalize: constructing assignments" << std::endl;
      std::vector<size_t> counts(rpc.numprocs());
      typedef boost::unordered_map<vertex_id_type, vertex_buffer_record> 
        vbuffer_map_type;
      typedef typename vbuffer_map_type::value_type vbuffer_pair_type;
      // Loop over all vertices and the vertex buffer
      foreach(vbuffer_map_type& vbuffer, vertex_buffers) {
        foreach(vbuffer_pair_type& pair, vbuffer) {
          const vertex_id_type vid = pair.first;
          vertex_buffer_record& record = pair.second;
          // Find the best (least loaded) processor to assign the vertex.
          std::pair<size_t, procid_t> 
            best_asg(counts[record.mirrors[0]], record.mirrors[0]);
          foreach(procid_t proc, record.mirrors)
            best_asg = std::min(best_asg, std::make_pair(counts[proc], proc));
          record.owner = best_asg.second;
          counts[record.owner]++;
          // Notify all machines of the new assignment
          foreach(procid_t proc, record.mirrors) {
            if(proc == rpc.procid()) update_vrecord(vid, record); 
            else rpc.remote_call(proc, 
                                 &distributed_random_ingress::update_vrecord,
                                 vid, record);
          }
        }
        // destroy the buffer
        vbuffer_map_type().swap(vbuffer);
      } // end of loop over 
    
      logstream(LOG_INFO) << "Finished sending all vertex data." << std::endl;
      rpc.full_barrier();
      logstream(LOG_INFO) << "Finished receiving all vertex data." << std::endl;

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
    procid_t vertex_to_proc(const vertex_id_type vid) const { 
      return vid % rpc.numprocs();
    }        
    bool is_local(const vertex_id_type vid) const {
      return vertex_to_proc(vid) == rpc.procid();
    }

    procid_t edge_to_proc(const vertex_id_type source, 
                          const vertex_id_type target) const {
      typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
      boost::hash< edge_pair_type >  hash_function;
      const edge_pair_type edge_pair(std::min(source, target), 
                                     std::max(source, target));
      return hash_function(edge_pair) % rpc.numprocs();
    }    
    
    bool is_local(const vertex_id_type source,
                  const vertex_id_type target) const {
      return edge_to_proc(source, target) == rpc.procid();
    }

    size_t vertex_to_buffer(const vertex_id_type vid) const {
      return (vid * 131071) % vertex_buffers.size();
    }


  }; // end of distributed_random_ingress

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
