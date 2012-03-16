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

#ifndef GRAPHLAB_DISTRIBUTED_LBFS_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_LBFS_INGRESS_HPP


#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
    class distributed_graph;

  template<typename VertexData, typename EdgeData>
  class distributed_local_bfs_ingress: 
    public idistributed_ingress<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData  edge_data_type;

    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;
    /// Type for vertex colors 
    typedef graphlab::vertex_color_type vertex_color_type;

    typedef typename graph_type::lvid_type  lvid_type;
    typedef typename graph_type::vertex_record vertex_record;

    // typedef typename graph_type::SizeType SizeType;
    typedef typename graph_type::mirror_type mirror_type;

    dc_dist_object<distributed_local_bfs_ingress> rpc;
    graph_type& graph;

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

    /** The map from proc_id to num_edges on that proc */
    std::vector<size_t> proc_num_edges;

    /** The map from vertex id to pairs of <pid, local_degree_of_v> */
    typedef typename boost::unordered_map<vertex_id_type, std::vector<size_t> > degree_hash_table_type;
    degree_hash_table_type dht;

    class edge_type {
      vertex_id_type _source;
      vertex_id_type _target;
      size_t _eid;
      public:
        edge_type(vertex_id_type source, vertex_id_type target, size_t eid):
         _source(source), _target(target), _eid(eid)  {}
        vertex_id_type source() const { return _source;}
        vertex_id_type target() const { return _target;}
        size_t edge_id() const { return _eid;}
    }; // end of edge_type

    class bfs_buffer_type {
      public:
        typedef std::pair<vertex_id_type, size_t> adj_entry_type;
        typedef typename boost::unordered_map<vertex_id_type, std::list<adj_entry_type> > adj_list_type;
        typedef typename adj_list_type::iterator adj_iter_type;
      public:
        bfs_buffer_type(){ }
        void add_edge(vertex_id_type source, vertex_id_type target, const EdgeData& edata) {
          edata_list.push_back(edata);
          edge_list[source].push_back(std::make_pair(target, edata_list.size()-1));
        }
        EdgeData& edge_data (edge_type edge) { 
          return edata_list[edge.edge_id()];
        }
        void clear() {
          adj_iter_type iter = edge_list.begin(); 
          while(iter != edge_list.end())
            iter->second.clear();
          edge_list.clear();
          edata_list.clear();
        }
        // Return an list of edge_type ordered by traversing BFS starting from 
        // the input seed vertices.
        std::vector<edge_type> bfs_order(const std::vector<vertex_id_type>& seed){
          std::vector<edge_type> ret;
          std::deque<edge_type> queue; 
          int covered_edge = 0;

          foreach(vertex_id_type startv, seed) {
            adj_iter_type iter = edge_list.find(startv);
            if (iter != edge_list.end()) {
              foreach(adj_entry_type item, iter->second) { 
                queue.push_back(edge_type(startv, item.first, item.second));
                ++covered_edge;
              }
              edge_list[startv].clear();
              edge_list.erase(iter);

              // Fill in the subgraph spanned by the seed vertices.
              while(!queue.empty()) {
                edge_type e = queue.front();
                adj_iter_type iter = edge_list.find(e.target());
                if (iter != edge_list.end()) {
                  foreach(adj_entry_type item, iter->second) { 
                    queue.push_back(edge_type(iter->first, item.first, item.second));
                  }
                  iter->second.clear();
                  edge_list.erase(iter);
                }
                queue.pop_front();
                ret.push_back(e);
                ++covered_edge;
              }
            }
          }

          logstream(LOG_DEBUG) << "Number of covered edges from BFS: " << covered_edge<< std::endl;

          // Fill in the rest of the dangling edges.
          size_t dangling_edge = 0;
          while(!edge_list.empty()) {
            adj_iter_type iter = edge_list.begin();
            foreach(adj_entry_type item, iter->second) {
              ret.push_back(edge_type(iter->first, item.first, item.second));
              ++dangling_edge;
            }
            iter->second.clear();
            edge_list.erase(iter);
          }
          logstream(LOG_DEBUG) << "Number of dangling edges from BFS: " << dangling_edge << std::endl;

          return ret;
        }
      private:
        adj_list_type edge_list;
        std::vector<EdgeData> edata_list;
    }; // end of bfs_buffer
    bfs_buffer_type bfs_buffer;

    vertex_id_type max_vid;
    double seed_percent;
    

    PERMANENT_DECLARE_DIST_EVENT_LOG(eventlog);
    DECLARE_TRACER(ob_ingress_add_edge);
    DECLARE_TRACER(ob_ingress_compute_assignments);
    DECLARE_TRACER(ob_ingress_finalize);

    enum {
      EVENT_EDGE_SEEN_NONE_UNIQUE = 0,
      EVENT_EDGE_SEEN_NONE_TIE = 1,
      EVENT_EDGE_SEEN_ONE_UNIQUE = 2,
      EVENT_EDGE_SEEN_ONE_TIE = 3,
      EVENT_EDGE_SEEN_BOTH_UNIQUE =4,
      EVENT_EDGE_SEEN_BOTH_TIE = 5
    };

  public:
    distributed_local_bfs_ingress(distributed_control& dc, graph_type& graph, double seed_percent = 5) :
      rpc(dc, this), graph(graph), vertex_exchange(dc), edge_exchange(dc),
      proc_num_edges(dc.numprocs()), seed_percent(seed_percent) { 
      rpc.barrier(); 
      ASSERT_GT(seed_percent, 0);
      max_vid = 0;


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


      INITIALIZE_TRACER(ob_ingress_add_edge, "Time spent in add edge");
      INITIALIZE_TRACER(ob_ingress_compute_assignments, "Time spent in compute assignment");
      INITIALIZE_TRACER(ob_ingress_finalize, "Time spent in finalize");
     }

    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      bfs_buffer.add_edge(source, target, edata);
      max_vid = std::max(std::max(source, target), max_vid);
    } // end of add_edge

    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      procid_t owning_proc = vertex_to_proc(vid);
      const vertex_buffer_record record(vid, vdata);
      vertex_exchange.send(owning_proc, record);
    } // end of add vertex

    void flush_bfs_buffer() {
      logstream(LOG_DEBUG) << "Flushing bfs buffer..." << std::endl;
      std::vector<vertex_id_type> seeds = get_seed_vertices(max_vid);
      std::vector<edge_type> edge_list =  bfs_buffer.bfs_order(seeds);

      std::cout << "First 100 bfs edges: " << std::endl;
      for (size_t i = 0; i < 100; ++i) {
        std::cout << "(" << edge_list[i].source() << ", " << edge_list[i].target()
        << ")\t";
      }
        std::cout << std::endl;

      foreach(const edge_type& edge, edge_list) {
        const procid_t owning_proc = edge_to_proc(edge.source(), edge.target());
        const edge_buffer_record record(edge.source(), 
            edge.target(), bfs_buffer.edge_data(edge));
        edge_exchange.send(owning_proc, record);
      }
      bfs_buffer.clear();
    }

    std::vector<vertex_id_type> get_seed_vertices(vertex_id_type max_vid) {
      std::vector<vertex_id_type> ret;
      for (size_t i = 0; i <= max_vid; ++i) {
        if (is_seed(i))
          ret.push_back(i);
      }
      return ret;
    }

    void finalize() { 
      flush_bfs_buffer();
      edge_exchange.flush(); vertex_exchange.flush();
      rpc.full_barrier();
      BEGIN_TRACEPOINT(ob_ingress_finalize);

      // add all the edges to the local graph --------------------------------
      {
        typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
          edge_buffer_type;
        edge_buffer_type edge_buffer;
        procid_t proc;
        while(edge_exchange.recv(proc, edge_buffer)) {
          foreach(const edge_buffer_record& rec, edge_buffer) {
            // Get the source_vlid;
            lvid_type source_lvid(-1);
            if(graph.vid2lvid.find(rec.source) == graph.vid2lvid.end()) {
              source_lvid = graph.vid2lvid.size();
              graph.vid2lvid[rec.source] = source_lvid;
              graph.local_graph.resize(source_lvid + 1);
              graph.lvid2record.push_back(vertex_record(rec.source));
            } else source_lvid = graph.vid2lvid[rec.source];
            // Get the target_lvid;
            lvid_type target_lvid(-1);
            if(graph.vid2lvid.find(rec.target) == graph.vid2lvid.end()) {
              target_lvid = graph.vid2lvid.size();
              graph.vid2lvid[rec.target] = target_lvid;
              graph.local_graph.resize(target_lvid + 1);
              graph.lvid2record.push_back(vertex_record(rec.target));
            } else target_lvid = graph.vid2lvid[rec.target];
            // Add the edge data to the graph
            graph.local_graph.add_edge(source_lvid, target_lvid, rec.edata);          
          } // end of loop over add edges
        } // end for loop over buffers
      }

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
          // negotiator_rec.mirrors.push_back(proc);
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

      END_TRACEPOINT(ob_ingress_finalize);
    } // end of finalize


  private:

    // HELPER ROUTINES =======================================================>    
    procid_t vertex_to_proc(vertex_id_type vid) const { 
      return vid % rpc.numprocs();
    }    

    bool is_local(vertex_id_type vid) const {
      return vertex_to_proc(vid) == rpc.procid();
    }

    int try_edge_hash (vertex_id_type source, vertex_id_type target) const {
      int hashproc = -1;
      bool source_hashed = is_seed(source);
      bool target_hashed = is_seed(target);

      typedef std::pair<vertex_id_type, vertex_id_type> edge_pair_type;
      boost::hash< edge_pair_type >  hash_function;
      const edge_pair_type edge_pair(std::min(source, target), 
                                     std::max(source, target));

      if (source_hashed | target_hashed){
        return hash_function(edge_pair) % rpc.numprocs();
      }
      return hashproc;
    }
    
    bool is_seed (vertex_id_type v) const {
      boost::hash< vertex_id_type>  hash_function;
      return double(hash_function(v) % 100000)/1000.0 <  seed_percent;
    }

    procid_t edge_to_proc(vertex_id_type src, vertex_id_type dst) {
     std::vector<size_t>& src_degree = dht[src];
     std::vector<size_t>& dst_degree = dht[dst];
     if (src_degree.size() == 0)
       src_degree.resize(rpc.numprocs(), 0);
     if (dst_degree.size() == 0)
       dst_degree.resize(rpc.numprocs(), 0);

     procid_t best_proc = -1; 
     double maxscore = 0.0;
     double epsilon = 0.01;
     std::vector<double> proc_score(rpc.numprocs()); 

     int seed_hash = try_edge_hash(src, dst);
     if (seed_hash >= 0) {
       best_proc = (procid_t)seed_hash;
       ++src_degree[best_proc];
       ++dst_degree[best_proc];
       ++proc_num_edges[best_proc];
       return best_proc;
     }

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
     // best_proc = top_procs[std::max(src, dst) % top_procs.size()];
     
     ASSERT_LT(best_proc, rpc.numprocs());
     ++src_degree[best_proc];
     ++dst_degree[best_proc];
     ++proc_num_edges[best_proc];
     return best_proc;
   }
  }; // end of distributed_ob_ingress

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
