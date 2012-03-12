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

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/idistributed_ingress.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <google/malloc_extension.h>

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


    struct vertex_negotiator_record {
      vertex_id_type vid;
      vertex_id_type num_in_edges, num_out_edges;
      procid_t owner;
      std::vector<procid_t> mirrors;
      vertex_data_type vdata;
      vertex_negotiator_record() : 
        vid(-1), num_in_edges(0), num_out_edges(0), owner(-1) { }
      void load(iarchive& arc) { 
        arc >> vid >> num_in_edges >> num_out_edges >> owner >> mirrors >> vdata;
      }
      void save(oarchive& arc) const { 
        arc << vid << num_in_edges << num_out_edges << owner << mirrors << vdata;
      }
    };

    DECLARE_TRACER(random_ingress);
    DECLARE_TRACER(random_ingress_add_edge);
    DECLARE_TRACER(random_ingress_add_vertex);
    DECLARE_TRACER(random_ingress_recv_edges);
    DECLARE_TRACER(random_ingress_compute_assignments);
    DECLARE_TRACER(random_ingress_finalize);

  public:

    distributed_random_ingress(distributed_control& dc, graph_type& graph) :
      rpc(dc, this), graph(graph), vertex_exchange(dc), edge_exchange(dc) {
      rpc.barrier();
      INITIALIZE_TRACER(random_ingress, "Time spent in random ingress");
      INITIALIZE_TRACER(random_ingress_add_edge, "Time spent in add edge");
      INITIALIZE_TRACER(random_ingress_add_vertex, "Time spent in add vertex" );
      INITIALIZE_TRACER(random_ingress_recv_edges, "Time spent in recv edges");
      INITIALIZE_TRACER(random_ingress_compute_assignments, "Compute Assignment");
      INITIALIZE_TRACER(random_ingress_finalize, "Time spent in finalize");

    } // end of constructor

    ~distributed_random_ingress() { }

    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      BEGIN_TRACEPOINT(random_ingress_add_edge);
      const procid_t owning_proc = edge_to_proc(source, target);
      const edge_buffer_record record(source, target, edata);
      edge_exchange.send(owning_proc, record);
      END_TRACEPOINT(random_ingress_add_edge);
    } // end of add edge

    void add_vertex(vertex_id_type vid, const VertexData& vdata)  { 
      BEGIN_TRACEPOINT(random_ingress_add_vertex);
      const procid_t owning_proc = vertex_to_proc(vid);
      const vertex_buffer_record record(vid, vdata);
      vertex_exchange.send(owning_proc, record);
      END_TRACEPOINT(random_ingress_add_vertex);
    } // end of add vertex


    void finalize() {
      BEGIN_TRACEPOINT(random_ingress_finalize);
      BEGIN_TRACEPOINT(random_ingress_recv_edges);
      edge_exchange.flush(); vertex_exchange.flush();

       size_t value;
       MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
       logstream(LOG_DEBUG) << "Heap Size (entering finalize): " << (double)value/(1024*1024) << "MB" << "\n";
       MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
       logstream(LOG_DEBUG) << "Allocated Size (entering finalize): " << (double)value/(1024*1024) << "MB" << "\n";


      // add all the edges to the local graph --------------------------------
      logstream(LOG_DEBUG) << "Finalize: constructing local graph" << std::endl;
      size_t nedges = edge_exchange.size()+1;
      graph.local_graph.reserve_edge_space(nedges + 1);
      logstream(LOG_INFO) << "Finalize: number of edges adding to the local graph" << nedges << std::endl;
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
        } // end for loop over buffers
      }

      logstream(LOG_INFO) << "Finalizing local graph" << std::endl;
      END_TRACEPOINT(random_ingress_recv_edges);
      edge_exchange.clear();

       MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
       logstream(LOG_DEBUG) << "Heap Size (local graph constructed): " << (double)value/(1024*1024) << "MB" << "\n";
       MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
       logstream(LOG_DEBUG) << "Allocated Size (local graph constructed): " << (double)value/(1024*1024) << "MB" << "\n";



      // Finalize local graph
      graph.local_graph.finalize();
      logstream(LOG_INFO) << "Local graph info: " << std::endl
                          << "\t nverts: " << graph.local_graph.num_vertices()
                          << std::endl
                          << "\t nedges: " << graph.local_graph.num_edges()
                          << std::endl;

       MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
       logstream(LOG_DEBUG) << "Heap Size (local graph finalized): " << (double)value/(1024*1024) << "MB" << "\n";
       MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
       logstream(LOG_DEBUG) << "Allocated Size (local graph finalized): " << (double)value/(1024*1024) << "MB" << "\n";



      BEGIN_TRACEPOINT(random_ingress_compute_assignments);
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

      MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
      logstream(LOG_DEBUG) << "Heap Size (vertex record shuffle): " << (double)value/(1024*1024) << "MB" << "\n";
      MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
      logstream(LOG_DEBUG) << "Allocated Size (vertex record shuffle): " << (double)value/(1024*1024) << "MB" << "\n";


   
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
      END_TRACEPOINT(random_ingress_compute_assignments);

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

      MallocExtension::instance()->GetNumericProperty("generic.heap_size", &value);
      logstream(LOG_DEBUG) << "Heap Size (vertex record finalize): " << (double)value/(1024*1024) << "MB" << "\n";
      MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &value);
      logstream(LOG_DEBUG) << "Allocated Size (vertex record finalize): " << (double)value/(1024*1024) << "MB" << "\n";



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

      const size_t min_ecount = 
        *(std::min_element(swap_counts.begin(), swap_counts.end()));
      const size_t max_ecount = 
        *(std::max_element(swap_counts.begin(), swap_counts.end()));
      std::cout << "BALANCE: " << (double(max_ecount)/min_ecount) << std::endl;

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
      END_TRACEPOINT(random_ingress_finalize);
    } // end of finalize

  private:


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


  }; // end of distributed_random_ingress

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
