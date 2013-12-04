

#ifndef GRAPHLAB_DISTRIBUTED_BIPARTITE_AFFINITY_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BIPARTITE_AFFINITY_INGRESS_HPP

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
#include <graphlab/parallel/pthread_tools.hpp>

namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object with data affinity.
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
    typedef typename buffered_exchange<edge_buffer_record>::buffer_type 
        edge_buffer_type;
    typedef typename base_type::vertex_buffer_record vertex_buffer_record;
    typedef typename buffered_exchange<vertex_buffer_record>::buffer_type 
        vertex_buffer_type;

    
    graph_type& graph;
    dc_dist_object<distributed_bipartite_affinity_ingress> bipartite_rpc;
    
        
    typedef typename boost::unordered_map<vertex_id_type, procid_t> 
        master_hash_table_type;
    typedef typename std::pair<vertex_id_type, procid_t> 
        master_pair_type;
    typedef typename buffered_exchange<master_pair_type>::buffer_type 
        master_buffer_type;

    master_hash_table_type mht;
    

    buffered_exchange<edge_buffer_record> bipartite_edge_exchange;
    buffered_exchange<master_pair_type> master_exchange;

    bool source_is_special;
    simple_spinlock bipartite_vertex_lock;
    std::vector<vertex_buffer_record> bipartite_vertexs;
    simple_spinlock bipartite_edge_lock;
    std::vector<edge_buffer_record> bipartite_edges;

  public:
    distributed_bipartite_affinity_ingress(distributed_control& dc, graph_type& graph,const std::string& specialvertex):
      base_type(dc, graph),graph(graph),bipartite_rpc(dc, this),
      bipartite_edge_exchange(dc),master_exchange(dc)
    {
      if(specialvertex=="source")
        source_is_special=true;
      else
        source_is_special=false;            
    } // end of constructor

    ~distributed_bipartite_affinity_ingress() { 
    }

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
      mht[vid]=bipartite_rpc.procid();
      bipartite_vertex_lock.unlock();
    } // end of add vertex

    void finalize() {

      {
        size_t changed_size = bipartite_vertexs.size() + bipartite_edges.size();
        bipartite_rpc.all_reduce(changed_size);
        if (changed_size == 0) {
          logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
          return;
        }
      }
      // directly add the vertices loaded from local file to the local graph
      graph.lvid2record.resize(bipartite_vertexs.size());
      graph.local_graph.resize(bipartite_vertexs.size());
      lvid_type lvid = 0;
      foreach(const vertex_buffer_record& rec,bipartite_vertexs){
        graph.vid2lvid[rec.vid]=lvid;
        graph.local_graph.add_vertex(lvid,rec.vdata);
        vertex_record& vrec = graph.lvid2record[lvid];
        vrec.gvid = rec.vid;
        vrec.owner=bipartite_rpc.procid();
        lvid++;
      }

      // use master_exchange to exchange the mapping table
      for (typename master_hash_table_type::iterator it = mht.begin() ; it != mht.end(); ++it) {
        for (procid_t i = 0; i < bipartite_rpc.numprocs(); ++i) {
          if (i != bipartite_rpc.procid())
            master_exchange.send(i, master_pair_type(it->first, it->second));
        }
      }
      master_exchange.flush();
      master_buffer_type master_buffer;
      procid_t proc = -1;
      while(master_exchange.recv(proc, master_buffer)) {
        foreach(const master_pair_type& pair, master_buffer) {
          mht[pair.first] = pair.second;
        }
      }
      master_exchange.clear();

      // resend all the edges 
      foreach(const edge_buffer_record& rec,bipartite_edges){
        if(source_is_special){
          if(mht.find(rec.source)==mht.end())
            mht[rec.source]=graph_hash::hash_vertex(rec.source) % bipartite_rpc.numprocs();
          const procid_t owning_proc = mht[rec.source] ;
          base_type::edge_exchange.send(owning_proc, rec);
        }
        else{
          if(mht.find(rec.target)==mht.end())
            mht[rec.target]=graph_hash::hash_vertex(rec.target) % bipartite_rpc.numprocs();
          const procid_t owning_proc = mht[rec.target] ;
          base_type::edge_exchange.send(owning_proc, rec);
        }
      }
      bipartite_vertexs.clear();
      bipartite_edges.clear();

      base_type::finalize();

    } // end of finalize

  }; // end of distributed_bipartite_affinity_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
