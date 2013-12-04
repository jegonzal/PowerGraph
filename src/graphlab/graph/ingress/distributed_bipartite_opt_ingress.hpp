

#ifndef GRAPHLAB_DISTRIBUTED_BIPARTITE_OPT_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_BIPARTITE_OPT_INGRESS_HPP

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

namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object benefit for bipartite graph.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_bipartite_opt_ingress : 
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
    dc_dist_object<distributed_bipartite_opt_ingress> bipartite_rpc;

    buffered_exchange<vertex_buffer_record> bipartite_vertex_exchange;
    buffered_exchange<edge_buffer_record> bipartite_edge_exchange;
    
    
    std::vector<edge_buffer_record> bipartite_edges;

    bool source_is_special;



  public:
    distributed_bipartite_opt_ingress(distributed_control& dc, graph_type& graph,const std::string& specialvertex) :
      base_type(dc, graph),graph(graph),bipartite_rpc(dc, this),
      bipartite_vertex_exchange(dc),bipartite_edge_exchange(dc)
    {

      if(specialvertex=="source")
        source_is_special=true;
      else
        source_is_special=false;         
      
    } // end of constructor

    ~distributed_bipartite_opt_ingress() { 
      
    }

    /** Add an edge to the ingress object to the favorite subset. */
    void add_edge(vertex_id_type source, vertex_id_type target,
                  const EdgeData& edata) {
      const edge_buffer_record record(source, target, edata);
      if(source_is_special){  
        const procid_t owning_proc = 
          graph_hash::hash_vertex(source) % bipartite_rpc.numprocs();
        bipartite_edge_exchange.send(owning_proc, record);
      }
      else{
        const procid_t owning_proc = 
          graph_hash::hash_vertex(target) % bipartite_rpc.numprocs();
        bipartite_edge_exchange.send(owning_proc, record);
      }
    } // end of add edge
    
    void add_vertex(vertex_id_type vid, const VertexData& vdata) { 
      const vertex_buffer_record record(vid, vdata);
      const procid_t owning_proc = 
        graph_hash::hash_vertex(vid) % bipartite_rpc.numprocs();        
      bipartite_vertex_exchange.send(owning_proc, record);
    } // end of add vertex

    void finalize() {
      

      edge_buffer_type edge_buffer;
      
      procid_t proc;

      bipartite_edge_exchange.flush();
      bipartite_vertex_exchange.flush();
      boost::unordered_map<vertex_id_type,std::vector<int> > count_map;
      boost::unordered_map<vertex_id_type,procid_t> owner_map;
      std::set<vertex_id_type> master_set;

      std::vector<int> proc_num_edges;

      { //skipping finalization when no changes hapened
        size_t changed_size = bipartite_edge_exchange.size() + bipartite_vertex_exchange.size();
        bipartite_rpc.all_reduce(changed_size);
        if (changed_size == 0) {
          logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
          return;
        }
      }

      // recv the edges 
      // count_map is used to calculate the distribution of a favorite vertice 's neighbors  
      // proc_num_edges contains the number of edges in all machines sent from this machine.
      proc_num_edges.resize(bipartite_rpc.numprocs());
      proc = -1;
      while(bipartite_edge_exchange.recv(proc, edge_buffer)) {
        foreach(const edge_buffer_record& rec, edge_buffer) {
          bipartite_edges.push_back(rec);
          if(source_is_special){
            if(count_map.find(rec.source)==count_map.end()){
              count_map[rec.source].resize(bipartite_rpc.numprocs());
            }
            const procid_t target_proc = 
              graph_hash::hash_vertex(rec.target) % bipartite_rpc.numprocs();
            count_map[rec.source][target_proc]+=1;
          } else {
            if(count_map.find(rec.target)==count_map.end()){
              count_map[rec.target].resize(bipartite_rpc.numprocs());
            }
            const procid_t source_proc = 
              graph_hash::hash_vertex(rec.source) % bipartite_rpc.numprocs();
            count_map[rec.target][source_proc]+=1;
          }
        }
      }
      bipartite_edge_exchange.clear();
      
      // a loop to calculate where should each favorite vertice go.
      // and tell that machine with  vid_buffer . 
      proc = -1;
      buffered_exchange<vertex_id_type> vid_buffer(bipartite_rpc.dc());
      for (typename boost::unordered_map<vertex_id_type,std::vector<int> >::iterator it = count_map.begin();
               it != count_map.end(); ++it){
        procid_t maxproc=bipartite_rpc.procid();
        double maxscore =(it->second)[maxproc]-sqrt(1.0*proc_num_edges[maxproc]);
        for(size_t i=0; i < bipartite_rpc.numprocs();i++){
          double current_score=(it->second)[i]-sqrt(1.0*proc_num_edges[i]);
          if(current_score > maxscore){
            maxscore =current_score;
            maxproc  =i;
          }
        }
        for(size_t i=0; i < bipartite_rpc.numprocs();i++){
          proc_num_edges[maxproc]+=(it->second)[i];
        }
        owner_map[it->first]=maxproc;
        vid_buffer.send(maxproc, it->first);
      }

      // find all the favorite vertices this machine own
      vid_buffer.flush();
      {
        typename buffered_exchange<vertex_id_type>::buffer_type buffer;
        procid_t recvid = -1;
        while(vid_buffer.recv(recvid, buffer)) {
          foreach(const vertex_id_type vid, buffer) {
            if(source_is_special)
              master_set.insert(vid);
            else
              master_set.insert(vid);
          }
        }
      }
      vid_buffer.clear();

      // send the all the edges out using  base_type::edge_exchange
      for (size_t i = 0; i < bipartite_edges.size(); i++) {
        edge_buffer_record& rec = bipartite_edges[i];
        procid_t owner_proc ;
        if(source_is_special)
          owner_proc = owner_map[rec.source];
        else
          owner_proc = owner_map[rec.target];
        base_type::edge_exchange.send(owner_proc,rec);   
      }
      
      // add the favorite vertices to the local graph.
      graph.lvid2record.resize(master_set.size());
      graph.local_graph.resize(master_set.size());
      lvid_type lvid = 0;
      foreach(const vertex_id_type& vid,master_set){
        graph.vid2lvid[vid]=lvid;
        vertex_record& vrec = graph.lvid2record[lvid];
        vrec.gvid = vid;
        vrec.owner=bipartite_rpc.procid();
        lvid++;
      }

      { // resend all vertex data .
        vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
        while(bipartite_vertex_exchange.recv(sending_proc, vertex_buffer)) {
          foreach(const vertex_buffer_record& rec, vertex_buffer) {
            if(owner_map.find(rec.vid) !=owner_map.end()){
              base_type::vertex_exchange.send(owner_map[rec.vid],rec);
            }
            else{
              base_type::vertex_exchange.send(bipartite_rpc.procid(),rec);
            }
          }
        }
        bipartite_vertex_exchange.clear();
      } 

      // reuse the finalize of base_type
      base_type::finalize();
    } // end of finalize

  }; // end of distributed_bipartite_opt_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
