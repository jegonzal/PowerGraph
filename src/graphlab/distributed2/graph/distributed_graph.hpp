#ifndef DISTRIBUTED_GRAPH_HPP
#define DISTRIBUTED_GRAPH_HPP
#include <graphlab/distributed2/graph/graph_local_store.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/caching_dht.hpp>
#include <graphlab/logger/assertions.hpp>
namespace graphlab {
  
/**
 * \brief Distributed Graph Implementation.
 * 
 * A fully distributed implementation of the \ref graph object.
 * Vertices are partitioned across machines. Each vertex is owned by a unique
 * machine. Each edge is owned by its destination vertex.
 * Each machine stores all vertex data for vertices within its partition,
 * as well as vertex data for vertices/edges on the boundary of the partition.
 * Each vertex data instance is therefore replicated as many times as the number
 * of distinct machines owning neighbors of the vertex in question.
 * 
 * Formally, where \f$ \Gamma(v)\f$ is the set of neighbors of \f$ v\f$ and \f$o(v)\f$ 
 * is the owner of vertex v, vertex v is replicated
 * \f$ \mbox{DISTINCT} \left( \left\{ o(u), u \in \left\{ v,\Gamma(v) \right\} \right\} \right) \f$
 * times.
 * 
 * Each edge is replicated a maximum of 2 times.
 * 
 * To standardize on terminology, we call the set of vertices and edges owned by 
 * a machine, the machine's <b> partition </b>. We call the set of vertices 
 * and edges adjacent to the partition (but not in the partition), the 
 * <b> boundary </b>. Finally, we will call a machine's local copy of the 
 * partition + boundary, the machine's <b> fragment </b>.
 * 
 * 
 *
 * Vertex / Edge IDs: 
 * 
 * Every vertex/edge in the graph has a uniquely assigned global vertex/edge ID.
 * The task of guaranteeing unique sequential assignment is currently managed by
 * machine 0. 
 * 
 * \note If this is a performance issue in the future, the sequential assignment
 * guarantee could be dropped by having either a post-graph-construction 
 * renumbering scheme, or by building the rest of the components of 
 * GraphLab to not depend on sequential numbering. The user should not 
 * expect the sequential numbering property to be preserved in future versions. 
 * 
 * Each machine has a local representation for its fragment of the 
 * graph. Within the local fragment, each vertex/edge has a local vertex/edge 
 * ID. This local ID is hidden and abstracted from the user. Implementors 
 * however should keep in mind the following requirements for the local 
 * representation:
 * <ul>
 * <li> Local vertex / edge IDs are unique and sequentially assigned </li>
 * <li> Sorting all vertices/edges in the local fragment must
 *      produce the same sequence whether or not we sort by global IDs 
 *      or Local IDs. </li>
 * </ul>
 * 
 * Consistency: 
 * 
 * The object guarantees global sequential consistency of the graph structure 
 * over add_vertex / add_edge operations. Consistency of graph data however,
 * is not managed and must be done manually through the various synchronize()
 * operations. All data reads will be accessed through the local fragment if
 * the local fragment contains the data. Otherwise, it will be requested from
 * the owner of the data. All data writes will be sent to the owner of the data.
 * The writes may not however, update all fragments unless explicitly requested.
 * 
 * Data references:
 * 
 * Now, unlike in the \ref graph datastructure, vertex_data(v) and edge_data(v)
 * cannot return a direct reference to the data since it could be on a remote
 * machine. It must therefore return a proxy object 
 * much like how std::vector<bool> works. This means that the following code 
 * will no longer compile:
 * 
 * \code
 * vertex_data &v = graph.vertex_data(v);
 * dosomething(v);
 * v = somethingelse;
 * \endcode
 * 
 * instead, the user should call graph.vertex(v) directly to both read and write.
  * \code
 * dosomething(graph.vertex_data(v)); // dosomething must not take by reference!
 * graph.vertex_data(v) = somethingelse;
 * \endcode
 * 
 */
template<typename VertexData, typename EdgeData> 
class distributed_graph {
 
 public:
  void distributed_graph(distributed_control &dc, std::string atomidxfile):
                              rmi(dc, this),
                              globalvid2owner(dc, 65536) {
    // read the atom index.
    atom_index_file atomindex = read_atom_index(atomidxfile);
    // store the graph size
    numglobalverts = atomindex.nverts;
    numglobaledges = atomindex.nedges;
    // machine 0 partitions it
    std::vector<std::vector<size_t> > partitions;
    if (dc.procid() == 0) {
      partitions = partition_atoms(atomindex, dc.numprocs());
    }
    dc.services().broadcast(partitions, dc.procid() == 0);
    construct_local_subgraph(atomindex, partitions, rmi.procid());
  }

  /**
   * Returns the number of vertices in the graph.
   * If machine 0, we can return this immediately since machine 0
   * manages the global vid allocation. Otherwise, a remote call is necessary.
   */
  size_t num_vertices() const{
    if (rmi.dc().procid() == 0) {
      return numglobalverts.value;
    }
    else {
      return rmi.fast_remote_request(0,
                 distributed_graph<VertexData,EdgeData>::num_vertices);                    
    }
  }

  /**
   * Returns the number of edges in the graph.
   * If machine 0, we can return this immediately since machine 0
   * manages the global eid allocation. Otherwise, a remote call is necessary.
   */  
  size_t num_edges() const{
    if (rmi.dc().procid() == 0) {
      return numglobaledges.value;
    }
    else {
      return rmi.fast_remote_request(0,
                 distributed_graph<VertexData,EdgeData>::num_edges);                    
    }
  }
 
 private:
  /// RMI object
  mutable dc_dist_object<distributed_graph<VertexData, EdgeData> > rmi;
  
  /** Protects structural modifications of the graph.
   * Modifications to the data store and to the local<->global mappings
   * must lock this.
   */
  mutex alldatalock;
  
  /// stores the local fragment of the graph
  graph_local_store<VertexData, EdgeData> localstore;


  /** all the mappings requried to move from global to local vid/eids
   *  We only store mappings if the vid/eid is in the local fragment 
   */
  boost::unordered_map<vertex_id_t, vertex_id_t> global2localvid, local2globalvid;
  boost::unordered_map<vertex_id_t, vertex_id_t> global2localeid, local2globaleid;
   
  /** To avoid requiring O(V) storage on each maching, the 
   * global_vid -> owner mapping cannot be stored in its entirely locally
   * instead, we store it in a DHT. \see globaleid2owner
   */
  caching_dht<vertex_id_t, procid_t> globalvid2owner;
  
  /** To avoid requiring O(E) storage on each maching, the 
   * global_eid -> owner mapping cannot be stored in its entirely locally
   * instead, we store it in a DHT \see globalvid2owner
   */
  caching_dht<vertex_id_t, procid_t> globaleid2owner;
  
  /** This provides a fast mapping from the local vids in the fragment
   * to its owner. Since this operation is quite frequently needed.
   */
  boost::unordered_map<vertex_id_t, procid_t> localvid2owner;
  
  /**
   * The number of vertices and edges in the entire graph so far.
   * Currently only consistent on machine 0 since machine 0 manages 
   * the allocation of global VIDs and local VIDs.
   */
  size_t numglobalverts, numglobaledges;

  
  
  /**
   * Returns true if the global vid is in the local fragment
   * This is not synchronized. Caller must lock if there is a risk
   * of the structure changing while this check is performed.
   */
  bool global_vid_in_local_fragment(vertex_id_t globalvid) {
    // easiest way to check is to see if it is in the global2localvid mapping
    return global2localvid.find(globalvid) != global2localvid.end();
  }
  
  /**
   * Returns true if the global eid is in the local fragment
   * This is not synchronized. Caller must lock if there is a risk
   * of the structure changing while this check is performed.
   */
  bool global_eid_in_local_fragment(edge_id_t globaleid) {
    // easiest way to check is to see if it is in the global2localvid mapping
    return global2localeid.find(globaleid) != global2localeid.end();
  }
  
  
  void construct_local_subgraph(atom_index_file atomindex,
                                std::vector<std::vector<size_t> > partitiontoatom,
                                size_t curpartition) {
    
    // merge all the atoms in partitiontoatom[curpartition]
    
    for (size_t i = 0;i < partitiontoatom[curpartition].size(); ++i) {
      
    }
  }
};

}
#endif
