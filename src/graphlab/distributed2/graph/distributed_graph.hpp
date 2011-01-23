#ifndef DISTRIBUTED_GRAPH_HPP
#define DISTRIBUTED_GRAPH_HPP
#include <graphlab/graph/graph.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
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
 * partition + boundary, the machine's <b> replica </b>.
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
 * Each machine has a local representation for its replica of the 
 * graph. Within the local replica, each vertex/edge has a local vertex/edge 
 * ID. This local ID is hidden and abstracted from the user. Implementors 
 * however should keep in mind the following requirements for the local 
 * representation:
 * <ul>
 * <li> Local vertex / edge IDs are unique and sequentially assigned </li>
 * <li> Sorting all vertices/edges in the local replica must
 *      produce the same sequence whether or not we sort by global IDs 
 *      or Local IDs. </li>
 * </ul>
 * 
 * Consistency: 
 * 
 * The object guarantees global sequential consistency of the graph structure 
 * over add_vertex / add_edge operations. Consistency of graph data however,
 * is not managed and must be done manually through the various synchronize()
 * operations. All data reads will be accessed through the local replica if
 * the local replica contains the data. Otherwise, it will be requested from
 * the owner of the data. All data writes will be sent to the owner of the data.
 * The writes may not however, update all replicas unless explicitly requested.
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
  void distributed_graph(distributed_control &dc):rmi(dc, this),
                                                  globalvid2owner(dc, 65536) { 
    init_mmap_store();
    numglobalverts.value = 0;
    numglobaledges.value = 0;
    numlocalverts = 0;
    numlocaledges = 0;
  }
  void distributed_graph(distributed_control &dc, size_t nverts):rmi(dc, this),
                                                                 globalvid2owner(dc, 65536) { 
    init_mmap_store();
    // ok. every machine takes an equal partition of nverts
    numglobalverts.value = nverts;
    numglobaledges.value = 0;
    numlocaledges = 0;

    size_t firstvid = (dc.procid() * nverts) / dc.numprocs();
    size_t lastvid = (((dc.procid() + 1) * nverts) / dc.numprocs()) - 1;
    // generate mappings
    localvdata->resize(lastvid - firstvid + 1);
    // fill the vid mappings
    for (size_t i = firstvid;i <= lastvid; ++i) {
      global2localvid[i] = i - firstvid;
      globalvid2owner.set(i, dc.procid());
      local2globalvid[i - firstvid] = i;
      localvid2owner[i - firstvid] = dc.procid();
    }
    rmi.barrier();
  }
  
  size_t num_vertices() const{
    return numglobalverts;
  }
  
  size_t num_edges() const{
    return numglobaledges;
  }
  
  /**
   * Add vertex picks a random processor and call add_vertex_impl on it
   */
  vertex_id_t add_vertex() {
    vertex_id_t globalvid = get_new_vertex_id();
    // find a machine to store it
    rmi.remote_call(random::rand_int(dc.numprocs() - 1), 
                    &distributed_graph<VertexData, EdgeData>::add_vertex_impl,
                    globalvid);
                    
  }
  
  /**
   * Add vertex picks a random processor and call add_vertex_impl on it
   */
  vertex_id_t add_vertex(const VertexData &vdata) {
    vertex_id_t globalvid = get_new_vertex_id();
    rmi.remote_call(random::rand_int(dc.numprocs() - 1), 
                    &distributed_graph<VertexData, EdgeData>::add_vertex_impl_with_data,
                    globalvid,
                    vdata);
  }
  
  

  edge_id_t add_edge(vertex_id_t source, vertex_id_t target) {
    edge_id_t globaleid = get_new_edge_id();
    
  }
 
 private:
  mutable dc_dist_object<distributed_graph<VertexData, EdgeData> > rmi;
  
  mutex alldatalock;
  // all the mappings requried to move from global to local vid/eids
  // We only store mappings if the vid/eid is stored on this machine
  boost::unordered_map<vertex_id_t, vertex_id_t> global2localvid, local2globalvid;
  boost::unordered_map<vertex_id_t, vertex_id_t> global2localeid, local2globaleid;
  caching_dht<vertex_id_t, procid_t> globalvid2owner;
  boost::unordered_map<vertex_id_t, procid_t> localvid2owner;
  
  atomic<size_t> numglobalverts, numglobaledges;

  
  // --------------  local graph store -----------------------
  
  size_t numlocalverts, numlocaledges;  
  
  /** An internal class describing an edge */
  class edge {
    vertex_id_t _source;
    vertex_id_t _target;
  public:
    edge() : _source(-1), _target(-1) { }
    edge(const edge& other) :
      _source(other.source()), _target(other.target()) { }
    edge(vertex_id_t source, vertex_id_t target) :
      _source(source), _target(target)  { }
    edge(vertex_id_t source, vertex_id_t target) : 
      _source(source), _target(target) {}

    bool operator<(const edge& other) const {
      return (_source < other._source) || 
        (_source == other._source && _target < other._target); 
    }
    
    inline vertex_id_t source() const { return _source; }
    inline vertex_id_t target() const { return _target; }   
    
    void load(iarchive& arc) {
      arc >> _source
          >> _target;
    }
    
    void save(oarchive& arc) const {
      arc << _source
          << _target;
    }
  }; // end of edge
  
  
  
  // vector containing the vertex data. TODO: shift to mmap
  std::vector<VertexData> *localvdata;
  // vector containing the edge data. TODO: shift to mmap
  std::vector<EdgeData> *localedata;
  
  /** A map from local_src_vertex to local incoming edge indices */   
  std::vector< std::vector<edge_id_t> >  in_edges;
  /** A map from local_dest_vertex to local outgoing edge indices */
  std::vector< std::vector<edge_id_t> >  out_edges;

  std::vector<edge> localedges;
  
 
  void init_mmap_store() {
    localvdata = new std::vector<VertexData>;
    localedata = new std::vector<EdgeData>;
  }
  
  vertex_id_t get_new_vertex_id() {
    if (rmi.procid() == 0) return numglobalverts.inc();
    else return rmi.fast_remote_request(0,
                &distributed_graph<VertexData, EdgeData>::get_new_vertex_id);

  }

  edge_id_t get_new_edge_id() {
    if (rmi.procid() == 0) return numglobaledges.inc();
    else return rmi.fast_remote_request(0,
                &distributed_graph<VertexData, EdgeData>::get_new_vertex_id);
  }

  /**
   * Gets a new vertex ID from machine 0, locks the data
   * and updates all the mappings.
   */
  void add_vertex_impl_with_data(size_t globalvid, const VertexData &vdata) {
    
    alldatalock.lock();
    // get the local vertex id. This is just the number of local vertices
    size_t localvid = numlocalverts;
    numlocalverts++;
    // update all mappings
    global2localvid[globalvid] = localvid;
    local2globalvid[localvid] = globalvid;
    globalvid2owner.set(globalvid, rmi.procid());
    localvid2owner[localvid] = rmi.procid();
    // store the data
    localvdata->push_back(vdata);
    alldatalock.unlock();
  }
  
  void add_vertex_impl(size_t globalvid) {
    alldatalock.lock();
    // get the local vertex id. This is just the number of local vertices
    size_t localvid = numlocalverts;
    numlocalverts++;
    // update all mappings
    global2localvid[globalvid] = localvid;
    local2globalvid[localvid] = globalvid;
    globalvid2owner.set(globalvid, rmi.procid());
    localvid2owner[localvid] = rmi.procid();
    // store the data
    localvdata->push_back(VertexData());
    alldatalock.unlock();
  }
  
  /**
   * Gets a new vertex ID from machine 0, locks the data
   * and updates all the mappings.
   */
  void add_edge_impl_with_data(size_t eid, const EdgeData &edata, 
                               size_t globalsrc, size_t globaldest) {
    size_t globaleid = rmi.fast_remote_request(0,
                         &distributed_graph<VertexData, EdgeData>::get_new_edge_id);
    
    alldatalock.lock();
    // get the local vertex id. This is just the number of local vertices
    size_t localeid = numlocaledges;
    numlocaledges++;
    // update all mappings
    global2localeid[globalvid] = localeid;
    local2globaleid[localvid] = globaleid;
    // store the data
    localedata->push_back(edata);
    alldatalock.unlock();
  }
  
  void add_edge_impl(size_t globalsrc, size_t globaldest) {
    size_t globaleid = rmi.fast_remote_request(0,
                         &distributed_graph<VertexData, EdgeData>::get_new_edge_id);
    alldatalock.lock();
    // get the local vertex id. This is just the number of local vertices
    size_t localeid = numlocaledges;
    numlocaledges++;
    // update all mappings
    global2localeid[globalvid] = localeid;
    local2globaleid[localvid] = globaleid;
    // store the data
    localedata->push_back(EdgeData());
    alldatalock.unlock();
  }
};

}
#endif