#ifndef DISTRIBUTED_GRAPH_HPP
#define DISTRIBUTED_GRAPH_HPP
#include <algorithm>
#include <graphlab/distributed2/graph/graph_local_store.hpp>
#include <graphlab/distributed2/graph/atom_index_file.hpp>
#include <graphlab/distributed2/graph/atom_file.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/caching_dht.hpp>
#include <graphlab/util/stl_util.hpp>
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
 * Consistency of graph data, is not managed and must be done manually 
 * through the various synchronize() operations. All data reads will be accessed 
 * through the local fragment if the local fragment contains the data. Otherwise, 
 * it will be requested from the owner of the data. All data writes will be sent 
 * to the owner of the data. The writes may not however, update all fragments 
 * unless explicitly requested.
 * 
 * 
 */
template<typename VertexData, typename EdgeData> 
class distributed_graph {
 
 public:
  distributed_graph(distributed_control &dc, std::string atomidxfile):
                              rmi(dc, this),
                              globalvid2owner(dc, 65536),
                              globaleid2owner(dc, 65536){
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
    construct_local_fragment(atomindex, partitions, rmi.procid());
  }

  /**
   * Returns the number of vertices in the graph.
   */
  size_t num_vertices() const{
      return numglobalverts;
  }

  /**
   * Returns the number of edges in the graph.
   */  
  size_t num_edges() const{
      return numglobaledges;
  }
 
  void print(std::ostream &out) const {
    out << localstore;
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
  dist_graph_impl::graph_local_store<VertexData, EdgeData> localstore;


  /** all the mappings requried to move from global to local vid/eids
   *  We only store mappings if the vid/eid is in the local fragment 
   */
  boost::unordered_map<vertex_id_t, vertex_id_t> global2localvid;
  std::vector<vertex_id_t> local2globalvid;
  boost::unordered_map<vertex_id_t, vertex_id_t> global2localeid;
  std::vector<edge_id_t> local2globaleid;
   
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
  std::vector<procid_t> localvid2owner;
  
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
  
  
  /**
   * From the atoms listed in the atom index file, construct the local store
   * using all the atoms in the current partition.
   */
  void construct_local_fragment(const atom_index_file &atomindex,
                                std::vector<std::vector<size_t> > partitiontoatom,
                                size_t curpartition) {
    // first make a map mapping atoms to machines
    // we will need this later
    std::vector<procid_t> atom2machine;
    for (size_t i = 0 ;i< partitiontoatom.size(); ++i) {
      for (size_t j = 0 ; j < partitiontoatom[i].size(); ++j) {
        if (atom2machine.size() <= partitiontoatom[i][j]) atom2machine.resize(partitiontoatom[i][j] + 1);
        atom2machine[partitiontoatom[i][j]] = i;
      }
    }
    
    
    // the atomfiles for the local fragment
    std::vector<atom_file<VertexData, EdgeData>* > atomfiles;
    // for convenience take a reference to the list of atoms in this partition
    std::vector<size_t>& atoms_in_curpart = partitiontoatom[curpartition];
    
    // create the atom file readers.
    // and load the vid / eid mappings
    atomfiles.resize(atoms_in_curpart.size());
    for (size_t i = 0;i < atoms_in_curpart.size(); ++i) {
      atomfiles[i] = new atom_file<VertexData, EdgeData>;
      atomfiles[i]->input_filename(atomindex.atoms[atoms_in_curpart[i]].protocol,
                                atomindex.atoms[atoms_in_curpart[i]].file);
      atomfiles[i]->load_id_maps();
    }
    
    // Lets first construct the global/local vid/eid mappings by merging
    // the mappings in each atom
    // cat all the globalvids and globaleids into a single big list
    // and sort it
    for (size_t i = 0;i < atomfiles.size(); ++i) {
      std::copy(atomfiles[i]->globalvids().begin(), atomfiles[i]->globalvids().end(),
                std::back_inserter(local2globalvid));
      std::copy(atomfiles[i]->globaleids().begin(), atomfiles[i]->globaleids().end(),
                std::back_inserter(local2globaleid));
    }
    
    // Find only unique occurances of each vertex, by sorting, unique, and resize
    std::sort(local2globalvid.begin(), local2globalvid.end());
    std::vector<vertex_id_t>::iterator uviter = std::unique(local2globalvid.begin(), 
                                                            local2globalvid.end());
    local2globalvid.resize(uviter - local2globalvid.begin());
    
    // do the same thing to each edge    
    std::sort(local2globaleid.begin(), local2globaleid.end());
    std::vector<edge_id_t>::iterator ueiter = std::unique(local2globaleid.begin(), 
                                                            local2globaleid.end());
    local2globaleid.resize(ueiter - local2globaleid.begin());
    
    //construct the reverse maps
    for (size_t i = 0;i < local2globalvid.size(); ++i) global2localvid[local2globalvid[i]] = i;
    for (size_t i = 0;i < local2globaleid.size(); ++i) global2localeid[local2globaleid[i]] = i;
    
    
    // now lets construct the graph structure
    localstore.create_store(local2globalvid.size(), local2globaleid.size(),
                            "vdata." + tostr(curpartition),
                            "edata." + tostr(curpartition));
                            
    std::vector<bool> eidloaded(local2globaleid.size(), false);
    for (size_t i = 0;i < atomfiles.size(); ++i) {
      atomfiles[i]->load_structure();
      for (size_t j = 0;j < atomfiles[i]->edge_src_dest().size(); ++j) {
        // convert from the atom's local eid, to the global eid, then to the fragment localeid
        size_t localeid = global2localeid[atomfiles[i]->globaleids()[j]];
        if (eidloaded[localeid] == false) {
          std::pair<vertex_id_t, vertex_id_t> srcdest = atomfiles[i]->edge_src_dest()[j];
          localstore.add_edge(localeid, srcdest.first, srcdest.second);
          eidloaded[localeid] = true;
        }
      }
      
      // set the color
      for (size_t j = 0; j < atomfiles[i]->vcolor().size(); ++j) {
        // convert from the atom's local vid, to the global vid, then to the fragment localvid
        size_t localvid = global2localvid[atomfiles[i]->globalvids()[j]];
        localstore.color(localvid) = atomfiles[i]->vcolor()[j];
      }
    }
    
    
    // done! structure constructed!
    // now for the data!
    // load atoms one at a time, don't keep more than one atom in memor at any one time
    for (size_t i = 0;i < atomfiles.size(); ++i) {
      atomfiles[i]->load_all();
      for (size_t j = 0; j < atomfiles[i]->vdata().size(); ++j) {
        // convert from the atom's local vid, to the global vid, then to the fragment localvi
        size_t localvid = global2localvid[atomfiles[i]->globalvids()[j]];
        localstore.vertex_data(localvid) = atomfiles[i]->vdata()[j];
      }
      for (size_t j = 0; j < atomfiles[i]->edata().size(); ++j) {
        // convert from the atom's local vid, to the global vid, then to the fragment localvi
        size_t localeid = global2localeid[atomfiles[i]->globaleids()[j]];
        localstore.edge_data(localeid) = atomfiles[i]->edata()[j];
      }
      atomfiles[i]->clear();
      delete atomfiles[i];
    }
    // flush the store
    localstore.flush();
    localstore.compute_minimal_prefetch();
  }
};

template<typename VertexData, typename EdgeData>
std::ostream& operator<<(std::ostream& out,
                           const distributed_graph<VertexData, EdgeData>& graph) {
  graph.print(out);
  return out;
}


}
#endif
