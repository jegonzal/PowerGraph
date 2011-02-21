#ifndef DISTRIBUTED_CLONED_GRAPH_HPP
#define DISTRIBUTED_CLONED_GRAPH_HPP
#include <cassert>
#include <cmath>
#include <vector>

#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/util/generics/blob.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed/distributed_control.hpp>




namespace graphlab { 
  

  template<typename VertexData, typename EdgeData> class cloned_graph;
  typedef cloned_graph<blob, blob> blob_cloned_graph;

  
  /**
   * \class cloned_graph
   * This implements a simple distributed graph datastructure which clones
   * the graph on all machines when distributed_partition() and distribute() are 
   * called.
   * 
   * The cloned_graph does not support structural mutation once the graph has 
   * been distributed. Data mutation is managed in an unsynchronized
   * fashion and must be manually triggered through the update_vertex(),
   * and the update_edge() function. This will broadcast the changes
   * to all machines.
   *
   * The standard use case is:
   *  - Create cloned_graph
   *  - create structure and/or initial data
   *  - call distributed_partition() and distribute()
   *  - perform only data mutation. 
   
   * Note that distributed_partition() and distribute() must be called 
   * simultaneously by all processors.
   */
  template<typename VertexData, typename EdgeData>
  class cloned_graph {
  public:

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
      
    typedef edge_list edge_list_type;

  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    cloned_graph():distributed(false), constant_edges(false) {    }


    cloned_graph(size_t nverts) : 
      mgraph(nverts), distributed(false), constant_edges(false) { }

    cloned_graph(const graph<VertexData, EdgeData> &g): distributed(false), constant_edges(false) { 
      mgraph = g;
      
    }


    // DISTRIBUTION FUNCTIONS ================================================
    /** Generates a graph partitioning. All processes should call this 
    simultaneously.  */
    void distributed_partition(distributed_control &dc, 
                               partition_method::partition_method_enum p, 
                               const size_t overpartition) {
      ASSERT_FALSE(distributed);
      finalize();
      // only node 0 does partitioning for 
      if (dc.procid() == 0) {
        size_t numparts = overpartition * dc.numprocs();
        std::vector<uint32_t> parts;
        partition(p, numparts, parts);

        // prepare the owner2vertex / vertex2owner structures
        vertex2owner.resize(num_vertices());
        owner2vertex.resize(dc.numprocs());
        for (vertex_id_t i = 0;i < num_vertices(); ++i) {
          vertex2owner[i] = parts[i] % dc.numprocs();
          owner2vertex[parts[i] % dc.numprocs()].push_back(i);
        }
        for (procid_t i = 0;i < dc.numprocs(); ++i) {
          logstream(LOG_INFO) << "Part " << i << ": "
                              << owner2vertex[i].size() << " vertices" << std::endl;
        }
      }
      // done
    }

    edge_id_t global_to_local_eid(edge_id_t e) {
      return e;
    }

    edge_id_t local_to_global_eid(edge_id_t e) {
      return e;
    }
    size_t num_local_edges() const{
      return num_edges();
    }
    /** Generates a graph partitioning. All processes should call this
    simultaneously.  */
    void set_partition(distributed_control &dc,
                       const std::vector<procid_t> &vtx2owner) {
      ASSERT_FALSE(distributed);
      finalize();
      // only node 0 does partitioning for
      if (dc.procid() == 0) {
        ASSERT_EQ(vtx2owner.size(), num_vertices());
        vertex2owner = vtx2owner;
        owner2vertex.resize(dc.numprocs());
        for (vertex_id_t i = 0;i < num_vertices(); ++i) {
          owner2vertex[vertex2owner[i] % dc.numprocs()].push_back(i);
        }
      }
      // done
    }
   /** All processes must call "distribute()" simultaneously.
    * This essentially does a "broadcast" of the graph to all machines.
    * This function blocks until the graph has been transmitted to all 
    * machines. 
    */
    void distribute(distributed_control &dc) {
      logger(LOG_INFO, "%d: entering distribute()", dc.procid());
      ASSERT_FALSE(distributed);
      // remember the controller
      dcontrol = &dc;
      myprocid = dc.procid();
      // set the receive target of the handler
      receive_target = this;
      dc.barrier();
      // if I am process 0, send the graph
      if (dc.procid() == 0) {
        for (size_t i = 1;i < dc.numprocs(); ++i) {
          dc.remote_callxs(i, cloned_graph<VertexData, EdgeData>::receive_graph_handler,
                           NULL, 0, mgraph, owner2vertex, vertex2owner);
          logger(LOG_INFO, "sending Graph to %d", i);
        }
      }
      else {
        // otherwise, I shall wait until the graph is received
        while(!distributed) {         
           sched_yield();
        }
      }
      logger(LOG_INFO, "%d: receiving graph BARRIER", dc.procid());
      dc.barrier();
    }

    static cloned_graph<VertexData, EdgeData>* receive_target;
    // message handler which receives the graph
    static void receive_graph_handler(distributed_control& dc, size_t source, 
                      void* ptr, size_t len, 
                      graph<VertexData, EdgeData> &g,                              
                      std::vector<std::vector<vertex_id_t> > owner2vertex,
                      std::vector<procid_t> vertex2owner) {
      logger(LOG_INFO, "%d: receiving graph", dc.procid());
      receive_target->mgraph = g;
      receive_target->owner2vertex = owner2vertex;
      receive_target->vertex2owner= vertex2owner;
      receive_target->distributed = true;
      logger(LOG_INFO, "%d: receiving graph done", dc.procid());
    }

    // message handler which receives the graph
    static void update_vertex_handler(distributed_control& dc, size_t source, 
                      void* ptr, size_t len, 
                      vertex_id_t v, VertexData &vdata) {
      receive_target->vertex_data(v) = vdata;
    }

    static void update_edge_handler(distributed_control& dc, size_t source, 
                      void* ptr, size_t len, 
                      edge_id_t e, EdgeData &edata) {
      receive_target->edge_data(e) = edata;
    }

    void update_vertex(vertex_id_t v) const {
      DCHECK_LT(v, num_vertices());
      for (procid_t i = 0;i < dcontrol->numprocs(); ++i) {
        if (i != myprocid) {
          dcontrol->remote_callxs(i, cloned_graph<VertexData, EdgeData>::update_vertex_handler,
                          NULL, 0, v, vertex_data(v));
        }
      }
    }

    void update_edge(edge_id_t e) const {
      DCHECK_LT(e, num_edges());
      for (procid_t i = 0;i < dcontrol->numprocs(); ++i) {
        if (i != myprocid) {
          dcontrol->remote_callxs(i, cloned_graph<VertexData, EdgeData>::update_edge_handler,
                            NULL, 0, e, edge_data(e));
        }
      }
    }

    // Structural Mutators================================================
    // These functions are not callable once graph is distributed
    /**
     * Clear all internal data
     */
    void clear() {
      mgraph.clear();
    }
    

    void finalize() {
      ASSERT_FALSE(distributed);
      mgraph.finalize();
    }
   
    /** 
     * Creates a vertex containing the vertex data and returns the id
     * of the new vertex id.
     */
    vertex_id_t add_vertex(const VertexData& vdata = VertexData() ) {
      ASSERT_FALSE(distributed);
      return mgraph.add_vertex(vdata);
    }

    
    /**
     * Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_t add_edge(vertex_id_t source, vertex_id_t target, 
                          const EdgeData& edata = EdgeData()) {
      ASSERT_FALSE(distributed);
      return mgraph.add_edge(source,target,edata);
    } 
        
    cloned_graph<VertexData,EdgeData>&
      operator=(const graph<VertexData, EdgeData> &g) { 
      ASSERT_FALSE(distributed);
      mgraph = g;
      return *this;
    }
     
     
    // Structural Query Functions =============================================


    bool has_constant_edges() {
    	return constant_edges;
    }
    
    void set_constant_edges(bool b) {
    	constant_edges = b;
    }
    

    /** Get the number of vetices */
    size_t num_vertices() const {
      return mgraph.num_vertices();
    }

    /** Get the number of edges */
    size_t num_edges() const {
      return mgraph.num_edges();
    } 


    /** Get the number of in edges */
    size_t num_in_neighbors(vertex_id_t v) const {
      return mgraph.num_in_neighbors(v);
    } 
    
    /** get the number of out edges */
    size_t num_out_neighbors(vertex_id_t v) const  {
      return mgraph.num_out_neighbors(v);
    } 

    /** Find an edge */
    std::pair<bool, edge_id_t>
    find(vertex_id_t source, vertex_id_t target) const {
      return mgraph.find(source,target);
    } 

    
    /** get the edge id for the edge */
    edge_id_t edge_id(vertex_id_t source, vertex_id_t target) const {
      return mgraph.edge_id(source,target);
    } 

    
    /** get the reverser edge id for the edge */
    edge_id_t rev_edge_id(edge_id_t eid) const {
      return mgraph.rev_edge_id(eid);
    }


    /** get the source of the edge */
    vertex_id_t source(edge_id_t edge_id) const {
      return mgraph.source(edge_id);
    }

    /** get the dest of the edge */
    vertex_id_t target(edge_id_t edge_id) const {
      return mgraph.target(edge_id);
    }
    
        
    /** Get the ids of the in edges */
    edge_list in_edge_ids(vertex_id_t v) const {
      return mgraph.in_edge_ids(v);
    } 

    /** Get the ids of the out edges */
    edge_list out_edge_ids(vertex_id_t v) const {
      return mgraph.out_edge_ids(v);
    } 

    void metis_partition(size_t nparts, std::vector<uint32_t>& ret_part) {
      mgraph.metis_partition(nparts, ret_part);
    } 

    void bfs_partition(size_t nparts, std::vector<uint32_t> &vertex2part) {
      mgraph.bfs_partition(nparts, vertex2part);
    } 

    void random_partition(size_t nparts, std::vector<uint32_t> &vertex2part) {
      mgraph.random_partition(nparts, vertex2part);
    }

    void partition(partition_method::partition_method_enum partmethod,
                  size_t nparts, std::vector<uint32_t> &vertex2part) {
      mgraph.partition(partmethod, nparts, vertex2part);
    }

    bool topological_sort(std::vector<vertex_id_t> &topsort) const {
      return mgraph.topological_sort(topsort);
    }

    const std::vector<vertex_id_t>& my_vertices() const{
      return owner2vertex[myprocid];
    }

    const procid_t owner(const vertex_id_t &v) const{
      return vertex2owner[v];
    }

    // Data Functions =============================================
    // These functions will only work on local data
    
    /** Get the vertex data */
    VertexData& vertex_data(vertex_id_t v) {
      return mgraph.vertex_data(v);
    } 
    
    /** Get the vertex data */
    const VertexData& vertex_data(vertex_id_t v) const {
      return mgraph.vertex_data(v);
    } 

    /** Get the edge_data */
    EdgeData& edge_data(vertex_id_t source, vertex_id_t target) {
      return mgraph.edge_data(source, target);
    } 
    
    /** Get the edge_data */
    const EdgeData& edge_data(vertex_id_t source, vertex_id_t target) const {
      return mgraph.edge_data(source, target);
    } 

    /** Get the edge_data */
    EdgeData& edge_data(edge_id_t edge_id) { 
      return mgraph.edge_data(edge_id);
    }
    
    /** Get the edge_data */
    const EdgeData& edge_data(edge_id_t edge_id) const {
      return mgraph.edge_data(edge_id);
    }



    // Data Functions =============================================
    // serialization functions. These only serialize the graph. Load may
    // not be called once distributed


    /** Load the graph from an archive */
    void load(iarchive &arc) {
      ASSERT_FALSE(distributed);
      mgraph.load(arc);
    } 

    /** Save the graph to an archive */
    void save(oarchive &arc) const {
      mgraph.save(arc);
    }
    

    /** Load the graph from a file */
    void load(const std::string& filename) {
      mgraph.load(filename);
    }

    
    void save(const std::string& filename) const {
      mgraph.save(filename);
    }

    
    void save_adjacency(const std::string& filename) const {
      mgraph.save_adjacency(filename);
    }

    
    
  private:    
    graph<VertexData, EdgeData> mgraph;
    std::vector<std::vector<vertex_id_t> > owner2vertex;
    std::vector<procid_t> vertex2owner;
    
    // remember the pointer to the distributed controller as well as my procid
    procid_t myprocid;
    distributed_control *dcontrol;
    
    // whether I am distributed;
    volatile bool distributed;
    
    bool constant_edges;

  }; // End of graph

  template<typename VertexData, typename EdgeData> 
  cloned_graph<VertexData, EdgeData>* cloned_graph<VertexData, EdgeData>::receive_target = NULL;
} // end of namespace graphlab


#endif
