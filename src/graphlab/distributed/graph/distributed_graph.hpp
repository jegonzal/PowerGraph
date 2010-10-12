#ifndef DISTRIBUTED_PARTITIONED_GRAPH_HPP
#define DISTRIBUTED_PARTITIONED_GRAPH_HPP
#include <cassert>
#include <cmath>
#include <vector>

#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/util/generics/blob.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed/distributed_control.hpp>

#include <graphlab/macros_def.hpp>



namespace graphlab { 
  

  template<typename VertexData, typename EdgeData> class distributed_graph;
  typedef distributed_graph<blob, blob> blob_distributed_graph;

  
  /**
   * \class distributed_graph
   *   
   * Distributed graph. Each partition is created by each processor
   * separately.  Vertices are cloned, but edges are NOT.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_graph {
  public:

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;
      
    
  public:

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    // TODO: remove code duplication
    distributed_graph() : 
      mgraph(0),global_id_counter(0)  { 
      receive_target = this;
      constant_edges = false;
      only_local_edges = false;
      dcontrol = NULL;
      checksum = size_t(-1);
      logger(LOG_INFO, "%d: initialized distributed graph (wout dc)", myprocid);
    }
   
   
    distributed_graph(distributed_control &dc) : 
      mgraph(0), global_id_counter(0)  { 
      dcontrol = &dc;
      myprocid = dc.procid();
      receive_target = this; 
      constant_edges = false;
      only_local_edges = (dc.numprocs() > 1);
      checksum = 0;
      for(int i=0; i<dc.numprocs(); i++) {
        receiverlist.push_back(i);
      }
      std::random_shuffle(receiverlist.begin(), receiverlist.end());
      logger(LOG_INFO, "%d: initialized distributed graph", myprocid);
    }

    
    

    static distributed_graph<VertexData, EdgeData>* receive_target;
    
    // checksum check
    static void check_checksum(distributed_control& dc, size_t source, void * ptr, size_t len, size_t check) {
      printf("%d Checksum check: %ld / %ld\n", (int) source, check, receive_target->checksum);

      ASSERT_MSG(check == receive_target->checksum, "Proc id %d has different idea of vertex ownerships than proc id %d, checksums  %ld vs. %ld", 
                 (int) source, 0, check, receive_target->checksum);
      printf("%d Checksum check passed: %ld\n", (int) source, check);
    }
    
    
    // message handler which receives the graph
    static void dist_update_vertex_handler(distributed_control& dc, size_t source, 
                                           void* ptr, size_t len, 
                                           vertex_id_t v, VertexData vdata) {
      receive_target->vertex_data(v) = vdata;
      if (v == 141400) printf("%d received 141400\n", (int) dc.procid());
    }

    static void dist_update_edge_handler(distributed_control& dc, size_t source, 
                                         void* ptr, size_t len, 
                                         edge_id_t e, EdgeData edata) {
      receive_target->edge_data(e) = edata;
    }

    void update_vertex(vertex_id_t v) const {
      DCHECK_LT(v, num_vertices());
      
      std::vector<bool> send_to_proc(dcontrol->numprocs(), false);
      // Always send to 0
      // send_to_proc[0] = true;
      
      // Loop thru out edges to see if need update
      foreach(edge_id_t e, out_edge_ids(v)) {
        send_to_proc[vertex2owner[target(e)]] = true;
      }
      foreach(edge_id_t e, in_edge_ids(v)) {
        send_to_proc[vertex2owner[source(e)]] = true;
      }
     
      for(int j=0; j<dcontrol->numprocs(); ++j) {
        size_t i = receiverlist[j];
        if (i != myprocid && send_to_proc[i]) {
          dcontrol->remote_callxs(i, distributed_graph<VertexData, EdgeData>::dist_update_vertex_handler,
                                  NULL, 0, v, vertex_data(v));
        }
      }
    }
    
    void broadcast_vertex_data(vertex_id_t v) const {
      DCHECK_LT(v, num_vertices());    
      // Vertex data is broadcasted to everyone. TODO: send only to owner?
      for(int j=0; j<dcontrol->numprocs(); ++j) {
        size_t i = receiverlist[j];
        if (i != myprocid) {
          dcontrol->remote_callxs(i, distributed_graph<VertexData, EdgeData>::dist_update_vertex_handler,
                                  NULL, 0, v, vertex_data(v));
        }
      }
    }
    void send_vertices_to_proczero() {
      if (myprocid != 0) {
        foreach(vertex_id_t v, myvertices) {
          dcontrol->remote_callxs(0, distributed_graph<VertexData, EdgeData>::dist_update_vertex_handler,
                                  NULL, 0, v, vertex_data(v));
        }
      }
    }

    void update_edge(edge_id_t e) const {
    
      vertex_id_t _source = source(e);
      vertex_id_t _target = target(e);
      

      // Send updated edge data only to source and target.
      if (vertex2owner[_target] != myprocid) {
        dcontrol->remote_callxs(vertex2owner[_target], distributed_graph<VertexData, EdgeData>::dist_update_edge_handler,
                                NULL, 0, e, edge_data(e));
      }
      if (vertex2owner[_source] != myprocid) {
        dcontrol->remote_callxs(vertex2owner[_source], distributed_graph<VertexData, EdgeData>::dist_update_edge_handler,
                                NULL, 0, e, edge_data(e));
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
    
    struct edge_id_less_functor {
      distributed_graph<VertexData, EdgeData> * g_ptr;
      edge_id_less_functor(distributed_graph<VertexData, EdgeData> * g_ptr) : g_ptr(g_ptr) { }
      bool operator()(edge_id_t a, edge_id_t b) {
        if(g_ptr->source(a) == g_ptr->source(b)) return g_ptr->target(a) < g_ptr->target(b);
        else return g_ptr->source(a) < g_ptr->source(b);
      }
    };


    
    void finalize_dist(bool donotsort=false) {
      ASSERT_EQ(receive_target, this);
      ASSERT_EQ(dcontrol->procid(), myprocid);
      
      // We want every node to send updates in different order
      receiverlist.clear();
      for(int i=0; i<dcontrol->numprocs(); i++) {
        receiverlist.push_back((i+dcontrol->procid())%dcontrol->numprocs());
      }
      std::random_shuffle(receiverlist.begin(), receiverlist.end());

      
      // Compute checksum of vertex ownerships
      size_t cs = 0;
      for(size_t i=0; i<mgraph.num_vertices(); i++) {
        cs += (i % 31) * vertex2owner[i];
      }
      checksum = cs;
      dcontrol->mpi_barrier();
      printf("%d My checksum :%ld %ld\n", dcontrol->procid(), checksum, receive_target->checksum);
      if (dcontrol->procid() != 0)
        dcontrol->remote_callxs(0, distributed_graph::check_checksum, NULL, 0, checksum);
   
      long int overhead = 0;
  	 
      if (g_inedges.size() > 0) return;
      // Hack
		  
      if (!only_local_edges) {
        g_inedges.resize(mgraph.num_vertices());
        g_outedges.resize(mgraph.num_vertices());
	
        for(size_t i=0; i<mgraph.num_vertices(); i++) {
          if (vertex2owner[i] == myprocid) {
            edge_list ine = mgraph.in_edge_ids(i);
            for(size_t j=0; j<ine.size(); j++) {
              g_inedges[i].push_back(eid_local_to_global[ine[j]]);
              overhead += sizeof(edge_id_t);
            }
				 
            edge_list oute = mgraph.out_edge_ids(i);
            for(size_t j=0; j<oute.size(); j++) {
              g_outedges[i].push_back(eid_local_to_global[oute[j]]);
              overhead += sizeof(edge_id_t);
            }
          }
        }
    
      	edge_id_less_functor less_functor(this);
      	if (donotsort == false) {
          	for(size_t i=0; i<mgraph.num_vertices(); i++) {
              std::sort(g_inedges[i].begin(), g_inedges[i].end(), less_functor);
              std::sort(g_outedges[i].begin(), g_outedges[i].end(), less_functor);
      	    }
      	}
      }
      
      logger(LOG_INFO, "Finalized distributed graph.. Overhead: %ld bytes", overhead);
      logger(LOG_INFO, "Size of global-local mappings: %d edges to care", 
             (int) eid_local_to_global.size());
    }
    
    
    void finalize() {
      mgraph.finalize();
      finalize_dist();
    }
   
    /** 
     * Creates a vertex containing the vertex data and returns the id
     * of the new vertex id.
     * Note: each processor needs to create the same vertices in same order!!
     */
    vertex_id_t add_vertex(int procid, const VertexData& vdata = VertexData()) {
      if (dcontrol != NULL) ASSERT_TRUE(procid >= 0 && procid < dcontrol->numprocs());
    
      // Create new vertex
      vertex_id_t vid = mgraph.add_vertex(vdata);
      // Store vertex owner.
      vertex2owner.push_back(procid);
      if (procid == myprocid) {
        myvertices.push_back(vid);
      }
      return vid;
    }

    
    /**
     * Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared. Note, edge is created locally only
     * if either target or source is this processor id.
     */
    edge_id_t add_edge(vertex_id_t source, vertex_id_t target, 
                       const EdgeData& edata = EdgeData()) {
      ASSERT_TRUE(vertex2owner.size() > source && vertex2owner.size() > target);
      if (vertex2owner[source] == myprocid || vertex2owner[target] == myprocid) {
        edge_id_t local_id = mgraph.add_edge(source,target,edata);
        if (!only_local_edges) {
          eid_local_to_global.push_back(global_id_counter);
          eid_global_to_local[global_id_counter] = local_id;
        } else {
          return local_id;
        }	
      }
      return global_id_counter++;
    } 
    
    
    edge_id_t local_to_global_eid(const edge_id_t &e) const{
      return eid_local_to_global[e];
    }
    
    edge_id_t global_to_local_eid(const edge_id_t &e) const{
      edge_id_t loceid = (eid_global_to_local.find(e))->second;
      return loceid;
    }
    
  
    distributed_graph<VertexData,EdgeData>&
    operator=(const graph<VertexData, EdgeData> &g) { 
      mgraph = g;
      return *this;
    }
     
     
    // Structural Query Functions =============================================


    /** Get the number of vetices */
    size_t num_vertices() const {
      return mgraph.num_vertices();
    }

    /** Get the number of edges */
    size_t num_edges() const {
      DASSERT_MSG(false, "Querying edge count on distributed graph not allowed.");
    } 
    
    size_t num_local_edges() const {
      return mgraph.num_edges();
    } 

    /** Get the number of in edges */
    size_t num_in_neighbors(vertex_id_t v) const {
      ASSERT_EQ(vertex2owner[v], myprocid);
      return mgraph.num_in_neighbors(v);
    } 
    
    /** get the number of out edges */
    size_t num_out_neighbors(vertex_id_t v) const  {
      ASSERT_EQ(vertex2owner[v], myprocid);
      return mgraph.num_out_neighbors(v);
    } 

    /** Find an edge */
    std::pair<bool, edge_id_t>
    find(vertex_id_t source, vertex_id_t target) const {
      std::pair<bool, edge_id_t> loc = mgraph.find(source,target);
      if (loc.first == true && !only_local_edges) 
      	loc.second = eid_local_to_global[loc.second];
      return loc;
    } 

    
    /** get the edge id for the edge */
    edge_id_t edge_id(vertex_id_t source, vertex_id_t target) const {
      return (only_local_edges ? mgraph.edge_id(source,target) :
              eid_local_to_global[mgraph.edge_id(source,target)]);
    } 

    
    /** get the reverser edge id for the edge */
    edge_id_t rev_edge_id(edge_id_t eid) const {
      edge_id_t loceid = (only_local_edges ? eid : (eid_global_to_local.find(edge_id))->second);
      return eid_local_to_global[mgraph.rev_edge_id(loceid)];
    }


    /** get the source of the edge */
    vertex_id_t source(edge_id_t edge_id) const {
      edge_id_t loceid = (only_local_edges ? edge_id : (eid_global_to_local.find(edge_id))->second);
      return mgraph.source(loceid);
    }

    /** get the dest of the edge */
    vertex_id_t target(edge_id_t edge_id) const {
      edge_id_t loceid = (only_local_edges ? edge_id : (eid_global_to_local.find(edge_id))->second);
      return mgraph.target(loceid);
    }


    /** Get the ids of the in edges */
    edge_list in_edge_ids(vertex_id_t v) const {
      ASSERT_EQ(vertex2owner[v], myprocid);
      return (only_local_edges ? mgraph.in_edge_ids(v) : g_inedges[v]);
    } 

    /** Get the ids of the out edges */
    edge_list out_edge_ids(vertex_id_t v) const {
      ASSERT_EQ(vertex2owner[v], myprocid);
      return (only_local_edges ? mgraph.out_edge_ids(v) : g_outedges[v]);
    } 


    const std::vector<vertex_id_t>& my_vertices() const{
      return myvertices;
    }

    const procid_t owner(const vertex_id_t &v) const{
      return vertex2owner[v];
    }
    
    procid_t myproc() {
      return myprocid;
    }

    // Data Functions =============================================
    // These functions will only work on local data
    
    bool has_constant_edges() {
      return constant_edges;
    }
    
    void set_constant_edges(bool b) {
      constant_edges = b;
    }
    
    void set_local_edges(bool b) {
      only_local_edges = b;
    }
    
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
      ASSERT_TRUE(vertex2owner[source] == myprocid || vertex2owner[target] == myprocid);
      return mgraph.edge_data(source, target);
    } 
    
    /** Get the edge_data */
    const EdgeData& edge_data(vertex_id_t source, vertex_id_t target) const {
      ASSERT_TRUE(vertex2owner[source] == myprocid || vertex2owner[target] == myprocid);
      return mgraph.edge_data(source, target);
    } 

    /** Get the edge_data */
    EdgeData& edge_data(edge_id_t edge_id) { 
      edge_id_t loceid = (only_local_edges ? edge_id : eid_global_to_local[edge_id]);
      ASSERT_TRUE(vertex2owner[mgraph.source(loceid)] == myprocid || vertex2owner[mgraph.target(loceid)] == myprocid);
      return mgraph.edge_data(loceid);
    }
    
    /** Get the edge_data */
    const EdgeData& edge_data(edge_id_t edge_id) const {
      edge_id_t loceid  =  (only_local_edges ? edge_id : (eid_global_to_local.find(edge_id))->second);
      ASSERT_TRUE(vertex2owner[mgraph.source(loceid)] == myprocid || vertex2owner[mgraph.target(loceid)] == myprocid);
      return mgraph.edge_data(loceid);
    }


    // Data Functions =============================================
    // serialization functions. These only serialize the graph. Load may
    // not be called once distributed

    /** Load the graph from an archive */
    void load(iarchive &arc) {
      // read the vertices and colors
      arc >> mgraph
          >> vertex2owner
          >> eid_local_to_global
          >> myvertices
          >> myprocid;
          
      // Create global_to_local mapping
      
      for(unsigned int i=0; i<eid_local_to_global.size(); i++) {
        eid_global_to_local[eid_local_to_global[i]] = i;
      }
      
      // Optimization
      if (eid_local_to_global.size() == 0 && mgraph.num_edges() > 0) {
        only_local_edges = true;
      }
      
      logstream(LOG_INFO) << "Number of edges: " << mgraph.num_edges() << ", loctoglobal: " <<
        eid_global_to_local.size() << " / " << eid_local_to_global.size() << std::endl;
      finalize_dist();
    } // end of load

    /** Save the graph to an archive */
    void save(oarchive &arc) const {
      // Write the number of edges and vertices
      arc << mgraph
          << vertex2owner
          << eid_local_to_global
          << myvertices
          << myprocid;
    } // end of save
    
  	
    // Dirty, lazy..
    static char * filename_for_part(const std::string& basename, int procid, int numofprocs) {
      static char fname[255];
      sprintf(fname, "%s_%dof%d.gpart", basename.c_str(), procid+1, numofprocs);
      printf("%s\n", fname);
      return fname;
    }	
    

    /** Load the graph from a file. Loads only the part for this procid. */
    void load(const std::string& filename, distributed_control &_dc) {
      dcontrol = &_dc;
      std::ifstream fin(filename_for_part(filename, _dc.procid(), _dc.numprocs()));
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();

    } // end of load


    /**
     * \brief save the graph to the file given by the filename
     * 
     */    
    void save(const std::string& filename, int procid, int numofprocs) const {
      std::ofstream fout(filename_for_part(filename, procid, numofprocs));
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save
    
    
    static void partition_graph_tofile(graph<VertexData, EdgeData> &_mgraph, int nparts, partition_method::partition_method_enum partmethod,
                                       std::string &basename) {
      std::vector<uint32_t> vertex2part;
      logger(LOG_INFO, "Going to partition graph to %d parts", nparts);
      _mgraph.partition(partmethod, nparts, vertex2part);
      logger(LOG_INFO, "Partition ready");

      // Create partition files
      for(int partnum=0; partnum<nparts; partnum++) {
        logger(LOG_INFO, "Storing partnum %d", partnum);

        distributed_graph<VertexData, EdgeData> dgraph;
        dgraph.myprocid = partnum;
        for (vertex_id_t vid=0; vid<_mgraph.num_vertices(); vid++) {
          vertex_id_t nvid = dgraph.add_vertex(vertex2part[vid], _mgraph.vertex_data(vid));
          ASSERT_EQ(vid, nvid);
        }
        for(edge_id_t e=0; e<_mgraph.num_edges(); e++) {
          edge_id_t ne = dgraph.add_edge(_mgraph.source(e), _mgraph.target(e), _mgraph.edge_data(e));
          ASSERT_EQ(e, ne);
        }
    	 	 
        // Save
        dgraph.save(basename, partnum, nparts);
      }
    }

    graph<VertexData, EdgeData> mgraph;

  private:    
    std::vector<procid_t> vertex2owner;
    std::vector<vertex_id_t> myvertices;
    
    bool constant_edges;
    bool only_local_edges;
        
    std::vector< std::vector< edge_id_t > >  g_inedges;    
    std::vector< std::vector< edge_id_t > >  g_outedges;    
        
    size_t checksum;
    
    // remember the pointer to the distributed controller as well as my procid
    procid_t myprocid;
    distributed_control *dcontrol;
    
    // local to global edge id mapping
    boost::unordered_map<edge_id_t, edge_id_t> eid_global_to_local;
    std::vector<edge_id_t> eid_local_to_global;
    
    edge_id_t global_id_counter;

    std::vector<size_t> receiverlist;

  }; // End of graph

  template<typename VertexData, typename EdgeData> 
  distributed_graph<VertexData, EdgeData>* distributed_graph<VertexData, EdgeData>::receive_target = NULL;

} // end of namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif
