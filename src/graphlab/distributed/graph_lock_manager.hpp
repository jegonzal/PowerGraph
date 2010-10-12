#ifndef DISTRIBUTED_GRAPH_LOCK_MANAGER
#define DISTRIBUTED_GRAPH_LOCK_MANAGER
#include <vector>
#include <queue>
#include <boost/unordered_map.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope_factory.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/lock_manager.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/blocking_queue.hpp>
#include <graphlab/util/synchronized_unordered_map.hpp>
#include <graphlab/util/synchronized_unordered_map2.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {
#define MAX_TOTAL_REFERENCES 5

class dist_scope_request{
 public:
  vertex_id_t vertex;
  scope_range::scope_range_enum scoperange;
  dist_scope_request():vertex(-1), scoperange(scope_range::USE_DEFAULT) { }
  
  dist_scope_request(vertex_id_t vertex,
                     scope_range::scope_range_enum scoperange):
            vertex(vertex),scoperange(scoperange){ }
 private:
  size_t requestid;       // the lock request id which satisfied this scope request
  template <typename Graph>
  friend class graph_lock_manager;
};


/**
Manages the locks on a distributed graph. Allows thread/distributed safe
acquiring and releasing of graph data.

\TODO: Cache size control
*/
template<typename Graph>
class graph_lock_manager{
 public:
  /******************    Internal Structs ******************/

    struct request_descriptor{
      size_t requestid;
      vertex_id_t vertex;
      scope_range::scope_range_enum scoperange;

      size_t blockid;
      bool ready;
      
      procid_t acquiringfromproc;
      std::vector<bool> needvertexdata;
      std::vector<std::vector<lockrequest> > proc2lockrequests;
      std::vector<std::vector<edge_id_t> > edgerequests;
      
      // the lock only locks these structures below
      mutex lock;
      size_t referencecounter;
      std::queue<size_t> associatedids;
      size_t childof;
      size_t totalreferences; // the total number of references ever acquired
    };

    struct vertex_data_wrapper {
      typename Graph::vertex_data_type vdata;
      atomic<int> numlocks;
      bool modified;
    };

    struct edge_data_wrapper {
      typename Graph::edge_data_type edata;
      atomic<int> numlocks;
      bool modified;
    };

    struct block_request {
      size_t blockid;
      std::queue<dist_scope_request> newsatisfied;
      std::queue<dist_scope_request> pendingrequests;
      size_t numreqs; // total number of requests
      size_t numreleased; // number of requests released
      size_t numsatisfied;  // number if requests satisfied
      size_t numissued;  // number of requests where the lock request has been issued
      size_t throttle_rate;
      std::set<size_t> reqids;
/*      std::vector<size_t> taken;
      std::vector<size_t> returned;*/
      mutex lock;
    };
    
    
    class lock_send_request_thread:public thread {
      graph_lock_manager<Graph> &owner;
     public: 
      lock_send_request_thread(graph_lock_manager<Graph> &glm):thread(NULL, 0), 
                                                               owner(glm) { }
      void run() {
        logger(LOG_INFO, "background locking thread started");
        while(1) {
          std::pair<size_t, bool> ret = owner.progress_requests.dequeue();
          if (ret.second == false) break;
          owner.progress_lock(ret.first);
        }
      }
    };
    class background_unlock_thread:public thread {
      graph_lock_manager<Graph> &owner;
     public: 
      background_unlock_thread(graph_lock_manager<Graph> &glm):thread(NULL, 0), 
                                                               owner(glm) { }
      void run() {
        logger(LOG_INFO, "background unlocking thread started");
        while(1) {
          std::pair<size_t, bool> ret = owner.unlock_requests.dequeue();
          if (ret.second == false) break;
          owner.release_blocking(ret.first);
        }
      }
    };
 private:
    distributed_control &dc;
    distributed_lock_manager<Graph> &dlm;
    Graph &graph;
    
    atomic<size_t> lastreqid; /// Keep an incremental counter over request ids
    
    typedef synchronized_unordered_map<request_descriptor> reqdesc_container;
    reqdesc_container activerequests;

    typedef synchronized_unordered_map<vertex_data_wrapper> remote_vdata_type;
    remote_vdata_type remote_vdata;
    
    typedef synchronized_unordered_map<edge_data_wrapper> remote_edata_type;
    remote_edata_type remote_edata;

    // block requests
    atomic<size_t> blockreqid; /// Keep an incremental counter over block ids
    boost::unordered_map<size_t, block_request> blockreqs;
    rwlock blockreqslock;
    
    blocking_queue<size_t> progress_requests;
    blocking_queue<size_t> unlock_requests;
    lock_send_request_thread *locking_thread;
    background_unlock_thread *unlocking_thread;
    std::vector<mutex> referencecounter_locks;
    synchronized_unordered_map2<std::set<size_t> > vertex2reqids;  // only hold parent requests
    
    bool caching;
    bool const_edges;
    bool use_adjacent_vertices;
    bool vertex_scope_pushed_updates;
    dense_bitset dirtyvertices;
    dense_bitset dirtyedges;
    atomic<size_t> vertexpacks;
    atomic<size_t> vertexunpacks;
    atomic<size_t> edgepacks;
    atomic<size_t> edgeunpacks;
    atomic<size_t> reusedscopes;
 public:
  /** Constructor. Creates the set of locks associated with this machine.
  Graph should already be partitioned and distributed at this point   */
  graph_lock_manager(distributed_control &dc,
                     distributed_lock_manager<Graph> &dlm,
                     Graph &graph) :
           dc(dc), dlm(dlm), graph(graph), lastreqid(0), activerequests(11), remote_vdata(11), remote_edata(17),
           vertex2reqids(131071){
    referencecounter_locks.resize(131071);  //2^17 - 1
    graph_lock_manager_target = this;
    locking_thread = new lock_send_request_thread(*this);
    locking_thread->start();
    dirtyvertices.resize(graph.num_vertices());
    dirtyvertices.clear();
    dirtyedges.resize(graph.num_local_edges());
    dirtyedges.clear();
    unlocking_thread = new background_unlock_thread(*this);
    unlocking_thread->start();
    caching = false;
    dlm.set_caching(caching);
    const_edges = false;
    use_adjacent_vertices = true;
    vertex_scope_pushed_updates = false;
    vertexpacks.value = 0;
    edgepacks.value = 0;
    vertexunpacks.value = 0;
    edgeunpacks.value = 0;
    reusedscopes.value = 0;
  }
  
  void set_caching(bool _caching) {
    caching = _caching;
    dlm.set_caching(caching);
  }

  void set_vertex_scope_pushed_updates(bool _vertex_scope_pushed_updates) {
    vertex_scope_pushed_updates = _vertex_scope_pushed_updates;
    dlm.set_vertex_scope_pushed_updates(vertex_scope_pushed_updates, pushed_cache_update);
  }

  void set_constant_edges(bool _const_edges) {
    const_edges = _const_edges;
  }

  void set_use_adjacent_vertices(bool _use_adjacent_vertices) {
    use_adjacent_vertices = _use_adjacent_vertices;
  }

  void init() {
    dc.barrier();
    // fill caches
    if (caching) {
      logger(LOG_INFO, "Setting up caches...");
      for (size_t i = 0;i < graph.my_vertices().size(); ++i) {
        fixed_dense_bitset<32> vertexsent2procs;
        vertexsent2procs.clear();
        // send this vertex to all machines who might need it
        foreach(edge_id_t eid, graph.in_edge_ids(graph.my_vertices()[i])) {
          procid_t remotemachine = graph.owner(graph.source(eid));
          if (remotemachine != dc.procid() && 
              !vertexsent2procs.get(remotemachine)) {
            dc.remote_callxs(remotemachine, 
                           pushed_cache_update,
                           NULL,
                           0,
                           graph.my_vertices()[i],
                           graph.vertex_data(graph.my_vertices()[i]));
          }
          if (remotemachine != dc.procid()) {
            // always send the edge
            dc.remote_callxs(remotemachine, 
                  pushed_cache_update_edge,
                  NULL,
                  0,
                  eid,
                  graph.edge_data(eid));
          }
        }
        
        // do the same for outedges
        foreach(edge_id_t eid, graph.out_edge_ids(graph.my_vertices()[i])) {
          procid_t remotemachine = graph.owner(graph.target(eid));
          if (remotemachine != dc.procid() && 
              !vertexsent2procs.get(remotemachine)) {
            dc.remote_callxs(remotemachine, 
                           pushed_cache_update,
                           NULL,
                           0,
                           graph.my_vertices()[i],
                           graph.vertex_data(graph.my_vertices()[i]));
          }
        }
      }
    }
    logger(LOG_INFO, "Waiting for barrier.");
    dc.barrier();
    logger(LOG_INFO,"Cache set up complete.");
    dlm.clear_dirty_flags();
  }
  void print_stats() {
    std::cout << "Scopes issued: " << lastreqid.value << std::endl;
    std::cout << "Scopes reused: " << reusedscopes.value << std::endl;
  }
  ~graph_lock_manager() {
    progress_requests.stop_blocking();
    locking_thread->join();
    unlock_requests.stop_blocking();
    unlocking_thread->join();
    delete locking_thread;
    delete unlocking_thread;
  }
  /**
    Requests a collection of locks
  */
  size_t block_deferred_lock(const std::vector<dist_scope_request>& reqs, size_t throttle_rate = 500) {
    // find a new block id
    size_t newblockid = blockreqid.inc();
    blockreqslock.writelock();
    block_request &b = blockreqs[newblockid];
    blockreqslock.unlock();

    // create the block request
    b.numsatisfied = 0;
    b.numreqs = 0;
    b.numissued = 0;
    b.numreleased = 0;
    b.blockid = newblockid;
    b.throttle_rate = throttle_rate;
    block_add_deferred_lock(newblockid, reqs);
    return newblockid;
  }
  /** Adds requests to a block */
  void block_add_deferred_lock(size_t blockid, const std::vector<dist_scope_request>& reqs) {
    if (reqs.size() == 0) return;
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "adding requests to non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    blockreqslock.unlock();
    
    blockreq.lock.lock();
    blockreq.numreqs += reqs.size();
    blockreq.lock.unlock();
    // issue all the block requests
    for (size_t i = 0; i < reqs.size(); ++i) {
      // if we are holding more than BLOCK_MAX_HELD_LOCKS locks in this block
      if (blockreq.numissued - blockreq.numreleased >= blockreq.throttle_rate) {
        // hold it in the pending requests queue
        blockreq.lock.lock();
        blockreq.pendingrequests.push(reqs[i]);
        blockreq.lock.unlock();
      }
      else {
        size_t newreqid = lastreqid.inc();
        blockreq.lock.lock();
        blockreq.numissued++;
        blockreq.reqids.insert(newreqid);
        blockreq.lock.unlock();
        deferred_lock_private(reqs[i].vertex,
                              reqs[i].scoperange,
                              blockid, newreqid);
      }
    }
  }
  
  /** Adds requests to a block */
  void block_add_deferred_lock(size_t blockid, const dist_scope_request& req) {
//     std::cout << "Request: " << req.vertex << std::endl;
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "adding requests to non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    blockreqslock.unlock();
    
    blockreq.lock.lock();
    blockreq.numreqs++;
    // if we are holding more than BLOCK_MAX_HELD_LOCKS locks in this block
    if (blockreq.numissued - blockreq.numreleased >= blockreq.throttle_rate) {
      // hold it in the pending requests queue
      blockreq.pendingrequests.push(req);
      blockreq.lock.unlock();
    }
    else {
      // issue all the block requests
      size_t newreqid = lastreqid.inc();
      blockreq.reqids.insert(newreqid);
      blockreq.numissued++;
      blockreq.lock.unlock();
      deferred_lock_private(req.vertex,
                            req.scoperange,
                            blockid, newreqid);
    }
  }
  /**
    Returns the new set of locks acquired through newlocks.
    Return value contains the number of locks remaining that
    has not been acquired yet.
  */
  size_t block_status(size_t blockid, std::vector<dist_scope_request> &newlocks) {
    newlocks.clear();
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "Requesting status of non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    blockreqslock.unlock();
    newlocks.reserve(blockreq.newsatisfied.size());
    blockreq.lock.lock();
    while(!blockreq.newsatisfied.empty()){
      newlocks.push_back(blockreq.newsatisfied.front());
      blockreq.newsatisfied.pop();
    }
    size_t retval = blockreq.numreqs - blockreq.numsatisfied;
    blockreq.lock.unlock();
    return retval;
  }

  /**
    Returns a newly acquired lock in newlock.
    If no locks are available, newlock.vertex == -1
    Return value contains the number of locks remaining that
    has not been acquired yet.
  */
  size_t block_status(size_t blockid, dist_scope_request &newlock) {
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "Requesting status of non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    blockreqslock.unlock();
    blockreq.lock.lock();
    if (blockreq.newsatisfied.empty() == false) {
      newlock = blockreq.newsatisfied.front();
//       blockreq.taken.push_back(newlock.vertex);
      blockreq.newsatisfied.pop();
    }
    else {
      newlock.vertex = vertex_id_t(-1);
    }
    size_t retval = blockreq.numreqs - blockreq.numsatisfied;
    blockreq.lock.unlock();
    return retval;
  }
  
  /** Release part of a block that is complete. This must
     contain dist_scope_requests returned by block_status */

  void block_release_partial(size_t blockid,
                             const std::vector<dist_scope_request> &releases) {
    // search for the block
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "Requesting release partial of non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    blockreqslock.unlock();
    // remove this set of requests from the set of requests the block
    // is managing
    blockreq.lock.lock();
    // 
    blockreq.numreleased += releases.size();
    for (size_t i = 0;i < releases.size(); ++i) {
      std::set<size_t>::iterator reqiter = blockreq.reqids.find(releases[i].requestid);
      if (reqiter == blockreq.reqids.end()) {
        logstream(LOG_FATAL) << "Requesting partial block release of request id "
                            << releases[i].requestid
                            <<  " which does not belong to block "
                            << blockid << std::endl;
        ASSERT_TRUE(false);
      }
      blockreq.reqids.erase(reqiter);
    }
    blockreq.lock.unlock();

    // issue the release
    for (size_t i = 0;i < releases.size(); ++i) {
      release(releases[i].requestid);
    }
    block_try_acquire_from_pending_queue(blockreq);
  }


  void block_release_partial(size_t blockid,
                             const dist_scope_request &releasereq) {
//     std::cout << "Release: " << releasereq.vertex << std::endl;
    // search for the block
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "Requesting release partial of non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    
    blockreqslock.unlock();
    // remove this set of requests from the set of requests the block
    // is managing
    blockreq.lock.lock();
//     blockreq.returned.push_back(releasereq.vertex);
    blockreq.numreleased++;
    std::set<size_t>::iterator reqiter = blockreq.reqids.find(releasereq.requestid);
    if (reqiter == blockreq.reqids.end()) {
      logstream(LOG_FATAL) << "Requesting partial block release of request id "
                          << releasereq.requestid
                          <<  " which does not belong to block "
                          << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    blockreq.reqids.erase(reqiter);

    blockreq.lock.unlock();

    // issue the release
    release(releasereq.requestid);
    block_try_acquire_from_pending_queue(blockreq);
  }
  
  /** Releases a block */
  void block_release(size_t blockid) {
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "Requesting release of non-existant block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }
    block_request& blockreq = iter->second;
    blockreqslock.unlock();
    // lock the block
    blockreq.lock.lock();
    if (blockreq.numsatisfied != blockreq.numreqs) {
      logstream(LOG_FATAL) << "Cannot release incomplete block "
                           << blockid << std::endl;
      ASSERT_TRUE(false);
    }

    if (!blockreq.pendingrequests.empty()) {
      logstream(LOG_FATAL) << "Block " << blockid << " still has pending requests!" << std::endl;
      ASSERT_TRUE(false);
    }
    // release everything
    blockreq.numreleased += blockreq.reqids.size();
    foreach(size_t reqid, blockreq.reqids){
      release(reqid);
    }
    blockreq.lock.unlock();
    
    blockreqslock.writelock();
    blockreqs.erase(blockid);
    blockreqslock.unlock();
  }
  
  size_t deferred_lock(const vertex_id_t &vertex,
                               scope_range::scope_range_enum scoperange) {
    return deferred_lock_private(vertex, scoperange, -1, lastreqid.inc());
  }

  /** Returns true once all locks in the request has been acquired. 
      False otherwise */
  bool locks_ready(const size_t &reqid) const {
    typename reqdesc_container::const_datapointer i = activerequests.find(reqid);
    assert(i.first);
    return i.second->ready;
  }

  /** Gets a constant reference to the vertex data associated with a vertex.
      This vertex must already be acquired through a deferred lock.
      The lock must be ready before the data is available.
      This function will quit with a fatal exception if the conditions
      are not satisfied.*/
  const typename Graph::vertex_data_type& get_const_vertex(const vertex_id_t &v) const {
    if (graph.owner(v) == dc.procid()) return graph.vertex_data(v);
    else {
      typename remote_vdata_type::const_datapointer it = remote_vdata.find(v);
      DASSERT_TRUE(it.first);
      const vertex_data_wrapper &vdatawrapper = *(it.second);
      return vdatawrapper.vdata;
    }
  }
  
  const typename Graph::vertex_data_type& get_vertex(const vertex_id_t &v) const {
    return get_const_vertex(v);
  }
  /** Gets a constant reference to the edge data associated with an edge.
      The destination vertex of the edge must already be acquired through 
      a deferred lock. The lock must be ready before the data is available.
      This function might quit with a fatal exception if the conditions
      are not satisfied.*/
  const typename Graph::edge_data_type& get_const_edge(const edge_id_t &e) const {
    if (graph.owner(graph.target(e)) == dc.procid()) {
      return graph.edge_data(e); 
    }
    else {
      typename remote_edata_type::const_datapointer it = remote_edata.find(e);
      DASSERT_TRUE(it.first);
      const edge_data_wrapper &edatawrapper = *it.second;
      return edatawrapper.edata;
    }
  }


  const typename Graph::edge_data_type& get_edge(const edge_id_t &e) const {
    return get_const_edge(e);
  }


  /** Gets a modifiable reference to the vertex data associated with a vertex.*/
  typename Graph::vertex_data_type& get_vertex(const vertex_id_t &v) {
    if (graph.owner(v) == dc.procid()) {
      //logstream(LOG_INFO) << "local set dirty on " << v << std::endl;
      dirtyvertices.set_bit(v);
      return graph.vertex_data(v);
    }
    else {
      typename remote_vdata_type::datapointer it = remote_vdata.find(v);
      DASSERT_TRUE(it.first);
      vertex_data_wrapper &vdatawrapper = *(it.second);
      vdatawrapper.modified = true;
      return vdatawrapper.vdata;
    }
  }
  
  /** Gets a modifiable reference to the edge data associated with an edge.*/
  typename Graph::edge_data_type& get_edge(const edge_id_t &e) {
    if (graph.owner(graph.target(e)) == dc.procid()) {
      dirtyedges.set_bit(graph.global_to_local_eid(e));
      return graph.edge_data(e); 
    }
    else {      
      typename remote_edata_type::datapointer it = remote_edata.find(e);
      DASSERT_TRUE(it.first);
      edge_data_wrapper &edatawrapper = *it.second;
      edatawrapper.modified = true;
      return edatawrapper.edata;
    }
  }

  /** Releases the request, and commits all modifications. 
  All references acquired from get_vertex(...) or get_edge(...) associated
  with this request will be invalidated. This release is non-blocking */
  void release(const size_t &reqid) {
    return release_blocking(reqid);
  }
  
  /** Releases the request, and commits all modifications. 
  All references acquired from get_vertex(...) or get_edge(...) associated
  with this request will be invalidated */
  void release_blocking(const size_t &reqid) {
//    logstream(LOG_INFO) << "release of " << reqid << std::endl;

    typename reqdesc_container::datapointer it = activerequests.find(reqid);
    DASSERT_TRUE(it.first);
    request_descriptor &desc = *(it.second);
    
    // unlocking a child request
    // decrement the parent, and call for its release if its counter goes to 0
    if (desc.childof != size_t(-1)){
      typename reqdesc_container::datapointer parent = activerequests.find(desc.childof);
      DASSERT_TRUE(parent.first);
      request_descriptor &parentdesc = *(parent.second);
      // erase this child request
      activerequests.erase(reqid);
      // try to release my parent.
      release_blocking(parentdesc.requestid);
      return;
    }
    
    // try to remove myself from the vertex2reqids vector
    vertex2reqids.write_critical_section(desc.vertex);
    std::pair<bool, std::set<size_t>* > vt = vertex2reqids.find(desc.vertex);
    assert(vt.first);
    std::set<size_t>& curvt = *(vt.second);
    desc.lock.lock();
    desc.referencecounter--;
    if (desc.referencecounter == 0) {
        // we are clear for deletion!
      curvt.erase(reqid);
      desc.lock.unlock();
      if (curvt.size() == 0) {
        vertex2reqids.erase(desc.vertex);
      }
      vertex2reqids.release_critical_section(desc.vertex);
    }
    else {
      // we can't delete yet. someone still has a reference to me
      ASSERT_FALSE(desc.associatedids.empty());
      // wake up a child
      size_t child = desc.associatedids.front();
      desc.associatedids.pop();
      desc.lock.unlock();
      vertex2reqids.release_critical_section(desc.vertex);
      progress_lock(child); //progress lock on a child will always immediately succeed
      return;
    }
    
    // issue the unlocks in reverse order though I don't think this is too 
    // important
    //
    for (int i = desc.proc2lockrequests.size() - 1; i >= 0 ; --i) {
      if (desc.proc2lockrequests[i].size() > 0 || desc.edgerequests[i].size() > 0) {
        std::map<vertex_id_t, typename Graph::vertex_data_type> vdatamap;
        std::map<edge_id_t, typename Graph::edge_data_type> edatamap;
        if (i != dc.procid()) {
          // I only need to commit the writes on remote machines
          // grab the vertex and edge data we need to commit
          pack_vertex_writes(desc.proc2lockrequests[i], vdatamap);
          pack_edge_writes(desc.edgerequests[i], edatamap);
        }
        else {
          // this is myself
          for (size_t j = 0;j < desc.proc2lockrequests[i].size(); ++j) {
            if (dirtyvertices.clear_bit(desc.proc2lockrequests[i][j].vertex)) {
              dlm.set_dirty(desc.proc2lockrequests[i][j].vertex, dc.procid());
            }
          }
          for (size_t j = 0;j < desc.edgerequests[i].size(); ++j) {
            if (dirtyedges.clear_bit(desc.edgerequests[i][j])) {
              dlm.set_dirty_edge(desc.edgerequests[i][j], dc.procid());
            }
          }
        }
        // issue the unlocks
        dc.remote_callxs(i,
            distributed_lock_manager<Graph>::request_unlock_write_handler,
            NULL, 0, 
            desc.proc2lockrequests[i],    // vertices to un lock
            vdatamap,                     // vdata to commit
            edatamap);                    // edata to commit
      }
    }      
    activerequests.erase(reqid);
    
  }
  
  
  /*********        Message Handlers *****************/
  static graph_lock_manager* graph_lock_manager_target;
  
  static void pushed_cache_update(distributed_control& dc, size_t source,
                void* unused, size_t len, vertex_id_t vid,
                typename Graph::vertex_data_type &vdata) {
    graph_lock_manager_target->insert_vdata(vid, vdata, false);
  }
  
  static void pushed_cache_update_edge(distributed_control& dc, size_t source,
                void* unused, size_t len, edge_id_t eid,
                typename Graph::edge_data_type &edata) {
    graph_lock_manager_target->insert_edata(eid, edata, false);
  }
  
  static void lock_response(distributed_control& dc, size_t source,
                  void* unused, size_t len, handlerarg_t requestid,
                std::map<vertex_id_t, typename Graph::vertex_data_type> &vdata,
                std::map<edge_id_t, typename Graph::edge_data_type> &edata) {
    if (source != dc.procid()) {
      typename reqdesc_container::datapointer it = graph_lock_manager_target->activerequests.find(requestid);
      DASSERT_TRUE(it.first);
      request_descriptor &desc = *(it.second);

      // write the vdata/edata into the remote_vdata/remote_edata structures
      //ASSERT_EQ(desc.proc2lockrequests[source].size(), vdata.size());
      //ASSERT_EQ(desc.edgerequests[source].size(), edata.size());
      
      for (size_t i = 0; i < desc.proc2lockrequests[source].size(); ++i) {
        vertex_id_t vid = desc.proc2lockrequests[source][i].vertex;
        typename std::map<vertex_id_t, typename Graph::vertex_data_type>::iterator iter = 
                                                                          vdata.find(vid);
        if (iter != vdata.end()) {
          ASSERT_TRUE(desc.needvertexdata[source]);
          graph_lock_manager_target->insert_vdata(iter->first, iter->second);
          graph_lock_manager_target->vertexunpacks.inc();
        }
        // we need to increment the numlocks as well even if we did not get the data
        else {
          //std::cout << ".";
          ASSERT_TRUE(graph_lock_manager_target->caching);
         // ASSERT_FALSE(desc.needvertexdata[source]);
         size_t reflockid = vid % graph_lock_manager_target->referencecounter_locks.size();
          graph_lock_manager_target->referencecounter_locks[reflockid].lock();
          typename remote_vdata_type::datapointer iter = graph_lock_manager_target->remote_vdata.find(vid);
          
          ASSERT_TRUE(iter.first);
          
          iter.second->numlocks.inc();
          graph_lock_manager_target->referencecounter_locks[reflockid].unlock();

        }
      
      }

      
      // write the edgedata
      for (size_t i = 0; i < desc.edgerequests[source].size(); ++i) {
        edge_id_t eid = desc.edgerequests[source][i];;
        typename std::map<edge_id_t, typename Graph::edge_data_type>::iterator iter = 
                                                                  edata.find(eid);
        if (iter != edata.end()) {
          graph_lock_manager_target->insert_edata(iter->first, iter->second);
          graph_lock_manager_target->edgeunpacks.inc();
        }
        else {
          ASSERT_TRUE(graph_lock_manager_target->caching);
         // ASSERT_FALSE(desc.needvertexdata[source]);
          size_t reflockid = eid % graph_lock_manager_target->referencecounter_locks.size();
          graph_lock_manager_target->referencecounter_locks[reflockid].lock();
          typename remote_edata_type::datapointer iter = graph_lock_manager_target->remote_edata.find(eid);
          ASSERT_TRUE(iter.first);
          iter.second->numlocks.inc();
          graph_lock_manager_target->referencecounter_locks[reflockid].unlock();
        }
      }
    }
    // push the lock
    graph_lock_manager_target->progress_lock(requestid);
    //graph_lock_manager_target->progress_requests.enqueue(requestid);
  }
  


 private:


  inline bool can_avoid_vtx_request(vertex_id_t vid, lock_type locktype) {
    return locktype == NOLOCK && 
          (graph.owner(vid) == dc.procid() || 
          (caching && graph.owner(vid) != dc.procid() && remote_vdata.find(vid).first)); 
          
  }
  /** Requests a lock on a set of vertices and a set of edges.*/
  size_t deferred_lock_private(const vertex_id_t &vertex,
                               scope_range::scope_range_enum scoperange,
                               size_t blockid,
                               size_t newreqid) {
    // build the request datastructre
    request_descriptor &desc = activerequests.insert(newreqid, request_descriptor());
    
    // get an new requestid
    desc.blockid = blockid;
    desc.requestid = newreqid;
    desc.ready = false;
    desc.vertex = vertex;
    desc.scoperange = scoperange;
    desc.acquiringfromproc = 0;
    desc.proc2lockrequests.resize(dc.numprocs());
    desc.edgerequests.resize(dc.numprocs());
    desc.needvertexdata.resize(dc.numprocs());
    desc.childof = size_t(-1);
    desc.referencecounter = 1;
    desc.totalreferences = 1;
    
    
    // search for a request which will suit our needs
    // bool newready = false;
    vertex2reqids.write_critical_section(desc.vertex);
    std::pair<bool, std::set<size_t>*> vt = vertex2reqids.insert_with_failure_detect(desc.vertex, std::set<size_t>());
    std::set<size_t>& curvt = *(vt.second);
    foreach(size_t other, curvt) {
      typename reqdesc_container::datapointer it = activerequests.find(other);
      DASSERT_TRUE(it.first);
      request_descriptor &otherdesc = *(it.second);
      if (scope_is_subset_of(scoperange, otherdesc.scoperange)) {
         // we are in lock!
         // increment reference of parent
         otherdesc.lock.lock();
         if (otherdesc.totalreferences >= MAX_TOTAL_REFERENCES) {
           // this scope has been overused
           // disallow
           otherdesc.lock.unlock();
           continue;
         }
         otherdesc.referencecounter++;
         otherdesc.totalreferences++;
         otherdesc.associatedids.push(newreqid);
         // set myself as child
         desc.childof = other;
         // I am always done. (Rather I am done when my parent is done)
         desc.acquiringfromproc = dc.numprocs();
         otherdesc.lock.unlock();
         break;
      }
    }
    // if I got hooked up
    if (desc.childof != size_t(-1)) {
      vertex2reqids.release_critical_section(desc.vertex);
      reusedscopes.inc();
      return newreqid;
    }
    else {
      curvt.insert(desc.requestid);
      vertex2reqids.release_critical_section(desc.vertex);
    }
    // fill the lock requests
    // loop over all the in and out neighbors
    edge_list inedges =  graph.in_edge_ids(vertex);
    edge_list outedges = graph.out_edge_ids(vertex);
    // type of lock to use for neighbors
    lock_type nbrlocktype = RDLOCK;
    lock_type vtxlocktype = WRLOCK;
    switch(scoperange){
      case scope_range::VERTEX_CONSISTENCY:
        vtxlocktype = WRLOCK;
        nbrlocktype = NOLOCK;
        break;
      case scope_range::VERTEX_READ_CONSISTENCY:
        vtxlocktype = RDLOCK;
        nbrlocktype = NOLOCK;
        break;
      case scope_range::READ_CONSISTENCY:
        nbrlocktype = RDLOCK;
        vtxlocktype = RDLOCK;
        break;
      case scope_range::EDGE_CONSISTENCY:
        vtxlocktype = WRLOCK;
        nbrlocktype = RDLOCK;
        break;
      case scope_range::FULL_CONSISTENCY:
        vtxlocktype = WRLOCK;
        nbrlocktype = WRLOCK;
        break;
      case scope_range::NULL_CONSISTENCY:
        nbrlocktype = NOLOCK;
        vtxlocktype = NOLOCK;
        break;
      case scope_range::USE_DEFAULT:
        ASSERT_MSG(false, "No default consistency model set");
    }

    // this is annoying. Is there a shorter way to do this
    // build up the lock requests by looping through all the inedges/outedges
    bool curlocked = false;
    size_t inidx = 0;
    size_t outidx = 0;
    size_t numv = graph.num_vertices();
    vertex_id_t curv = vertex;
    vertex_id_t inv  = (inedges.size() > 0) ? graph.source(inedges[0]) : numv;
    vertex_id_t outv  = (outedges.size() > 0) ? graph.target(outedges[0]) : numv;
    // iterate both in order and lock
    // include the current vertex in the iteration
    // check
 /*   for (size_t i = 1;i < inedges.size(); ++i) {
      ASSERT_LE(graph.source(inedges[i-1]), graph.source(inedges[i]));
    }
    for (size_t i = 1;i < outedges.size(); ++i) {
      ASSERT_LE(graph.target(outedges[i-1]), graph.target(outedges[i]));
    }*/
    while (inidx < inedges.size() || outidx < outedges.size()) {
      if (!curlocked && curv < inv  && curv < outv) {
        procid_t owner = graph.owner(curv);
        if (!can_avoid_vtx_request(curv, vtxlocktype)) {
          desc.proc2lockrequests[owner].push_back(lockrequest(curv, vtxlocktype));
        }
        curlocked = true;
        curv = numv;
      } else if (inv < outv) {
        procid_t owner = graph.owner(inv);
        if (!can_avoid_vtx_request(inv, nbrlocktype)) {
          desc.proc2lockrequests[owner].push_back(lockrequest(inv, nbrlocktype));
        }
        ++inidx;
        inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
      } else if (outv < inv) {
        procid_t owner = graph.owner(outv);
        if (!can_avoid_vtx_request(outv, nbrlocktype)) {
          desc.proc2lockrequests[owner].push_back(lockrequest(outv, nbrlocktype));
        }
        ++outidx;
        outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
      } else if (inv == outv){
        procid_t owner = graph.owner(inv);
        if (!can_avoid_vtx_request(inv, nbrlocktype)) {
          desc.proc2lockrequests[owner].push_back(lockrequest(inv, nbrlocktype));
        }
        ++inidx; ++outidx;
        inv = (inedges.size() > inidx) ? graph.source(inedges[inidx]) : numv;
        outv= (outedges.size() > outidx) ? graph.target(outedges[outidx]) : numv;
      }
    }
    if (!curlocked) {
      procid_t owner = graph.owner(curv);
      if (!can_avoid_vtx_request(curv, vtxlocktype)) {
        desc.proc2lockrequests[owner].push_back(lockrequest(curv, vtxlocktype));
      }
    }
       
    // begin request
    // add the edges
    // if I am the owner of the root vertex, I do not need to requeset these added

    procid_t rootowner = graph.owner(vertex);
    if (rootowner != dc.procid()) {
      if (const_edges) {
        foreach(const edge_id_t& eid, inedges) {
          if (remote_edata.find(eid).first == false) {
            desc.edgerequests[rootowner].push_back(eid);
          }
        }  
      }
      else {
        inedges.fill_vector(desc.edgerequests[rootowner]);
      }
    }


    
    // add all the outedges
    for (size_t i = 0;i < outedges.size(); ++i) {
      vertex_id_t targetvertex = graph.target(outedges[i]);
      procid_t targetvertexowner = graph.owner(targetvertex);
      if (targetvertexowner != dc.procid()) {
        if (const_edges || can_avoid_vtx_request(targetvertex,nbrlocktype)) {
          if (remote_edata.find(outedges[i]).first == false) {
            desc.edgerequests[targetvertexowner].push_back(outedges[i]);
          }
        }
        else {
          desc.edgerequests[targetvertexowner].push_back(outedges[i]);
        }
      }
    }
    
    for (procid_t i = 0; i < dc.numprocs(); ++i) {
      if (use_adjacent_vertices) {
        // we only need it if it is not local
        desc.needvertexdata[i] = i != dc.procid();
      }
      else {
        // ok. I only need vertex data if 
        // this is a remote call, && I am requesting for the root vertex && and the root vertex is remote
        desc.needvertexdata[i] = (i != dc.procid()) && (rootowner == i);
      }
    }
    


    // find the first cpu we need to acquire from
    while(desc.acquiringfromproc < dc.numprocs() && 
          desc.proc2lockrequests[desc.acquiringfromproc].size() == 0 && 
          desc.edgerequests[desc.acquiringfromproc].size() == 0) { 
      ++desc.acquiringfromproc;
    }
    // push the lock
    //progress_lock(desc.requestid);
    // if I am self-acquiring. get it immediately
    // otherwise defer to the lock progress thread
    if (desc.acquiringfromproc == dc.procid()) {
      progress_lock(desc.requestid);
    }
    else {
      progress_requests.enqueue(desc.requestid);
    }
    return newreqid;
  }

  /// progress the lock acquisition
  void progress_lock(size_t requestid) {
    
    typename reqdesc_container::datapointer it = activerequests.find(requestid);
    DASSERT_TRUE(it.first);
    request_descriptor &desc = *(it.second);

    // find the next proc we need to acquire from
    while(desc.acquiringfromproc < dc.numprocs() && 
          desc.proc2lockrequests[desc.acquiringfromproc].size() == 0 &&
          desc.edgerequests[desc.acquiringfromproc].size() == 0) { 
      ++desc.acquiringfromproc;
    }

    // if I am out of processors, I am done.
    if (desc.acquiringfromproc == dc.numprocs()) {
      // Ready!
      desc.ready = true;
      update_block_request(desc);    
    }
    else {
      // issue the lock request
      int curproc = desc.acquiringfromproc;
      desc.acquiringfromproc++;

      dc.remote_callxs(curproc,
          distributed_lock_manager<Graph>::request_lock_handler,
          NULL, 0, 
          desc.proc2lockrequests[curproc],// vertices to lock
          requestid,                                     // local requestid
          (size_t)(graph_lock_manager<Graph>::lock_response),                
          desc.edgerequests[curproc],     // edges to read
          desc.needvertexdata[curproc]);     // get the edgedata if not local
    }
  }
  
  
  /** Pack the vertex writes into a map structure */
  void pack_vertex_writes(const std::vector<lockrequest> &lockreqs, 
              std::map<vertex_id_t, typename Graph::vertex_data_type> &vdatamap) {
    // for each vertex
    for (size_t j = 0;j < lockreqs.size(); ++j) {
      vertex_id_t vertex_to_unlock = lockreqs[j].vertex;
//      logger(LOG_INFO, "%d: Decrement of vertex %d", dc.procid(), vertex_to_unlock);
      typename remote_vdata_type::datapointer iter = remote_vdata.find(vertex_to_unlock);
      DASSERT_TRUE(iter.first);
      ASSERT_GT(iter.second->numlocks.value, 0);
      
      bool wasmodified = false;
      vertex_data_wrapper &vdatawrapper = *(iter.second);
      atomic_exchange(vdatawrapper.modified, wasmodified);
      
      if (wasmodified) {
        vdatamap[vertex_to_unlock] = vdatawrapper.vdata;
        vertexpacks.inc();
      }
      referencecounter_locks[vertex_to_unlock % referencecounter_locks.size()].lock();
      // decrement the lock count
      int remaininglocks = vdatawrapper.numlocks.dec();
      referencecounter_locks[vertex_to_unlock % referencecounter_locks.size()].unlock();
      //if I am the last one, try to delete it
      if (remaininglocks == 0) try_delete_vdata(vertex_to_unlock);
    }
  }

  /** Pack the edge writes into a map structure */
void pack_edge_writes(const std::vector<edge_id_t> &edgelist, 
              std::map<edge_id_t, typename Graph::edge_data_type> &edatamap) {
    // for each vertex
    for (size_t j = 0;j < edgelist.size(); ++j) {
      edge_id_t edge = edgelist[j];
      // look for it in remote_vdata
      // acquire a read lock to read the vertex data
      typename remote_edata_type::datapointer it = remote_edata.find(edge);
      DASSERT_TRUE(it.first);
      edge_data_wrapper& edw = *(it.second);
      ASSERT_GT(edw.numlocks.value, 0);
      bool wasmodified = false;
      atomic_exchange(edw.modified, wasmodified);
      if (wasmodified) {
        edatamap[edge] = edw.edata;
        edgepacks.inc();
      }
      // but happily references are always good.

      // decrement the lock count
      referencecounter_locks[edge % referencecounter_locks.size()].lock();
      int remaininglocks = edw.numlocks.dec();
      referencecounter_locks[edge % referencecounter_locks.size()].unlock();
      //if I am the last one, try to delete it
      if (remaininglocks == 0) try_delete_edata(edge);
    }
  }
  /** Inserts vdata into the remote_vdata structure. Incrementing
      the counter if it already exists */
  void insert_vdata(vertex_id_t vertex, 
                    const typename Graph::vertex_data_type &vdata,
                    bool increment = true) {
    referencecounter_locks[vertex % referencecounter_locks.size()].lock();

    typename remote_vdata_type::datapointer iter = remote_vdata.find(vertex);
    // if I am already in remote_vdata. just increment the locks counter
    if (iter.first) {
//      logger(LOG_INFO, "%d: Increment of vertex %d", dc.procid(), vertex);
      if (increment) iter.second->numlocks.inc();
      iter.second->vdata = vdata;
    }
    else {
      // try to insert into the vdata
      vertex_data_wrapper newvdatawrapper;
      newvdatawrapper.vdata = vdata;
      if (increment) newvdatawrapper.numlocks.value = 1;
      else newvdatawrapper.numlocks.value = 0;
      newvdatawrapper.modified = false;
      typename remote_vdata_type::datapointer newinsertion = remote_vdata.insert_with_failure_detect(vertex, 
                                                                                                  newvdatawrapper);
      if (newinsertion.first == false) {
        if (increment) newinsertion.second->numlocks.inc();
        newinsertion.second->vdata = vdata;
      }
    }
    referencecounter_locks[vertex % referencecounter_locks.size()].unlock();
  }
  /** Inserts vdata into the remote_vdata structure. Incrementing
      the counter if it already exists */
  void insert_edata(edge_id_t edge, 
                    const typename Graph::edge_data_type &edata,
                    bool increment = true) {

    referencecounter_locks[edge % referencecounter_locks.size()].lock();
    typename remote_edata_type::datapointer it = remote_edata.find(edge);
    // if I am already in remote_vdata. just increment the locks counter
    if (it.first) {
      if (increment) it.second->numlocks.inc();
      it.second->edata = edata;
    }
    else {
      edge_data_wrapper newedatawrapper;
      newedatawrapper.edata = edata;
      if (increment) newedatawrapper.numlocks.value = 1;
      else newedatawrapper.numlocks.value = 0;
      newedatawrapper.modified = false;
      typename remote_edata_type::datapointer newinsertion = remote_edata.insert_with_failure_detect(edge, 
                                                                                                  newedatawrapper);
      if (newinsertion.first == false) {
        newinsertion.second->numlocks.inc();
        newinsertion.second->edata = edata;
      }
    }
    referencecounter_locks[edge % referencecounter_locks.size()].unlock();
  }
  
  template <typename WrapperType>
  static bool FUNCTOR_if_numlocks_zero(const WrapperType &vdata) {
    return vdata.numlocks.value == 0;
  }
  /** Deletes vdata if the counter goes to 0 */
  void try_delete_vdata(vertex_id_t vertex) {
    if (caching) return;
    referencecounter_locks[vertex % referencecounter_locks.size()].lock();
    remote_vdata.erase_if(vertex, FUNCTOR_if_numlocks_zero<vertex_data_wrapper>);
    referencecounter_locks[vertex % referencecounter_locks.size()].unlock();
  }
  
  /** Deletes edata if the counter goes to 0*/
  void try_delete_edata(edge_id_t eid) {
    if (caching || const_edges) return;
    referencecounter_locks[eid % referencecounter_locks.size()].lock();
    remote_edata.erase_if(eid, FUNCTOR_if_numlocks_zero<edge_data_wrapper>);
    referencecounter_locks[eid % referencecounter_locks.size()].unlock();
  }
  
  void block_try_acquire_from_pending_queue(block_request & blockreq) {
// I released some stuff. see if I have any left in the pending queue
    // I can grab
    blockreq.lock.lock();
    while(blockreq.numissued - blockreq.numreleased < blockreq.throttle_rate
          && !blockreq.pendingrequests.empty()) {
      dist_scope_request req = blockreq.pendingrequests.front();
      blockreq.pendingrequests.pop();
    
      size_t newreqid = lastreqid.inc();
      blockreq.numissued++;
      blockreq.reqids.insert(newreqid);
      blockreq.lock.unlock();
      deferred_lock_private(req.vertex,
                              req.scoperange,
                              blockreq.blockid, newreqid);
      blockreq.lock.lock();
    }
    blockreq.lock.unlock();
  }

  void update_block_request(request_descriptor &desc) {
    if (desc.blockid == size_t(-1)) return;
    blockreqslock.readlock();
    typename boost::unordered_map<size_t, block_request>::iterator iter =
                                          blockreqs.find(desc.blockid);
    if (iter == blockreqs.end()) {
      logstream(LOG_FATAL) << "Successful lock meant for lock block "
                            << desc.blockid << ". But lock block is missing. " << std::endl;
      ASSERT_TRUE(false);
    }
    block_request &blockreq = iter->second;
    blockreqslock.unlock();
    dist_scope_request distreq(desc.vertex, desc.scoperange);
    distreq.requestid = desc.requestid;
    blockreq.lock.lock();
    blockreq.numsatisfied++;
    blockreq.newsatisfied.push(distreq);
    blockreq.lock.unlock();
  }
};

template<typename Graph>
graph_lock_manager<Graph>* graph_lock_manager<Graph>::graph_lock_manager_target;
}

#include <graphlab/macros_undef.hpp>
#endif
