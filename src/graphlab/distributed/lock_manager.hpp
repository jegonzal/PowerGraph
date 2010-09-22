#ifndef DISTRIBUTED_LOCK_MANAGER_HPP
#define DISTRIBUTED_LOCK_MANAGER_HPP
#include <vector>
#include <queue>
#include <map>
#include <logger/assertions.hpp>
#include <graphlab/graph/graph.hpp>
#include <serialization/serialization_includes.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/synchronized_unordered_map.hpp>
#include <graphlab/macros_def.hpp>

namespace graphlab {


enum lock_type{
  RDLOCK, WRLOCK, NOLOCK
};

/**
This struct describes a single request for a vertex lock
*/
struct lockrequest{
  vertex_id_t vertex;
  lock_type locktype;
  lockrequest() {}
  lockrequest(vertex_id_t vertex, lock_type locktype):
          vertex(vertex), locktype(locktype) {}

  void save(oarchive &arc) const{
    arc << vertex;
    arc << (size_t)locktype;
  }
  
  void load(iarchive &arc){
    arc >> vertex;
    size_t val;
    arc >> val;
    locktype = (lock_type)val;
  }

  bool operator<(const lockrequest &req) const {
    return vertex < req.vertex;
  }
  bool operator==(const lockrequest &req) const {
    return vertex == req.vertex;
  }
};

/**
This is a basic implementation of a distributed lock manager on graph data.
This distributed lock manager protects only the vertices assigned to this
machine (as read from the DistributedGraph passed in via the constructor).

To acquire locks, the caller must signal the message handler (request_lock_handler).
The caller must take the precautions to signal the right processor. If the caller
signals a processor who does not own the vertex/lock requested, it will
immediately raise a fatal exception.
*/
template <typename Graph>
class distributed_lock_manager {

 public:
  /************************ Internal Structures **************************/
  /**
  This is an "unsynchronized rwlock".  It requires support from the
  lock manager to manage its synchronization.
  */
  struct distrwlock{
    atomic<uint16_t> numreaders;
    bool writer;
    inline distrwlock():numreaders(0), writer(false){ }
    inline void readlock() {
      ASSERT_FALSE(writer);
      numreaders.inc();
    }
    inline void writelock() {
      ASSERT_EQ(numreaders.value, 0);
      ASSERT_FALSE(writer);
      writer = true;
    }
    inline void unlock() {
      // if I am the writer, I should be the only one here.
      if (writer) {
        ASSERT_EQ(numreaders.value, 0);
        writer = false;
      }
      else {
        ASSERT_GT(numreaders.value, 0);
        ASSERT_FALSE(writer);
        numreaders.dec();
      }
    }
  };

  /// this is the required type of the response handler
  typedef void (*lock_response_handler)(distributed_control& dc, size_t source,
                  void* unused, size_t len, size_t requestid,
                  std::map<vertex_id_t, typename Graph::vertex_data_type> &vdata,
                  std::map<edge_id_t, typename Graph::edge_data_type> &edata);

  typedef void (*pushed_cache_update_handler_type)(distributed_control& dc, size_t source,
                                    void* unused, size_t len, vertex_id_t vid,
                                    typename Graph::vertex_data_type &vdata);
  /// internal structure describing a lock request
  struct lock_descriptor{
    size_t requestid; /// an id of the request, set by the requester
    size_t sourcemachine; /// source machine which originated this request
    std::vector<lockrequest> lockrequests;  /// sequence of locks to fulfill
    size_t lockindex; /// inndex of the last unfulfilled lock
    mutex lock;  ///  the lock on this lockset
    lock_response_handler responsehandler;  /// the handler I reply to
    bool withdata;
    std::vector<edge_id_t> edgerequests;

  };

 private:
  /***********************   Private Variables  **************************/
  
  distributed_control &dc;  // reference to the communication later
  Graph &graph;             // reference to the graph I am protecting
  
  typedef fixed_dense_bitset<32> procid_set;
  std::map<vertex_id_t, size_t> vertex2locknumber;
  // the list of locks we are managing
  std::vector<distrwlock> locks;
  std::vector<dense_bitset> proc2dirty;
  std::vector<dense_bitset> proc2dirtyedge;
  // a queue over pending requests for each lock
  std::vector<std::list<lock_descriptor*> > writelockqueue;
  std::vector<std::list<lock_descriptor*> > readlockqueue; 
  std::vector<mutex> lockqueue_locks;
  synchronized_unordered_map<procid_set> vertexnbrs;
  bool vertex_scope_pushed_updates;
  pushed_cache_update_handler_type push_receiver_function;
  bool caching;
 public:

  /***********************   Constructor **************************/
  distributed_lock_manager(distributed_control &dc_, Graph &graph_):
                                                  dc(dc_),graph(graph_),
                                                  vertexnbrs(8191),
                                                  vertex_scope_pushed_updates(false),
                                                  caching(false){
    size_t numlocks = graph.my_vertices().size();
    locks.resize(numlocks);
    readlockqueue.resize(numlocks);
    writelockqueue.resize(numlocks);
    lockqueue_locks.resize(numlocks);
    // set the target of signals
    lock_manager_target = this;

    size_t locknumber = 0;
    for (size_t i = 0;i < graph.my_vertices().size(); ++i) {
      vertex2locknumber[graph.my_vertices()[i]] = locknumber;
      ++locknumber;
    }
    proc2dirty.resize(dc.numprocs());
    proc2dirtyedge.resize(dc.numprocs());
    for (size_t i = 0; i < dc.numprocs(); ++i) {
      proc2dirty[i].resize(numlocks);
      proc2dirty[i].fill();
      proc2dirtyedge[i].resize(graph.num_local_edges());
      proc2dirtyedge[i].fill();
    }
  }

  void clear_dirty_flags() {
    for (size_t i = 0; i < dc.numprocs(); ++i) {
      proc2dirty[i].clear();
      proc2dirtyedge[i].clear();
    }
  }
  void print_queues() {
    for (size_t i = 0;i < readlockqueue.size(); ++i) {
      if (readlockqueue[i].size() > 0 || writelockqueue[i].size() > 0 ||
          locks[i].numreaders.value > 0 || locks[i].writer) {
          std::cout <<dc.procid() <<">"  << i << ": " << locks[i].writer << " " << locks[i].numreaders.value
                    << ": " << readlockqueue[i].size() << " " << writelockqueue[i].size() << std::endl;
          // get the lock index of the heads of each queue
          if (!readlockqueue[i].empty()) std::cout << "\tr: " << readlockqueue[i].front()->lockindex << std::endl;
          if (!writelockqueue[i].empty()) std::cout << "\tw: " << writelockqueue[i].front()->lockindex << std::endl;
      }
    }
  }
  void set_vertex_scope_pushed_updates(bool _vertex_scope_pushed_updates,
                                       pushed_cache_update_handler_type _receiver_function) {
    vertex_scope_pushed_updates = _vertex_scope_pushed_updates;
    push_receiver_function = _receiver_function;
  }
  
  void set_caching(bool _caching) {
    caching = _caching;
  }
  
  void set_dirty_edge(edge_id_t edge, size_t exceptsource = -1) {
    if (!caching) return;
    size_t localedgeid = graph.global_to_local_eid(edge);
    for (size_t i = 0;i < dc.numprocs(); ++i) {
      if (i != exceptsource) proc2dirtyedge[i].set_bit(localedgeid);
    }
  }
    
  void set_dirty(vertex_id_t vertex, size_t exceptsource = -1) {
    if (!vertex_scope_pushed_updates) {
      if (!caching) return;
      size_t locknumber = vertex2locknumber[vertex];
      // if not pushing. just set cache to dirty
      for (size_t i = 0;i < dc.numprocs(); ++i) {
        if (i != exceptsource) proc2dirty[i].set_bit(locknumber);
      }

    }
    else {
      // if pushed. no dirty needed
      if (vertexnbrs.find(vertex).first == false) {
        // create the entry
        procid_set pset;
        foreach (edge_id_t inedge, graph.in_edge_ids(vertex)) {
          procid_t sourceowner = graph.owner(graph.source(inedge));
          if (sourceowner != dc.procid()) {
            pset.set_bit(sourceowner);
          }
        }
        
        foreach (edge_id_t outedge, graph.out_edge_ids(vertex)) {
          procid_t targetowner = graph.owner(graph.target(outedge));
          if (targetowner != dc.procid()) {
            pset.set_bit(targetowner);
          }
        }
        vertexnbrs.insert(vertex, pset);
      }
      // see who I need to send to
      std::pair<bool, procid_set*> iter = vertexnbrs.find(vertex);
      ASSERT_TRUE(iter.first);
      uint32_t b = 0;
      if (iter.second->first_bit(b)) {
        do {
          dc.remote_callxs(b, push_receiver_function,
                            NULL,
                            0,
                            vertex,
                            graph.vertex_data(vertex));
        }while(iter.second->next_bit(b));
      }
    }
  }


/*************            Message Handlers          *************************/


  /**
  This the primary conduit for requesting for locks and data.
  \param requests: A vector of all the vertices and the lock level desired
  \param requestid: A user defined value. This value will be returned to the
                    caller when the lock is complete
  \param responsehandlerptr: A pointer to a lock response handler
  \param edgedatareqs: Edge data requests. The destination of these edges
                      should be in the requests list
  \param withdata:  Whether vertex data should be transmitted along with the lock
                  completion signal */
  static void request_lock_handler(distributed_control& dc, size_t source,
                              void* unused, size_t len,
                              std::vector<lockrequest> &requests,
                              size_t requestid,
                              size_t responsehandlerptr,
                              std::vector<edge_id_t> &edgedatareqs,
                              bool withdata) {
    lock_manager_target->request_lock_deferred(source, requestid, requests,
                                    (lock_response_handler)responsehandlerptr,
                                    edgedatareqs, withdata);
  }

  /// The distributed_lock_manager used by the handlers
  static distributed_lock_manager* lock_manager_target;

  /**
    Unlocks a bunch of requested locks. the locktypes should match the original
    request.
  */
  static void request_unlock_handler(distributed_control& dc, size_t source,
                              void* unused, size_t len,
                              std::vector<lockrequest> &requests) {
    lock_manager_target->release_lock(requests);
  }

  /**
    Unlocks a bunch of requested locks and commits some modified data at the same
    time. The locktypes should match the original request.
  */
  static void request_unlock_write_handler(distributed_control& dc, size_t source,
                    void* unused, size_t len,
                    std::vector<lockrequest> &requests,
                    std::map<vertex_id_t, typename Graph::vertex_data_type> &vdata,
                    std::map<edge_id_t, typename Graph::edge_data_type> &edata) {
    lock_manager_target->release_lock(requests, vdata, edata, source);
  }

 private:
  /***********************   Private Functions  **************************/
  // this function accesses the unsynchronized distrwlock
  // and should not be called directly
  inline bool lock(const vertex_id_t &v, const lock_type &locktype) {
    ASSERT_LT(v, locks.size());
    if (locktype == NOLOCK) return true;
    // if there is a writer. fail
    if (locks[v].writer) return false;
    switch(locktype) {
     case WRLOCK:
      // if there is a reader, fail
      if (locks[v].numreaders.value > 0) {
        return false;
      }
      else locks[v].writelock();
      break;
     case RDLOCK:
        locks[v].readlock();
       break;
     case NOLOCK:
      ASSERT_TRUE(locktype != NOLOCK);
    }
    return true;
  }

  // this function accesses the unsynchronized distrwlock
  // and should not be called directly
  inline void unlock(const vertex_id_t &v) {
    ASSERT_LT(v, locks.size());
    locks[v].unlock();
  }

  /**
    Called by the lock request handler. Creates a lock_decscriptor
    and puts it on the queue. */
  void request_lock_deferred(size_t sourcemachine, size_t requestid,
                              const std::vector<lockrequest> &lockset,
                              lock_response_handler responsehandler,
                              const std::vector<edge_id_t> &edgerequests,
                              bool withdata) {
    // create a lock set. This lock set is implicitly managed
    // within the lockqueue datastructure
    // fill the lockset datastructure
    lock_descriptor *newlockset = new lock_descriptor;
    newlockset->sourcemachine = sourcemachine;
    newlockset->requestid = requestid;
    newlockset->lockindex = 0;
    newlockset->lockrequests = lockset;
    newlockset->responsehandler = responsehandler;
    newlockset->withdata = withdata;
    newlockset->edgerequests = edgerequests;
    acquire_as_much_as_possible(*newlockset);
  }

  /** releases the locks in the lockrequest array. locktypes should match
  the original request. */
  void release_lock(const std::vector<lockrequest> &vertexlists) {
    for (int i = vertexlists.size() - 1;i >= 0; --i) {
      vertex_id_t vertextounlock = vertexlists[i].vertex;
      
      if (vertexlists[i].locktype == NOLOCK) continue;

      // get the lock number to use
      DASSERT_MSG(vertex2locknumber.find(vertextounlock) != vertex2locknumber.end(),
                  "vertex %d not on processor %d", vertextounlock, dc.procid());
      size_t locknum = vertex2locknumber[vertextounlock];
      

      lockqueue_locks[locknum].lock();
      // we should be always able to unlock
      unlock(locknum);
      
        // writers first!
      if (!writelockqueue[locknum].empty()) {
        // and we have writers pending?
        if (locks[locknum].numreaders.value == 0) {
          // if readers are 0 , and we have writers pending, 
          // give it up to a writer
          lock_descriptor* writer = writelockqueue[locknum].front();
          writelockqueue[locknum].pop_front();
          ASSERT_TRUE(lock(locknum, WRLOCK));
          lockqueue_locks[locknum].unlock();
          writer->lockindex++;
          acquire_as_much_as_possible(*writer);
          continue;
        }
        else {
          // readers are not yet zero. wait till all readers release
          lockqueue_locks[locknum].unlock(); 
          continue;
        }
      }

      // if we reach here, we should try the readers
      
      // if we reach here, we just have either
      // there are no writers and we should clear all the readers!
      if (!readlockqueue[locknum].empty()) {
        std::list<lock_descriptor*> tmp;
        tmp.swap(readlockqueue[locknum]);
        foreach(lock_descriptor* cur, tmp) {
          ASSERT_TRUE(lock(locknum, RDLOCK));
          cur->lockindex++;
        }
        lockqueue_locks[locknum].unlock();

        foreach(lock_descriptor* cur, tmp) {
          acquire_as_much_as_possible(*cur);
        }
        continue;
      }
      else {
        lockqueue_locks[locknum].unlock();
        continue;
      }

      
    }
    
   
  }
  
  
  
  /** releases the locks in the lockrequest array and committing some
  data at the same time. locktypes should match the original request. */
  void release_lock(const std::vector<lockrequest> &vertexlists,
                  std::map<vertex_id_t, typename Graph::vertex_data_type> &vdata,
                  std::map<edge_id_t, typename Graph::edge_data_type> &edata,
                  size_t source = -1) {
    // write back the edge
    {
      typename std::map<edge_id_t,
                        typename Graph::edge_data_type>::iterator it = edata.begin();
      while(it != edata.end()) {
        graph.edge_data(it->first) = it->second;
        set_dirty_edge(it->first, source);
        ++it;
      }
    }
    // write back the vertex edge
    {
      typename std::map<vertex_id_t,
                        typename Graph::vertex_data_type>::iterator it = vdata.begin();
      while(it != vdata.end()) {
//        logstream(LOG_INFO) << "proc "<< dc.procid() << " update vertex " << it->first << std::endl;
        graph.vertex_data(it->first) = it->second;
        set_dirty(it->first, source);
        ++it;
      }
    }
    release_lock(vertexlists);
  }


  /** returns true if any locks were acquired
      this tries to acquire as many locks as possible.
      putting on the back of the queue of the last lock it failed to acquire */
  bool acquire_as_much_as_possible(lock_descriptor &lockset) {
    // there should only be one
    lockset.lock.lock();
    bool successfullocks = false;
    bool lockcomplete = true;
    vertex_id_t prevvertex = lockset.lockindex > 0 ? lockset.lockrequests[lockset.lockindex - 1].vertex : 0;


    for (size_t i = lockset.lockindex; i < lockset.lockrequests.size(); ++i) {
      vertex_id_t vertextolock = lockset.lockrequests[i].vertex;
      ASSERT_LE(prevvertex, vertextolock);
      // get the lock number to use
      DASSERT_MSG(vertex2locknumber.find(vertextolock) != vertex2locknumber.end(),
                  "vertex %d not on processor %d", vertextolock, dc.procid());
      size_t locknum = vertex2locknumber[vertextolock];
      lock_type locktype = lockset.lockrequests[i].locktype;
      
      if (locktype == NOLOCK) {
        // automatic success. don't even bother locking
        lockset.lockindex = i + 1;
        successfullocks = true;
      }
      else {
        lockqueue_locks[locknum].lock();
        // if all queues are empty and the lock is available, we can take it now
        if (writelockqueue[locknum].empty() && 
           lock(locknum, locktype)) {
          // lock succeeded!
          lockset.lockindex = i + 1;
          successfullocks = true;
          lockqueue_locks[locknum].unlock();
        }
        else {
          // put index at the last failure and put it into the queue
          lockset.lockindex = i;

          if (locktype == RDLOCK) {
            readlockqueue[locknum].push_back(&lockset);
          }
          else {
            writelockqueue[locknum].push_back(&lockset);
          }
          lockcomplete = false;
          lockqueue_locks[locknum].unlock();
          break;
        }
      }
    }
    // If I have successfully locked everything
    if (lockcomplete) {
      // reply lock success deletes the lockset
      reply_lock_success(lockset);
    }
    else {
      lockset.lock.unlock();
    }
    return successfullocks;
  }

  /**
  Replies to the lock request. Called when all locks have been successfully
  acquired.
  */
  void reply_lock_success(lock_descriptor &lockset) {
    std::map<vertex_id_t, typename Graph::vertex_data_type> vdata;
    std::map<edge_id_t, typename Graph::edge_data_type> edata;
    dense_bitset &machinedirty = proc2dirty[lockset.sourcemachine];    
    dense_bitset &edgedirty = proc2dirtyedge[lockset.sourcemachine];    
    if (lockset.withdata) {
      // build the datamaps;

      for (size_t i = 0;i < lockset.lockrequests.size(); ++i) {
        vertex_id_t vertex = lockset.lockrequests[i].vertex;
          // if it is dirty, insert it
        if (!caching || machinedirty.clear_bit(vertex2locknumber[vertex])) {
          vdata[vertex] = graph.vertex_data(vertex);
        }
     }
    }

    for (size_t i = 0;i < lockset.edgerequests.size(); ++i) {
      edge_id_t edge = lockset.edgerequests[i];
      if (!caching || edgedirty.clear_bit(graph.global_to_local_eid(edge))) {
        edata[edge] = graph.edge_data(edge);
      }
    }

//     printf("success %d %d\n", lockset.sourcemachine, lockset.requestid);

    lockset.lock.unlock();
    dc.remote_callxs(lockset.sourcemachine, lockset.responsehandler, NULL,
                    0, lockset.requestid,vdata,edata);
    delete(&lockset);
  }
};

template<typename Graph>
distributed_lock_manager<Graph>* distributed_lock_manager<Graph>::lock_manager_target;
}

#include <graphlab/macros_undef.hpp>
#endif
