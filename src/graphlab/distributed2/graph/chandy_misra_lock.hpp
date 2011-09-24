#ifndef CHANDY_MISRA_LOCK_HPP
#define CHANDY_MISRA_LOCK_HPP

#include <boost/function.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/parallel/parallel_includes.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/lazy_deque.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/graph/graph_lock.hpp>
#include <graphlab/macros_def.hpp>

//#define DISTRIBUTED_LOCK_DEBUG

namespace graphlab {

  template <typename GraphType>
  class shm_chandy_misra_lock {
   public:
    typedef typename GraphType::edge_list_type edge_list_type;
    typedef typename GraphType::local_edge_list_type local_edge_list_type;
    typedef typename GraphType::vertex_id_type vertex_id_type;
    typedef typename GraphType::edge_id_type   edge_id_type;
    
    struct vertex_lock_state {
      mutex lock; // protects this strict
      bool trying_to_lock;
      bool eating;  // true if we are eating
      bool callback_sent;  // only send callback once until a cancel
      int pending;  // counts the number of forks I still have not acquired
      boost::function<void(vertex_id_type)> ready_callback;
      vertex_lock_state():trying_to_lock(false), 
                           eating(false), callback_sent(false),pending(0), ready_callback(NULL){ }
    };

   private:
    GraphType& dgraph;
    dense_bitset forks_and_state;
    
    // forks_and_state[3 * forkid] is true if source own the fork
    // forks_and_state[3 * forkid + 1] is true if the fork is dirty
    // forks_and_state[3 * forkid + 2] is true if the fork is requested by the remote party
    
    inline size_t FORK_ID_OWNER(size_t forkid) {
      return 3 * forkid;
    }

    inline size_t FORK_ID_DIRTY(size_t forkid) {
      return 3 * forkid + 1;
    }

    inline size_t FORK_ID_REQUESTED(size_t forkid) {
      return 3 * forkid + 2;
    }

    vertex_id_t get_fork_owner(size_t forkid) {
      if (forks_and_state.get(FORK_ID_OWNER(forkid))) return dgraph.get_local_store().source(forkid);
      else return dgraph.get_local_store().target(forkid);
    }
    
    bool fork_owner_is_source(size_t forkid) {
      return forks_and_state.get(FORK_ID_OWNER(forkid));
    }    
    
    void swap_fork_owner(size_t forkid) {
      if (forks_and_state.get(FORK_ID_OWNER(forkid))) forks_and_state.clear_bit(FORK_ID_OWNER(forkid));
      else forks_and_state.set_bit(FORK_ID_OWNER(forkid));
    }
    
    void set_fork_owner(size_t forkid, vertex_id_t vid) {
      if (dgraph.get_local_store().source(forkid) == vid) forks_and_state.set_bit(FORK_ID_OWNER(forkid));
      else forks_and_state.clear_bit(FORK_ID_OWNER(forkid));
    }
    

    
    std::vector<vertex_lock_state> vertex_state; 
   public:
    shm_chandy_misra_lock(GraphType &g):dgraph(g),
                                        forks_and_state(3 * dgraph.get_local_store().num_edges()), 
                                        vertex_state(dgraph.get_local_store().num_vertices()) { 
      for (size_t i = 0;i < dgraph.get_local_store().num_edges(); ++i) {
        //all forks are dirty
        forks_and_state.set_bit(FORK_ID_DIRTY(i));
        //forks owned by lower ID
        vertex_id_t source_lvid = dgraph.get_local_store().source(i);
        vertex_id_t target_lvid = dgraph.get_local_store().target(i);
    
        vertex_id_t source_gvid = dgraph.localvid_to_globalvid(source_lvid);
        vertex_id_t target_gvid = dgraph.localvid_to_globalvid(target_lvid);
        if (source_gvid < target_gvid) {
          forks_and_state.set_bit(FORK_ID_OWNER(i));
        }
      }
    }
    // requesting vertex localvid for fork on edge forkid.
    // requestor is the other side of the edge    
    void request_for_fork(vertex_id_type localvid,
                          size_t forkid) {
#ifdef DISTRIBUTED_LOCK_DEBUG
//      logstream(LOG_DEBUG) << "Requesting For fork " << forkid << " from " << localvid << std::endl;
#endif        
      ASSERT_LT(localvid, vertex_state.size());
      vertex_state[localvid].lock.lock();
      ASSERT_EQ(get_fork_owner(forkid), localvid);
      // if I do not have the fork, we wait. 
      // If the fork is clean, we wait.
      // if we are eating, we wait
      if (vertex_state[localvid].eating ||
          forks_and_state.get(FORK_ID_DIRTY(forkid)) == false) {
        forks_and_state.set_bit(FORK_ID_REQUESTED(forkid));
        vertex_state[localvid].lock.unlock();
      }
      else {
        // all other cases, we can give away the fork

        bool rerequest_fork = false;
        // if we are trying to lock, increment pending and send
        // a request together with the fork forwarding
        if (vertex_state[localvid].trying_to_lock) {
          vertex_state[localvid].pending++;
          ASSERT_GE(vertex_state[localvid].pending, 0);
          rerequest_fork = true;
        }
        // clean fork
        forks_and_state.clear_bit(FORK_ID_DIRTY(forkid));
        
        // disown the fork
        swap_fork_owner(forkid);

        vertex_state[localvid].lock.unlock();
        // give away fork
        forward_fork(localvid, forkid, rerequest_fork);
      }
    }
    
    void finalize(vertex_id_type localvid, 
                  boost::function<void(vertex_id_type, bool)> finalize_callback) {
      vertex_state[localvid].lock.lock();
      bool response = (vertex_state[localvid].pending == 0);
      if (response) {
        vertex_state[localvid].trying_to_lock = false;
        vertex_state[localvid].eating = true;
      }
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Eating vid: " << dgraph.localvid_to_globalvid(localvid) << std::endl;
#endif
      vertex_state[localvid].lock.unlock();
      finalize_callback(dgraph.localvid_to_globalvid(localvid), response);
      
    }
    
    void cancel_finalize(vertex_id_type localvid) {
      bool ready = false;
      boost::function<void(vertex_id_type)> cb = NULL;
      vertex_state[localvid].lock.lock();
      vertex_state[localvid].trying_to_lock = true;
      vertex_state[localvid].eating = false;
      vertex_state[localvid].callback_sent = false;
      std::vector<edge_id_type> forks_to_forward;
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Cacncelling Eating vid: " << dgraph.localvid_to_globalvid(localvid) << std::endl;
#endif
      // this is quite expensive...
      foreach(edge_id_type forkid, dgraph.get_local_store().in_edge_ids(localvid)) {
        if (fork_owner_is_source(forkid) == false &&
            forks_and_state.get(FORK_ID_DIRTY(forkid)) &&
            forks_and_state.get(FORK_ID_REQUESTED(forkid))) {
          // someone asked for the fork. ok. lets disown the fork
          vertex_state[localvid].pending++;
          // clean fork
          forks_and_state.clear_bit(FORK_ID_REQUESTED(forkid));
          forks_and_state.clear_bit(FORK_ID_DIRTY(forkid));
          swap_fork_owner(forkid);
          forks_to_forward.push_back(forkid);
        }
      }
      // again for the out edges
      foreach(edge_id_type forkid, dgraph.get_local_store().out_edge_ids(localvid)) {
        if (fork_owner_is_source(forkid) &&
            forks_and_state.get(FORK_ID_DIRTY(forkid)) &&
            forks_and_state.get(FORK_ID_REQUESTED(forkid))) {
          // someone asked for the fork. ok. lets disown the fork
          vertex_state[localvid].pending++;
          // clean fork
          forks_and_state.clear_bit(FORK_ID_REQUESTED(forkid));
          forks_and_state.clear_bit(FORK_ID_DIRTY(forkid));
          swap_fork_owner(forkid);
          forks_to_forward.push_back(forkid);
        }
      }
      
      if (vertex_state[localvid].pending == 0) {
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Ready vid: " << dgraph.localvid_to_globalvid(localvid) << std::endl;
#endif
        // ready
        vertex_state[localvid].callback_sent = true;
        ready = true;
        cb = vertex_state[localvid].ready_callback;
      }

      vertex_state[localvid].lock.unlock();
      if (forks_to_forward.size() > 0) {
        for (size_t i = 0;i < forks_to_forward.size(); ++i) {
          forward_fork(localvid, forks_to_forward[i], true);
        }
      }
      if (ready) {
        cb(dgraph.localvid_to_globalvid(localvid));
      }
      
    }

    /**
      Issues a request for all forks in the scope
      The ready() callback is called once all forks are available; though at this point
      forks may still disappear. Once the client receives ready() calls from all its requests,
      it should then issue a finalize() request for the vertex which will lock down all forks.
      If finalize() returns false from any machine, the client must issue cancel_finalize()
      to all machines and wait for ready() again.
      
      Client pseudocode
      \code
      for s in all_scope_requests {
        issue_request_for_forks_in_scope(s, ready_callback);
      }

      while(1) {
        wait for all scope requests to return ready
        for s in all_scope_requests {
          finalize(s, finalize_callback)
        }
        wait for all finalize_callbacks
        if (not all finalize_callbacks return true) {
          for s in all_scope_requests {
            cancel_finalize(s)
          }        
        }
      }
      \endcode
    */
    void issue_request_for_forks_in_scope(vertex_id_type localvid, 
                                          boost::function<void(vertex_id_type)> ready_callback) {
      ASSERT_LT(localvid, vertex_state.size());
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Locking vid: " << dgraph.localvid_to_globalvid(localvid) << std::endl;
#endif
      vertex_state[localvid].lock.lock();
      vertex_state[localvid].trying_to_lock = true;
      vertex_state[localvid].eating = false;
      vertex_state[localvid].callback_sent = false;
      vertex_state[localvid].ready_callback = ready_callback;
      std::vector<std::pair<vertex_id_type, edge_id_type> > requests;
      foreach(edge_id_type eid, dgraph.get_local_store().in_edge_ids(localvid)) {
        // if we do not have the fork, we need to ask for it
        if (fork_owner_is_source(eid)) {
          requests.push_back(std::make_pair(dgraph.get_local_store().source(eid), eid));
        }
      }
      foreach(edge_id_type eid, dgraph.get_local_store().out_edge_ids(localvid)) {
        // if we do not have the fork, we need to ask for it
        if (!fork_owner_is_source(eid)) {
          requests.push_back(std::make_pair(dgraph.get_local_store().target(eid), eid));
        }
      }
      vertex_state[localvid].pending = requests.size();
      ASSERT_GE(vertex_state[localvid].pending, 0);
      if (vertex_state[localvid].pending == 0) {
        vertex_state[localvid].callback_sent = true;
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Ready vid: " << dgraph.localvid_to_globalvid(localvid) << std::endl;
#endif
        vertex_state[localvid].lock.unlock();
        // all locks are mine
        ready_callback(dgraph.localvid_to_globalvid(localvid));
      }
      else {
        // there are requests we need to fulfil 
        // note that the way it works is that without a fork request, we will never
        // get the forks. Therfore it is safe to release the lock here
        // then just loop over the request set
        vertex_state[localvid].lock.unlock();
        for (size_t i = 0;i < requests.size(); ++i) {
          request_for_fork(requests[i].first, requests[i].second);
        }
      }
    } 


    void release_forks_in_scope(vertex_id_type localvid) {
      ASSERT_LT(localvid, vertex_state.size());
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Releasing vid: " << dgraph.localvid_to_globalvid(localvid) << std::endl;
#endif

      vertex_state[localvid].lock.lock();
      ASSERT_TRUE(vertex_state[localvid].eating);
      ASSERT_FALSE(vertex_state[localvid].trying_to_lock);
      vertex_state[localvid].eating = false;
      vertex_state[localvid].ready_callback = NULL;
      std::vector<edge_id_type> forks_to_forward;
      foreach(edge_id_type forkid, dgraph.get_local_store().in_edge_ids(localvid)) {
        if (forks_and_state.get(FORK_ID_REQUESTED(forkid))) {
          // someone asked for the fork. ok. lets disown the fork
          // clean fork
          forks_and_state.clear_bit(FORK_ID_REQUESTED(forkid));
          forks_and_state.clear_bit(FORK_ID_DIRTY(forkid));
          swap_fork_owner(forkid);
          forks_to_forward.push_back(forkid);
        }
        else {
          // all used forks are dirty
          forks_and_state.set_bit(FORK_ID_DIRTY(forkid));
        }
      }
      // again for the out edges
      foreach(edge_id_type forkid, dgraph.get_local_store().out_edge_ids(localvid)) {
        if (forks_and_state.get(FORK_ID_REQUESTED(forkid))) {
          // someone asked for the fork. ok. lets disown the fork
          // clean fork
          forks_and_state.clear_bit(FORK_ID_REQUESTED(forkid));
          forks_and_state.clear_bit(FORK_ID_DIRTY(forkid));
          swap_fork_owner(forkid);
          forks_to_forward.push_back(forkid);
        }
        else {
          // all used forks are dirty
          forks_and_state.set_bit(FORK_ID_DIRTY(forkid));
        }
      }      
      vertex_state[localvid].lock.unlock();
      
      for (size_t i = 0;i < forks_to_forward.size(); ++i) {
        forward_fork(localvid, forks_to_forward[i], false);
      }

    }

    // Give fork 'forkid' to vertex localvid. If rerequest_fork is set
    // the fork is wanted by the other party
    void incoming_fork(vertex_id_t localvid, size_t forkid, bool rerequest_fork) {
      ASSERT_LT(localvid, vertex_state.size());
#ifdef DISTRIBUTED_LOCK_DEBUG
//      logstream(LOG_DEBUG) << "Fork " << localvid << " received by " << localvid << std::endl;
#endif

      bool ready = false;
      boost::function<void(vertex_id_type)> cb = NULL;
      vertex_state[localvid].lock.lock();
      // we now have the fork. make sure we actually do have the fork
      ASSERT_EQ(get_fork_owner(forkid), localvid);
      if (rerequest_fork) forks_and_state.set_bit(FORK_ID_REQUESTED(forkid));
      else forks_and_state.clear_bit(FORK_ID_REQUESTED(forkid));
      vertex_state[localvid].pending--;
      if (vertex_state[localvid].pending == 0) {
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Vertex " << dgraph.localvid_to_globalvid(localvid) << " Eating" << std::endl;
#endif
        // ready
        if (!vertex_state[localvid].callback_sent) {
          ready = true;
          vertex_state[localvid].callback_sent = true;
          cb = vertex_state[localvid].ready_callback;
        }
      }
#ifdef DISTRIBUTED_LOCK_DEBUG
      else {
        logstream(LOG_DEBUG) << "Vertex " << dgraph.localvid_to_globalvid(localvid) << ": " << vertex_state[localvid].pending << " Remaining" << std::endl;
      }
#endif

      vertex_state[localvid].lock.unlock();
      if (ready) {
        cb(dgraph.localvid_to_globalvid(localvid));
      }
    }


    // send fork to opposing party of the vid.
    void forward_fork(vertex_id_t localvid, size_t forkid, bool rerequest_from) {
      ASSERT_LT(localvid, vertex_state.size());
      vertex_id_t src_localvid = dgraph.get_local_store().source(forkid);
      vertex_id_t target_localvid = dgraph.get_local_store().target(forkid);
      
      vertex_id_t local_altvid = src_localvid!=localvid ? src_localvid : target_localvid;
      // fork should be cleaned
      ASSERT_FALSE(forks_and_state.get(FORK_ID_DIRTY(forkid)));
#ifdef DISTRIBUTED_LOCK_DEBUG
//      logstream(LOG_DEBUG) << "Tx Fork " << forkid << ": " << localvid << "-->" << local_altvid << std::endl;
#endif
      incoming_fork(local_altvid, forkid, rerequest_from);
    }


    
    void print_state() {
      for (size_t i = 0;i < vertex_state.size(); ++i) {
        std::cout << "vtx " << dgraph.localvid_to_globalvid(i) << ":" << vertex_state[i].trying_to_lock << " "
                           << vertex_state[i].eating << " "
                           << vertex_state[i].pending << " " 
                           << vertex_state[i].callback_sent << "\n";
                           
        foreach(edge_id_type eid, dgraph.get_local_store().in_edge_ids(i)) {
          std::cout << "\t" << dgraph.get_local_store().source(eid) << "->" << dgraph.get_local_store().target(eid)  << ": " 
                             << (get_fork_owner(eid) == i) << "," 
                            << (forks_and_state.get(FORK_ID_DIRTY(eid))?"D":"-" )
                            << (forks_and_state.get(FORK_ID_REQUESTED(eid))?"R":"-") << "\n";
        }
        foreach(edge_id_type eid, dgraph.get_local_store().out_edge_ids(i)) {
          std::cout << "\t" << dgraph.get_local_store().source(eid) << "->" << dgraph.get_local_store().target(eid)  << ": " 
                             << (get_fork_owner(eid) == i) << "," 
                            << (forks_and_state.get(FORK_ID_DIRTY(eid))?"D":"-") 
                            << (forks_and_state.get(FORK_ID_REQUESTED(eid))?"R":"-") << "\n";
        } 
      }
    }
  };


  /**
  An implementation of the Chandy/Misra solution to the dining philosophers 
  problem. Only works on edge and vertex scopes.
  
  \note
  This class is templatized over graph type, but really it requires relatively deep
  introspection into the graph and will only work with the distributed_graph implementation
  right now. This could be generalized.
  */
  template <typename GraphType>
  class chandy_misra_lock : public graph_lock<GraphType>{
   public:
    typedef typename GraphType::edge_list_type edge_list_type;
    typedef typename GraphType::local_edge_list_type local_edge_list_type;
    typedef typename GraphType::vertex_id_type vertex_id_type;
    typedef typename GraphType::edge_id_type   edge_id_type;
    
    
    struct vertex_lock_state {
      mutex lock; // protects this strict
      size_t pending;
      size_t finalize_count;
      bool has_failed_finalize;
      bool locked;
      bool locking;
      boost::function<void(vertex_id_type)> callback;
      vertex_lock_state():pending(0), has_failed_finalize(false), locked(false),locking(false),callback(NULL) { }
    };
    
   private:
    GraphType& dgraph;
    shm_chandy_misra_lock<GraphType> cmlock;
    dc_dist_object<chandy_misra_lock<GraphType> > rmi;
    bool synchronize_data;
    
    std::vector<vertex_lock_state> vertex_state; 

    // requesting for fork on edge forkid where t
    
   public:    
    void print_state() {
      cmlock.print_state();
      for (size_t i = 0;i < vertex_state.size(); ++i) {
        std::cout << "vtx " << dgraph.localvid_to_globalvid(i) << (vertex_state[i].locked?"L":"-")
                 << (vertex_state[i].locking?"K":"-")
                  << ":" << vertex_state[i].pending  << " "
                           << vertex_state[i].finalize_count << " "
                           << vertex_state[i].has_failed_finalize << std::endl;

      }
    }
    
    void local_scope_ready_callback(vertex_id_t globalvid, procid_t source) {
    #ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Local scope completion for "<< globalvid << " from " << source << std::endl;
    #endif
      if (source == rmi.procid()) {
        remote_scope_ready(globalvid);
      }
      else {
        rmi.remote_call(source, &chandy_misra_lock<GraphType>::remote_scope_ready, globalvid);
      }
    }
   
    void local_scope_request(vertex_id_t globalvid, procid_t source) {
    #ifdef DISTRIBUTED_LOCK_DEBUG
          logstream(LOG_DEBUG) << "Local scope request for "<< globalvid << " from " << source << std::endl;
    #endif
      cmlock.issue_request_for_forks_in_scope(dgraph.globalvid_to_localvid(globalvid),
                                      boost::bind(&chandy_misra_lock<GraphType>::local_scope_ready_callback,
                                                  this, _1, source));
    }
    
    void local_cancel_finalize(vertex_id_t globalvid) {
      cmlock.cancel_finalize(dgraph.globalvid_to_localvid(globalvid));
    
    }
    
    void local_finalize_reply_callback(vertex_id_t globalvid, bool ready, procid_t source) {
      if (source == rmi.procid()) {
        remote_finalize_reply_callback(globalvid, ready);
      }
      else {
        rmi.remote_call(source, &chandy_misra_lock<GraphType>::remote_finalize_reply_callback, globalvid, ready);
      }
    }
    
    void remote_finalize_reply_callback(vertex_id_t globalvid, bool ready) {
      size_t localvid = dgraph.globalvid_to_localvid(globalvid);
      bool done = false;
      bool hasfailed = false;
      vertex_state[localvid].lock.lock();
      ASSERT_GT(vertex_state[localvid].finalize_count, 0);
      vertex_state[localvid].finalize_count--;
      if (ready == false) {
        vertex_state[localvid].has_failed_finalize = true;
      }
      if (vertex_state[localvid].finalize_count == 0) {
        hasfailed = vertex_state[localvid].has_failed_finalize;
        if (hasfailed) {
          // revert back to the ready state
          vertex_state[localvid].pending = dgraph.localvid_to_replicas(localvid).size();
        }
        done = true;
      }
      vertex_state[localvid].lock.unlock();
      if (done) {
        if (hasfailed) {
          #ifdef DISTRIBUTED_LOCK_DEBUG
                logstream(LOG_DEBUG) << "Finalize Failed "<< globalvid << std::endl;
          #endif
          char prevkey = rmi.dc().set_sequentialization_key((globalvid % 254) + 1);
          const std::vector<procid_t>& procs = dgraph.localvid_to_replicas(localvid);
          for (size_t i = 0;i < procs.size(); ++i) {
            if (procs[i] == rmi.procid()) {
              local_cancel_finalize(globalvid);
            }
            else {
              rmi.remote_call(procs[i], 
                              &chandy_misra_lock<GraphType>::local_cancel_finalize,
                              globalvid);
            }
          }
          rmi.dc().set_sequentialization_key(prevkey);
        }
        else {
          #ifdef DISTRIBUTED_LOCK_DEBUG
                logstream(LOG_DEBUG) << "Finalize Success "<< globalvid << std::endl;
          #endif
          if (synchronize_data && dgraph.on_boundary(globalvid)) {
            // ok... we need to do this one more time.
            dgraph.async_synchronize_scope_callback(globalvid, 
                                                    boost::bind(&chandy_misra_lock<GraphType>::lock_completion, this, localvid));
          }
          else {
            lock_completion(localvid);
          }
        }
      }
    }
    
    void local_finalize_request(vertex_id_t globalvid, procid_t source) {
    #ifdef DISTRIBUTED_LOCK_DEBUG
          logstream(LOG_DEBUG) << "Local finalize request for "<< globalvid << " from " << source << std::endl;
    #endif
      cmlock.finalize(dgraph.globalvid_to_localvid(globalvid),
                      boost::bind(&chandy_misra_lock<GraphType>::local_finalize_reply_callback,
                                  this, _1, _2, source));

    }
    
    void remote_scope_ready(vertex_id_t globalvid) {
      bool done = false;
      size_t localvid = dgraph.globalvid_to_localvid(globalvid);
      vertex_state[localvid].lock.lock();
      ASSERT_GT(vertex_state[localvid].pending, 0);
      vertex_state[localvid].pending--;
      if (vertex_state[localvid].pending == 0) {
        done = true;
        vertex_state[localvid].finalize_count = dgraph.localvid_to_replicas(localvid).size();
        vertex_state[localvid].has_failed_finalize = false;
      }
      vertex_state[localvid].lock.unlock();
      
      if (done) {
        // issue finalize
        const std::vector<procid_t>& procs = dgraph.localvid_to_replicas(localvid);
        char prevkey = rmi.dc().set_sequentialization_key((globalvid % 254) + 1);
        for (size_t i = 0;i < procs.size(); ++i) {
          if (procs[i] == rmi.procid()) {
            local_finalize_request(globalvid, rmi.procid());
          }
          else {
            rmi.remote_call(procs[i], 
                            &chandy_misra_lock<GraphType>::local_finalize_request, 
                            globalvid, rmi.procid());
          }
        }
        rmi.dc().set_sequentialization_key(prevkey);
      }
    }
    
    void lock_completion(vertex_id_t localvid) {
    #ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "Lock completion for "<< dgraph.localvid_to_globalvid(localvid) << std::endl;
    #endif
      vertex_state[localvid].lock.lock();
      boost::function<void(vertex_id_type)> cb = vertex_state[localvid].callback;
      vertex_state[localvid].callback = NULL;
      vertex_state[localvid].locked = true;
      vertex_state[localvid].locking = false;
      vertex_state[localvid].lock.unlock();
      cb(dgraph.localvid_to_globalvid(localvid));
    }
    
    chandy_misra_lock(distributed_control &dc,
                      GraphType &dgraph, 
                      bool synchronize_data = false):dgraph(dgraph),cmlock(dgraph), 
                                       rmi(dc, this), 
                                       synchronize_data(synchronize_data),
                                       vertex_state(dgraph.get_local_store().num_vertices()){ }

    /**
       Requests a lock on the scope surrounding globalvid.
       This globalvid must be owned by the current machine.
       When lock is complete the handler is called.
    */
    void scope_request(vertex_id_type globalvid,
                       boost::function<void(vertex_id_type)> handler,
                       scope_range::scope_range_enum scopetype) {
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "scope request for "<< globalvid << std::endl;
#endif
      ASSERT_TRUE(scopetype == scope_range::EDGE_CONSISTENCY);
      // broadcast a request to all ghosts
      vertex_id_type localvid = dgraph.globalvid_to_localvid(globalvid);
      vertex_state[localvid].lock.lock();
      vertex_state[localvid].callback = handler;
      const std::vector<procid_t>& procs = dgraph.localvid_to_replicas(localvid);
      vertex_state[localvid].pending = procs.size();
      ASSERT_FALSE(vertex_state[localvid].locked);
      ASSERT_FALSE(vertex_state[localvid].locking);
      vertex_state[localvid].locked = false;
      vertex_state[localvid].locking = true;
      
      char prevkey = rmi.dc().set_sequentialization_key((globalvid % 254) + 1);
      vertex_state[localvid].lock.unlock();
      
      for (size_t i = 0;i < procs.size(); ++i) {
        if (procs[i] == rmi.procid()) {
          local_scope_request(globalvid, rmi.procid());
        }
        else {
          rmi.remote_call(procs[i], 
                          &chandy_misra_lock<GraphType>::local_scope_request, 
                          globalvid, rmi.procid());
        }
      }
      
      rmi.dc().set_sequentialization_key(prevkey);

    }


    void local_scope_unlock(vertex_id_t globalvid) {
      vertex_id_type localvid = dgraph.globalvid_to_localvid(globalvid);
      cmlock.release_forks_in_scope(localvid);
    }
   
    /**
       Isues an unlock on the scope surrounding globalvid.
       A lock on this scope MUST have been acquired before or
       very bad things will happen
    */
    void scope_unlock(vertex_id_type globalvid,
                      scope_range::scope_range_enum scopetype) {
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "scope release for "<< globalvid << std::endl;
#endif
      ASSERT_TRUE(scopetype == scope_range::EDGE_CONSISTENCY);
      
      vertex_id_type localvid = dgraph.globalvid_to_localvid(globalvid);
      
      
      unsigned char prevkey = 0;
      bool synchronizing = false;
      prevkey = rmi.dc().set_sequentialization_key((globalvid % 254) + 1);
      if (synchronize_data && dgraph.on_boundary(globalvid)) {
        synchronizing = true;
        dgraph.synchronize_scope(globalvid, true);
      }

      vertex_state[localvid].lock.lock();
      ASSERT_FALSE(vertex_state[localvid].locking);
      ASSERT_TRUE(vertex_state[localvid].locked);
      vertex_state[localvid].locked = false;
      const std::vector<procid_t>& procs = dgraph.localvid_to_replicas(localvid);
      ASSERT_EQ(vertex_state[localvid].pending, 0);
      vertex_state[localvid].lock.unlock();
      
      for (size_t i = 0;i < procs.size(); ++i) {
        if (procs[i] == rmi.procid()) {
          local_scope_unlock(globalvid);
        }
        else {
          rmi.remote_call(procs[i],&chandy_misra_lock<GraphType>::local_scope_unlock, globalvid);
        }
      }      
      rmi.dc().set_sequentialization_key(prevkey);
    }
  };
  
}

#include <graphlab/macros_undef.hpp>
#undef DISTRIBUTED_LOCK_DEBUG
#endif
