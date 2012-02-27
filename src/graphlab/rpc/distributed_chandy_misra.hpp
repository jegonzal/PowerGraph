#ifndef GRAPHLAB_LOCAL_CHANDY_MISRA_HPP
#define GRAPHLAB_LOCAL_CHANDY_MISRA_HPP
#include <vector>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/macros_def.hpp>
namespace graphlab {

template <typename GraphType>
class distributed_chandy_misra {
 public:
  typedef typename GraphType::local_graph_type local_graph_type;
  typedef typename local_graph_type::edge_type edge_type;
  typedef typename GraphType::vertex_record vertex_record;
  typedef distributed_chandy_misra<GraphType> dcm_type;
  dc_dist_object<dcm_type> rmi;
  GraphType &distgraph;
  local_graph_type &graph;
  boost::function<void(vertex_id_type)> callback;
  /*
   * Each "fork" is one character.
   * bit 0: owner. if 0 is src. if 1 is target
   * bit 1: clean = 0, dirty = 1
   * bit 2: owner 0 request
   * bit 3: owner 1 request
   */
  std::vector<unsigned char> forkset;
  enum { OWNER_BIT = 1,
         DIRTY_BIT = 2,
         REQUEST_0 = 4,
         REQUEST_1 = 8 };
  enum {OWNER_SOURCE = 0, OWNER_TARGET = 1};
  inline unsigned char request_bit(bool owner) {
    return owner ? REQUEST_1 : REQUEST_0;
  }

  struct philosopher {
    vertex_id_type num_edges;
    vertex_id_type forks_acquired;
    simple_spinlock lock;
    unsigned char state;
    unsigned char counter;
    bool cancellation_sent;
    bool lockid;
  };
  std::vector<philosopher> philosopherset;
  /*
   * Possible values for the philosopher state
   */
  enum {
    THINKING = 0,
    HUNGRY = 1,
    HORS_DOEUVRE = 2, 
    EATING = 3
  };

  /** Places a request for the fork. Requires fork to be locked */
  inline void request_for_fork(size_t forkid, bool nextowner) {
    forkset[forkid] |= request_bit(nextowner);
  }

  inline bool fork_owner(size_t forkid) {
    return forkset[forkid] & OWNER_BIT;
  }

  inline bool fork_dirty(size_t forkid) {
    return !!(forkset[forkid] & DIRTY_BIT);
  }

  inline void dirty_fork(size_t forkid) {
    forkset[forkid] |= DIRTY_BIT;
  }


  void compute_initial_fork_arrangement() {
    for (vertex_id_type i = 0;i < graph.num_vertices(); ++i) {
      philosopherset[i].num_edges = graph.num_in_edges(i) +
                                    graph.num_out_edges(i);
      philosopherset[i].state = THINKING;
      philosopherset[i].forks_acquired = 0;
      philosopherset[i].counter = 0;
      philosopherset[i].cancellation_sent = false;
      philosopherset[i].lockid = false;
    }
    for (vertex_id_type i = 0;i < graph.num_vertices(); ++i) {
      foreach(edge_type edge, graph.in_edges(i)) {
        if (distgraph.global_vid(edge.source()) >
            distgraph.global_vid(edge.target())) {
          forkset[graph.edge_id(edge)] = DIRTY_BIT | OWNER_TARGET;
          philosopherset[edge.target()].forks_acquired++;
        }
        else {
          forkset[graph.edge_id(edge)] = DIRTY_BIT | OWNER_SOURCE;
          philosopherset[edge.source()].forks_acquired++;
        }
      }
    }
  }

  /**
   * We already have v1, we want to acquire v2.
   * When this function returns, both v1 and v2 locks are acquired
   */
  void try_acquire_edge_with_backoff(vertex_id_type v1,
                                     vertex_id_type v2) {
    if (v1 < v2) {
      philosopherset[v2].lock.lock();
    }
    else if (!philosopherset[v2].lock.try_lock()) {
        philosopherset[v1].lock.unlock();
        philosopherset[v2].lock.lock();
        philosopherset[v1].lock.lock();
    }
  }
  
/****************************************************************************
 * Tries to move a requested fork
 *
 * Pseudocode:
 *  If current owner is hungry and fork is clean
 *    Ignore
 *  ElseIf current owner is Thinking
 *    Relinquish fork immediately and clear the request flag
 *  ElseIf current owner is hors_doeuvre and fork is clean
 *    Ignore
 *  ElseIf current owner is hors_doeuvre and fork is dirty
 *    Send cancellation message
 *    Set cancelsent
 *  End
 * Return true if changes were made
 ***************************************************************************/
  inline bool advance_fork_state_on_lock(size_t forkid,
                                        vertex_id_type source,
                                        vertex_id_type target) {
    unsigned char currentowner = forkset[forkid] & OWNER_BIT;
    if (currentowner == OWNER_SOURCE) {
      // if the current owner is not eating, and the
      // fork is dirty and other side has placed a request
      if (philosopherset[source].state != EATING &&
          (forkset[forkid] & DIRTY_BIT) &&
          (forkset[forkid] & REQUEST_1)) {

        if (philosopherset[source].state != HORS_DOEUVRE) {
          //  change the owner and clean the fork)
          forkset[forkid] = OWNER_TARGET;
          if (philosopherset[source].state == HUNGRY) {
            forkset[forkid] |= REQUEST_0;
          }
          philosopherset[source].forks_acquired--;
          philosopherset[target].forks_acquired++;
          return true;
        }
        else if (philosopherset[source].cancellation_sent == false) {
          philosopherset[source].cancellation_sent = true;
          bool lockid = philosopherset[source].lockid;
          philosopherset[source].lock.unlock();
          philosopherset[target].lock.unlock();
          issue_cancellation_request_unlocked(source, lockid);
          philosopherset[std::min(source, target)].lock.lock();
          philosopherset[std::max(source, target)].lock.lock();
        }
      }
    }
    else {
      // if the current owner is not eating, and the
      // fork is dirty and other side has placed a request
      if (philosopherset[target].state != EATING &&
          (forkset[forkid] & DIRTY_BIT) &&
          (forkset[forkid] & REQUEST_0)) {
        //  change the owner and clean the fork)
        if (philosopherset[target].state != HORS_DOEUVRE) {
          forkset[forkid] = OWNER_SOURCE;
          if (philosopherset[target].state == HUNGRY) {
            forkset[forkid] |= REQUEST_1;
          }
          philosopherset[source].forks_acquired++;
          philosopherset[target].forks_acquired--;
          return true;
        }
        else if (philosopherset[target].cancellation_sent == false) {
          philosopherset[target].cancellation_sent = true;
          bool lockid = philosopherset[target].lockid;
          philosopherset[source].lock.unlock();
          philosopherset[target].lock.unlock();
          issue_cancellation_request_unlocked(target, lockid);
          philosopherset[std::min(source, target)].lock.lock();
          philosopherset[std::max(source, target)].lock.lock();
        }
      }
    }
    return false;
  }


/****************************************************************************
 * Performs a cancellation on a vertex.
 * 
 * If lockIds do not match, ignore
 * If counter == 0 ignore
 * Otherwise, counter++ and reply cancellation accept.
 * Unfortunately, I cannot perform a local call here even if I am the
 * owner since this may produce a lock cycle. Irregardless of whether
 * the owner is local or not, this must be performed by a remote call
 ***************************************************************************/

  void cancellation_request_unlocked(vertex_id_type lvid, procid_t requestor, bool lockid) {
    philosopherset[lvid].lock.lock();
    
    if (philosopherset[lvid].lockid == lockid) {
      if (philosopherset[lvid].counter > 0) {
        ++philosopherset[lvid].counter;
        bool lockid = philosopherset[lvid].lockid;
        logstream(LOG_DEBUG) << rmi.procid() <<
            ": Cancellation accepted on " << distgraph.global_vid(lvid) <<
            "(" << (int)philosopherset[lvid].counter << ")" << std::endl;
        philosopherset[lvid].lock.unlock();
        const vertex_record& rec = distgraph.l_get_vertex_record(lvid);
        if (requestor != rmi.procid()) {
          unsigned char pkey = rmi.dc().set_sequentialization_key(rec.gvid % 254 + 1);
          rmi.remote_call(requestor,
                          &dcm_type::rpc_cancellation_accept,
                          rec.gvid,
                          lockid);
          rmi.dc().set_sequentialization_key(pkey);
        }
        else {
          cancellation_accept_unlocked(lvid, lockid);
        }
      }
      else {
        philosopherset[lvid].lock.unlock();
        logstream(LOG_DEBUG) << rmi.procid() <<
          ": Cancellation on " << distgraph.global_vid(lvid) <<
          " denied due to lock completion" << std::endl;
      }
    }
    else {
      philosopherset[lvid].lock.unlock();
      logstream(LOG_DEBUG) << rmi.procid() <<
        ": Cancellation on " << distgraph.global_vid(lvid) <<
        " denied to invalid lock ID" << std::endl;
    }
    
  }

  void rpc_cancellation_request(vertex_id_type gvid, procid_t requestor, bool lockid) {
    vertex_id_type lvid = distgraph.local_vid(gvid);
    cancellation_request_unlocked(lvid, requestor, lockid);
  }

  void issue_cancellation_request_unlocked(vertex_id_type lvid, bool lockid) {
    // signal the master
    logstream(LOG_DEBUG) << rmi.procid() <<
        ": Requesting cancellation on " << distgraph.global_vid(lvid) << std::endl;
    const vertex_record &rec = distgraph.l_get_vertex_record(lvid);
    
    if (rec.owner == rmi.procid()) {
      cancellation_request_unlocked(lvid, rmi.procid(), lockid);
    }
    else {
      unsigned char pkey = rmi.dc().set_sequentialization_key(rec.gvid % 254 + 1);
      rmi.remote_call(rec.owner,
                      &dcm_type::rpc_cancellation_request,
                      rec.gvid,
                      rmi.procid(), 
                      lockid);
      rmi.dc().set_sequentialization_key(pkey);

    }
  }


/****************************************************************************
 * Accepts a cancellation on a vertex.
 *
 * Pseudocode
 *  Change back to Hungry
 *  Releases all dirty forks
 ****************************************************************************/

  void rpc_cancellation_accept(vertex_id_type gvid, bool lockid) {
    vertex_id_type lvid = distgraph.local_vid(gvid);
    cancellation_accept_unlocked(lvid, lockid);
  }

  void cancellation_accept_unlocked(vertex_id_type p_id, bool lockid) {
    std::vector<vertex_id_type> retval;
    philosopherset[p_id].lock.lock();
    //philosopher is now hungry!
    ASSERT_EQ (lockid, philosopherset[p_id].lockid);
    ASSERT_EQ((int)philosopherset[p_id].state, (int)HORS_DOEUVRE);
    philosopherset[p_id].state = HUNGRY;
    philosopherset[p_id].cancellation_sent = false;
    logstream(LOG_DEBUG) << rmi.procid() <<
            ": Cancellation accept received on " << distgraph.global_vid(p_id) << " " <<
            philosopherset[p_id].state << std::endl;

    // for each fork I own, try to give it away
    foreach(edge_type edge, graph.in_edges(p_id)) {
      try_acquire_edge_with_backoff(edge.target(), edge.source());
      //std::cout << "\t" << graph.edge_id(edge) << ": " << edge.source() << "->" << edge.target() << std::endl;
      vertex_id_type other = edge.source();
      size_t edgeid = graph.edge_id(edge);
      if (fork_owner(edgeid) == OWNER_TARGET && fork_dirty(edgeid)) {
        
        if (advance_fork_state_on_lock(edgeid, edge.source(), edge.target()) &&
            philosopherset[other].state == HUNGRY &&
            philosopherset[other].forks_acquired == philosopherset[other].num_edges) {
          philosopherset[other].state = HORS_DOEUVRE;
          // signal eating on other
          retval.push_back(other);
        }
      }
      philosopherset[edge.source()].lock.unlock();
    }
    //std::cout << "out edges: " << std::endl;
    foreach(edge_type edge, graph.out_edges(p_id)) {
      //std::cout << "\t" << graph.edge_id(edge) << ": " << edge.source() << "->" << edge.target() << std::endl;
      try_acquire_edge_with_backoff(edge.source(), edge.target());

      vertex_id_type other = edge.target();
      size_t edgeid = graph.edge_id(edge);
      if (fork_owner(edgeid) == OWNER_SOURCE && fork_dirty(edgeid)) {
        if (advance_fork_state_on_lock(edgeid, edge.source(), edge.target()) &&
            philosopherset[other].state == HUNGRY &&
            philosopherset[other].forks_acquired == philosopherset[other].num_edges) {
          philosopherset[other].state = HORS_DOEUVRE;
          // signal eating on other
          retval.push_back(other);
        }
      }
      philosopherset[edge.target()].lock.unlock();
    }
    
    if (philosopherset[p_id].state == HUNGRY &&
        philosopherset[p_id].forks_acquired == philosopherset[p_id].num_edges) {
      retval.push_back(p_id);
    }
      
    philosopherset[p_id].lock.unlock();
    foreach(vertex_id_type lvid, retval) {
      enter_hors_doeuvre_unlocked(lvid);
    }

  }
  
/****************************************************************************
 * Make Philosopher Hungry.
 *
 * Pseudocode:
 * Set Philosopher to Hungry
 * For all edges adjacent to v with forks it does not own:
 *   Send request for fork to neighboring vertex
 *
 * Conditions:
 *   Must be Thinking
 *   New lock ID must not be the same as the old lock ID
 *
 * Possible Immediate Transitions:
 *   Current vertex may enter HORS_DOEUVRE
 ***************************************************************************/
  void rpc_make_philosopher_hungry(vertex_id_type gvid, bool newlockid) {
    vertex_id_type lvid = distgraph.local_vid(gvid);
    logstream(LOG_DEBUG) << rmi.procid() <<
          ": Local HUNGRY Philosopher  " << gvid << std::endl;
    philosopherset[lvid].lock.lock();

    ASSERT_EQ((int)philosopherset[lvid].state, (int)THINKING);
    philosopherset[lvid].state = HUNGRY;

    ASSERT_NE(philosopherset[lvid].lockid, newlockid);
    philosopherset[lvid].lockid = newlockid;

    philosopherset[lvid].lock.unlock();

    local_philosopher_grabs_forks(lvid);
  }

  void local_philosopher_grabs_forks(vertex_id_type p_id) {
    philosopherset[p_id].lock.lock();
    //philosopher is now hungry!
// now try to get all the forks. lock one edge at a time
    // using the backoff strategy
    //std::cout << "vertex " << p_id << std::endl;
    //std::cout << "in edges: " << std::endl;
    foreach(edge_type edge, graph.in_edges(p_id)) {
      try_acquire_edge_with_backoff(edge.target(), edge.source());
      //std::cout << "\t" << graph.edge_id(edge) << ": " << edge.source() << "->" << edge.target() << std::endl;
      size_t edgeid = graph.edge_id(edge);
      // if fork is owned by other edge, try to take it
      if (fork_owner(edgeid) == OWNER_SOURCE) {
        request_for_fork(edgeid, OWNER_TARGET);
        advance_fork_state_on_lock(edgeid, edge.source(), edge.target());
      }
      philosopherset[edge.source()].lock.unlock();
    }
    //std::cout << "out edges: " << std::endl;
    foreach(edge_type edge, graph.out_edges(p_id)) {
      //std::cout << "\t" << graph.edge_id(edge) << ": " << edge.source() << "->" << edge.target() << std::endl;
      try_acquire_edge_with_backoff(edge.source(), edge.target());
      size_t edgeid = graph.edge_id(edge);

      // if fork is owned by other edge, try to take it
      if (fork_owner(edgeid) == OWNER_TARGET) {
        request_for_fork(edgeid, OWNER_SOURCE);
        advance_fork_state_on_lock(edgeid, edge.source(), edge.target());
      }
      philosopherset[edge.target()].lock.unlock();
    }

    bool enter_hors = false;
    if (philosopherset[p_id].state == HUNGRY &&
        philosopherset[p_id].forks_acquired == philosopherset[p_id].num_edges) {
      philosopherset[p_id].state = HORS_DOEUVRE;
      enter_hors = true;
    }
    philosopherset[p_id].lock.unlock();
    if (enter_hors) enter_hors_doeuvre_unlocked(p_id);
  }

/************************************************************************
 *
 * Called when a vertex may be ready to enter hors dourre
 * Locks must be maintained. HORS_DOEUVRE must be set prior
 * to entering this function .
 *
 ***********************************************************************/
  void enter_hors_doeuvre_unlocked(vertex_id_type p_id) {
     // if I got all forks I can eat
    logstream(LOG_DEBUG) << rmi.procid() <<
            ": Local HORS_DOEUVRE Philosopher  " << distgraph.global_vid(p_id) << std::endl;
    // signal the master
    const vertex_record &rec = distgraph.l_get_vertex_record(p_id);
    if (rec.owner == rmi.procid()) {
      signal_ready_unlocked(p_id, philosopherset[p_id].lockid);
    }
    else {
      unsigned char pkey = rmi.dc().set_sequentialization_key(rec.gvid % 254 + 1);
      rmi.remote_call(rec.owner,
                      &dcm_type::rpc_signal_ready,
                      rec.gvid, philosopherset[p_id].lockid);
      rmi.dc().set_sequentialization_key(pkey);
    }
  }

/************************************************************************
 *
 * Called when a vertex enters HORS_DOEUVRE. Locks must be maintained.
 * 
 * Conditions:
 *   vertex must be in HUNGRY or HORS_DOEUVRE
 *   lock IDs must match
 *
 * Possible Immediate Transitions
 *   If counter == 0, transit to EATING
 ***********************************************************************/

  void signal_ready_unlocked(vertex_id_type lvid, bool lockid) {
    philosopherset[lvid].lock.lock();
    if(!(philosopherset[lvid].state == (int)HUNGRY ||
      philosopherset[lvid].state == (int)HORS_DOEUVRE)) {
      logstream(LOG_ERROR) << rmi.procid() <<
              ": Bad signal ready state!!!! : " << (int)philosopherset[lvid].state << std::endl;
      logstream(LOG_ERROR) << rmi.procid() <<
              " Lock IDs : " << (int)philosopherset[lvid].lockid << " " << (int)lockid << std::endl;
      logstream(LOG_ERROR) << rmi.procid() <<
            ": BAD Global HORS_DOEUVRE " << distgraph.global_vid(lvid)
            << "(" << (int)philosopherset[lvid].counter << ")" << std::endl;
  
      ASSERT_TRUE(philosopherset[lvid].state == (int)HUNGRY ||
                  philosopherset[lvid].state == (int)HORS_DOEUVRE);
    }
    
    ASSERT_EQ(philosopherset[lvid].lockid, lockid);
    philosopherset[lvid].counter--;

    logstream(LOG_DEBUG) << rmi.procid() <<
            ": Global HORS_DOEUVRE " << distgraph.global_vid(lvid)
            << "(" << (int)philosopherset[lvid].counter << ")" << " " << (int)(philosopherset[lvid].state) << std::endl;
  
    if(philosopherset[lvid].counter == 0) {
      philosopherset[lvid].lock.unlock();
      // broadcast EATING
      const vertex_record& rec = distgraph.l_get_vertex_record(lvid);
      unsigned char pkey = rmi.dc().set_sequentialization_key(rec.gvid % 254 + 1);
      rmi.remote_call(rec.mirrors().begin(), rec.mirrors().end(), 
                      &dcm_type::rpc_set_eating, rec.gvid, lockid);
      set_eating(lvid, lockid);
      rmi.dc().set_sequentialization_key(pkey);
    }
    else {
      philosopherset[lvid].lock.unlock();
    }
  }


  void rpc_signal_ready(vertex_id_type gvid, bool lockid) {
    vertex_id_type lvid = distgraph.local_vid(gvid);
    signal_ready_unlocked(lvid, lockid);
  }

  void set_eating(vertex_id_type lvid, bool lockid) {
    philosopherset[lvid].lock.lock();
    
    logstream(LOG_DEBUG) << rmi.procid() <<
            ": EATING " << distgraph.global_vid(lvid)
            << "(" << (int)philosopherset[lvid].counter << ")" << std::endl;
  
    ASSERT_EQ((int)philosopherset[lvid].state, (int)HORS_DOEUVRE);
    ASSERT_EQ(philosopherset[lvid].lockid, lockid);
    philosopherset[lvid].state = EATING;
    philosopherset[lvid].cancellation_sent = false;
    philosopherset[lvid].lock.unlock();
    if (distgraph.l_get_vertex_record(lvid).owner == rmi.procid()) {
    logstream(LOG_DEBUG) << rmi.procid() <<
            ": CALLBACK " << distgraph.global_vid(lvid) << std::endl;

      callback(lvid);
    }
  }

  void rpc_set_eating(vertex_id_type gvid, bool lockid) {

    logstream(LOG_DEBUG) << rmi.procid() <<
            ": Receive Set EATING " << gvid << std::endl;
  
    vertex_id_type lvid = distgraph.local_vid(gvid);
    set_eating(lvid, lockid);
  }
/************************************************************************
 *
 * Called when a vertex stops eating
 *
 ***********************************************************************/



  inline bool advance_fork_state_on_unlock(size_t forkid,
                                         vertex_id_type source,
                                         vertex_id_type target) {

    unsigned char currentowner = forkset[forkid] & OWNER_BIT;
    if (currentowner == OWNER_SOURCE) {
      // if the current owner is not eating, and the
      // fork is dirty and other side has placed a request
      if ((forkset[forkid] & DIRTY_BIT) &&
        (forkset[forkid] & REQUEST_1)) {
        //  change the owner and clean the fork)
        // keep my request bit if any
        forkset[forkid] = OWNER_TARGET;
        philosopherset[source].forks_acquired--;
        philosopherset[target].forks_acquired++;
        return true;
      }
    }
    else {
      // if the current owner is not eating, and the
      // fork is dirty and other side has placed a request
      if ((forkset[forkid] & DIRTY_BIT) &&
        (forkset[forkid] & REQUEST_0)) {
        //  change the owner and clean the fork)
        // keep my request bit if any
        forkset[forkid] = OWNER_SOURCE;
        philosopherset[source].forks_acquired++;
        philosopherset[target].forks_acquired--;
        return true;
      }
    }
    return false;
  }


  

  
  void local_philosopher_stops_eating(size_t p_id) {
    std::vector<vertex_id_type> retval;
    philosopherset[p_id].lock.lock();
    if (philosopherset[p_id].state != EATING) {
      std::cout << rmi.procid() << ": " << p_id << "FAILED!! Cannot Stop Eating!" << std::endl;
      ASSERT_EQ((int)philosopherset[p_id].state, (int)EATING);
    }
    philosopherset[p_id].state = THINKING;
    philosopherset[p_id].counter = 0;

    // now forks are dirty
    foreach(edge_type edge, graph.in_edges(p_id)) {
      try_acquire_edge_with_backoff(edge.target(), edge.source());
      size_t edgeid = graph.edge_id(edge);
      vertex_id_type other = edge.source();
      dirty_fork(edgeid);
      advance_fork_state_on_unlock(edgeid, edge.source(), edge.target());
      if (philosopherset[other].state == HUNGRY &&
            philosopherset[other].forks_acquired ==
                philosopherset[other].num_edges) {
        philosopherset[other].state = HORS_DOEUVRE;
        // signal eating on other
        retval.push_back(other);
      }
      philosopherset[other].lock.unlock();
    }

    foreach(edge_type edge, graph.out_edges(p_id)) {
      try_acquire_edge_with_backoff(edge.source(), edge.target());
      size_t edgeid = graph.edge_id(edge);
      vertex_id_type other = edge.target();
      dirty_fork(edgeid);
      advance_fork_state_on_unlock(edgeid, edge.source(), edge.target());
      if (philosopherset[other].state == HUNGRY &&
            philosopherset[other].forks_acquired ==
                philosopherset[other].num_edges) {
        philosopherset[other].state = HORS_DOEUVRE;
        // signal eating on other
        retval.push_back(other);
      }
      philosopherset[other].lock.unlock();
    }

    philosopherset[p_id].lock.unlock();
    foreach(vertex_id_type lvid, retval) {
      enter_hors_doeuvre_unlocked(lvid);
    }
  }

  void rpc_philosopher_stops_eating(vertex_id_type gvid) {
    logstream(LOG_DEBUG) << rmi.procid() << ": Receive STOP eating on " << gvid << std::endl;
    local_philosopher_stops_eating(distgraph.local_vid(gvid));
  }

 public:
  inline distributed_chandy_misra(distributed_control &dc,
                                  GraphType &distgraph,
                                  boost::function<void(vertex_id_type)> callback
                                  ):
                          rmi(dc, this),
                          distgraph(distgraph),
                          graph(distgraph.get_local_graph()),
                          callback(callback){
    forkset.resize(graph.num_edges(), 0);
    philosopherset.resize(graph.num_vertices());
    compute_initial_fork_arrangement();
    rmi.barrier();
  }

  inline const vertex_id_type invalid_vid() const {
    return (vertex_id_type)(-1);
  }

  void initialize_master_philosopher_as_hungry_locked(vertex_id_type p_id,
                                                      bool lockid) {
    philosopherset[p_id].lockid = lockid;
    philosopherset[p_id].state = HUNGRY;
    philosopherset[p_id].counter = 
      distgraph.l_get_vertex_record(p_id).num_mirrors() + 1;
  }
  
  void make_philosopher_hungry(vertex_id_type p_id) {
    const vertex_record &rec = distgraph.l_get_vertex_record(p_id);
    ASSERT_EQ(rec.get_owner(), rmi.procid());
    philosopherset[p_id].lock.lock();
    ASSERT_EQ((int)philosopherset[p_id].state, (int)THINKING);
    bool newlockid = !philosopherset[p_id].lockid;
    initialize_master_philosopher_as_hungry_locked(p_id, newlockid);
    
    logstream(LOG_DEBUG) << rmi.procid() <<
            ": Global HUNGRY " << distgraph.global_vid(p_id)
            << "(" << (int)philosopherset[p_id].counter << ")" << std::endl;
  
    philosopherset[p_id].lock.unlock();
    
    unsigned char pkey = rmi.dc().set_sequentialization_key(rec.gvid % 254 + 1);
    rmi.remote_call(rec.mirrors().begin(), rec.mirrors().end(), 
                    &dcm_type::rpc_make_philosopher_hungry, rec.gvid, newlockid);
    rmi.dc().set_sequentialization_key(pkey);
    local_philosopher_grabs_forks(p_id);
  }
  
  void philosopher_stops_eating(vertex_id_type p_id) {
    const vertex_record &rec = distgraph.l_get_vertex_record(p_id);
    ASSERT_EQ(rec.get_owner(), rmi.procid());


    logstream(LOG_DEBUG) << rmi.procid() <<
            ": Global STOP Eating " << distgraph.global_vid(p_id) << std::endl;

    philosopherset[p_id].lock.lock();
    ASSERT_EQ(philosopherset[p_id].state, (int)EATING);
    philosopherset[p_id].counter = 0;
    philosopherset[p_id].lock.unlock();
    unsigned char pkey = rmi.dc().set_sequentialization_key(rec.gvid % 254 + 1);
    rmi.remote_call(rec.mirrors().begin(), rec.mirrors().end(), 
                    &dcm_type::rpc_philosopher_stops_eating, rec.gvid);
    rmi.dc().set_sequentialization_key(pkey);
    local_philosopher_stops_eating(p_id);
  }
  
  void no_locks_consistency_check() {
    // make sure all forks are dirty
    for (size_t i = 0;i < forkset.size(); ++i) ASSERT_TRUE(fork_dirty(i));
    // all philosophers are THINKING
    for (size_t i = 0;i < philosopherset.size(); ++i) ASSERT_TRUE(philosopherset[i].state == THINKING);
  }

  void print_out() {
    std::cout << "Philosophers\n";
    std::cout << "------------\n";
    for (vertex_id_type v = 0; v < graph.num_vertices(); ++v) {
      std::cout << distgraph.global_vid(v) << ": " << (int)philosopherset[v].state << " " << (int)philosopherset[v].counter << "\n";
    }
  }
  
  void complete_consistency_check() {
    for (vertex_id_type v = 0; v < graph.num_vertices(); ++v) {
      // count the number of forks I own
      size_t numowned = 0;
      size_t numowned_clean = 0;
      foreach(edge_type edge, graph.in_edges(v)) {
        size_t edgeid = graph.edge_id(edge);
        if (fork_owner(edgeid) == OWNER_TARGET) {
          numowned++;
          if (!fork_dirty(edgeid)) numowned_clean++;
        }
      }
      foreach(edge_type edge, graph.out_edges(v)) {
        size_t edgeid = graph.edge_id(edge);
        if (fork_owner(edgeid) == OWNER_SOURCE) {
          numowned++;
          if (!fork_dirty(edgeid)) numowned_clean++;
        }
      }

      ASSERT_EQ(philosopherset[v].forks_acquired, numowned);
      if (philosopherset[v].state == THINKING) {
        ASSERT_EQ(numowned_clean, 0);
      }
      else if (philosopherset[v].state == HUNGRY) {
        ASSERT_NE(philosopherset[v].num_edges, philosopherset[v].forks_acquired);
        // any fork I am unable to acquire. Must be clean, and the other person
        // must be eating or hungry
        foreach(edge_type edge, graph.in_edges(v)) {
          size_t edgeid = graph.edge_id(edge);
          // not owned
          if (fork_owner(edgeid) == OWNER_SOURCE) {
            if (philosopherset[edge.source()].state != EATING) {
              if (fork_dirty(edgeid)) {
                std::cout << (int)(forkset[edgeid]) << " "
                          << (int)philosopherset[edge.source()].state
                          << "->" << (int)philosopherset[edge.target()].state
                          << std::endl;
                ASSERT_FALSE(fork_dirty(edgeid));
              }
            }
            ASSERT_NE(philosopherset[edge.source()].state, (int)THINKING);
          }
        }
        foreach(edge_type edge, graph.out_edges(v)) {
          size_t edgeid = graph.edge_id(edge);
          if (fork_owner(edgeid) == OWNER_TARGET) {
            if (philosopherset[edge.target()].state != EATING) {
              if (fork_dirty(edgeid)) {
                std::cout << (int)(forkset[edgeid]) << " "
                          << (int)philosopherset[edge.source()].state
                          << "->"
                          << (int)philosopherset[edge.target()].state
                          << std::endl;
                ASSERT_FALSE(fork_dirty(edgeid));
              }
            }
            ASSERT_NE(philosopherset[edge.target()].state, (int)THINKING);
          }
        }

      }
      else if (philosopherset[v].state == EATING) {
        ASSERT_EQ(philosopherset[v].forks_acquired, philosopherset[v].num_edges);
      }
    }
  }
};

}

#include <graphlab/macros_undef.hpp>

#endif
