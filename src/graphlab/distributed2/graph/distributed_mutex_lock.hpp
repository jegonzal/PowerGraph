/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


#ifndef GRAPHLAB_distributed_mutex_lock_HPP
#define GRAPHLAB_distributed_mutex_lock_HPP
#include <list>
#include <boost/function.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/parallel/deferred_rwlock.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/lazy_deque.hpp>
#include <graphlab/distributed2/graph/distributed_graph.hpp>
#include <graphlab/distributed2/graph/distributed_mutex_lock.hpp>
namespace graphlab {

  namespace gl_impl {
    void distributed_mutex_lock_callback(size_t ptr, vertex_id_t vid);
  }

  // #define COMPILER_WRITE_BARRIER asm volatile("":::"memory")
#define COMPILER_WRITE_BARRIER
  //#define DISTRIBUTED_LOCK_DEBUG
  /**

     The locking implementation is basically two families of continuations.
     The first family is called the scopelock_continuation
     This family completes the lock of a scope.
     It iterates over the owners of the replicas of the vertex, and issue remote
     calls to acquire locks on them.
  
     The second family is called partiallock_continuation
     It completes the lock on local vertices.
     It iterates over the owned vertices within the scope of the vertex, acquiring
     locks.
  
     \note
     This class is templatized over graph type, but really it requires relatively deep
     introspection into the graph and will only work with the distributed_graph implementation
     right now. This could be generalized.
  */
  template <typename GraphType>
  class distributed_mutex_lock : public graph_lock<GraphType>{
  public:
    typedef typename GraphType::edge_list_type edge_list_type;
    typedef typename GraphType::local_edge_list_type local_edge_list_type;
    typedef typename GraphType::vertex_id_type vertex_id_type;
    typedef typename GraphType::edge_id_type   edge_id_type;

    distributed_mutex_lock(distributed_control &dc,
               GraphType &dgraph, bool strict_scope = false, 
                bool synchronize_data = false):dgraph(dgraph), 
                                                rmi(dc, this), 
                                                synchronize_data(synchronize_data),
                                                strict_scope(strict_scope){
      locks.resize(dgraph.owned_vertices().size());
    }

    void print_state() { }

    /**
       Requests a lock on the scope surrounding globalvid.
       This globalvid must be owned by the current machine.
       When lock is complete the handler is called.
    */
    void scope_request(vertex_id_type globalvid,
                       boost::function<void(vertex_id_type)> handler,
                       scope_range::scope_range_enum scopetype,
                       bool priority) {
#ifdef DISTRIBUTED_LOCK_DEBUG
      logstream(LOG_DEBUG) << "scope request for "<< globalvid << std::endl;
#endif

      // construct the scope lock parameters
      scopelock_cont_params sparams;
      sparams.globalvid = globalvid;
    
      sparams.localvid = dgraph.globalvid_to_localvid(globalvid);
      sparams.nextowneridx = 0;
      dgraph.localvid_to_replicas(sparams.localvid).first_bit(sparams.nextowneridx);
      // the number of syncs I need to receive. This will be #ghosts (popcount - 1) + 1.
      // For simplicity, require the final caller to always decrement by one. Makes it 
      // easier to figure out if I have received everything or not
      sparams.num_syncs_pending.value = synchronize_data?dgraph.localvid_to_replicas(sparams.localvid).popcount():1;
      sparams.neighbor_dirty_bloom = dgraph.get_owned_scope_dirty_bloom_filter(globalvid);
      sparams.handler = handler;
      sparams.scopetype = scopetype;
      sparams.priority = priority;
      scopelock_lock.lock();
      typename lazy_deque<scopelock_cont_params>::value_type* 
        ptr = scopelock_continuation.push_anywhere(sparams);
      scopelock_lock.unlock();
    
      continue_scope_lock(ptr, 1);
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

      // check if I need to unlock neighbors
      if (adjacent_vertex_lock_type(scopetype) != scope_range::NO_LOCK) {
    
        unsigned char prevkey = rmi.dc().set_sequentialization_key((globalvid % 254) + 1);
        if (synchronize_data && dgraph.on_boundary(globalvid)) {
          dgraph.synchronize_scope(globalvid, true);
        }
        
        // the complicated case. I need to unlock on my neighbors
        // get all my replicas

        vertex_id_type localvid = dgraph.globalvid_to_localvid(globalvid);

        const fixed_dense_bitset<MAX_N_PROCS>& procs = dgraph.localvid_to_replicas(localvid);
        uint32_t p = 0;
        ASSERT_TRUE(procs.first_bit(p));
        do {
          if (p == rmi.procid()) {
            partial_unlock(globalvid, (size_t)scopetype);
          }
          else {
            rmi.remote_call(p,
                            &distributed_mutex_lock<GraphType>::partial_unlock,
                            globalvid,
                            (size_t)scopetype);
          }
        } while(procs.next_bit(p));
        rmi.dc().set_sequentialization_key(prevkey);
      }
      else {
        partial_unlock(globalvid, (size_t)scopetype);
      } 
    }

        
    /**
       The parameters passed on to the scope lock continuation
    */
    struct scopelock_cont_params {
      vertex_id_type globalvid;
      vertex_id_type localvid;
      uint32_t nextowneridx;
      uint64_t neighbor_dirty_bloom;
      atomic<uint32_t> num_syncs_pending;
      scope_range::scope_range_enum scopetype;
      boost::function<void(vertex_id_type)> handler;
      bool priority;
    };

    /**
       The parameters passed on to the partial lock continuation
    */  
    struct partiallock_cont_params {
      size_t inidx;       // the next in idx to consider in the in_edges/out_edges parallel iteration
      size_t outidx;      // the next out idx to consider in the in_edges/out_edges parallel iteration
      vertex_id_type localvid;
      procid_t srcproc;
      uint64_t neighbor_dirty_bloom;
      uint32_t num_unsent_syncs;
      fixed_dense_bitset<MAX_N_PROCS> procset;
      __attribute__((may_alias)) size_t src_tag;   // holds a pointer to the caller's scope lock continuation
      scope_range::scope_range_enum scopetype;
      bool curlocked;
      bool priority;
      deferred_rwlock::request req;
    } __attribute__ ((aligned (8))) ;
  
  private:
    /// The distributed graph we are locking over
    GraphType &dgraph;
    /// The RMI object
    dc_dist_object<distributed_mutex_lock<GraphType> > rmi;
  
    /** the set of deferred locks local to this machine.
     * lock i corresponds to local vertex i. (by construction, 
     * the owned vertices always come first in the local store)
     */
    std::vector<deferred_rwlock> locks;

    mutex scopelock_lock;
    lazy_deque<scopelock_cont_params> scopelock_continuation;
    mutex partiallock_lock;
    lazy_deque<partiallock_cont_params> partiallock_continuation;

    bool synchronize_data;
    bool strict_scope;
    
    /**
       partial lock request on the sending processor
       Requests a lock on the scope surrounding the vertex globalvid 
       on some destination processor. This call completes a lock which is
       purely local to the destination processor.
       This globalvid must be in the fragment of the destination processor, 
       (either owned or a ghost). When locks have been acquired the handler is
       called.
    */
    void partial_lock_request(procid_t srcproc,
                              procid_t destproc,
                              vertex_id_type globalvid,
                              uint64_t neighbor_dirty_bloom,
                              uint32_t num_unsent_syncs,
                              const fixed_dense_bitset<MAX_N_PROCS> &procset,
                              scope_range::scope_range_enum scopetype,
                              bool priority,
                              size_t scope_continuation_ptr) {
      // here I issue to call to the remote machine, 
      // but I must pass on enough information so that I can call the handler
      // on reply
      // fast track it if the destination processor is local
      if (destproc == rmi.procid()) {
        // fast track! If it is local, just call it directly
        partial_lock_request_impl(plock_req_impl_arguments(srcproc, 
                                  globalvid,
                                  neighbor_dirty_bloom,
                                  num_unsent_syncs),
                                  procset,
                                  (size_t)scopetype, 
                                  priority,
                                  scope_continuation_ptr);
      }
      else {
        rmi.remote_call(destproc,
                        &distributed_mutex_lock<GraphType>::partial_lock_request_impl,
                        plock_req_impl_arguments(srcproc,
                        globalvid,
                        neighbor_dirty_bloom,
                        num_unsent_syncs),
                        procset,
                        (size_t)scopetype,
                        priority,
                        scope_continuation_ptr);
      }
      
    }


    void partial_lock_completion(size_t scope_continuation_ptr, uint32_t num_unsent_syncs) {
      typename lazy_deque<scopelock_cont_params>::value_type* 
        ptr = (typename lazy_deque<scopelock_cont_params>::value_type*)(scope_continuation_ptr);    
#ifdef DISTRIBUTED_LOCK_DEBUG    
      logstream(LOG_DEBUG) << "Receiving successful remote lock of " << ptr->first.globalvid << std::endl;
#endif
      continue_scope_lock(ptr, num_unsent_syncs);
    }

    void data_synchronize_reply(typename lazy_deque<scopelock_cont_params>::value_type* ptr) {
      scopelock_cont_params& params = ptr->first;
      params.handler(params.globalvid);
      // finish the continuation by erasing the lazy_deque entry
      scopelock_lock.lock();
      scopelock_continuation.erase(ptr);
      scopelock_lock.unlock();
    }
  
    void push_synchronization(size_t scope_continuation_ptr, const std::string& data) {
      typename lazy_deque<scopelock_cont_params>::value_type* 
        ptr = (typename lazy_deque<scopelock_cont_params>::value_type*)(scope_continuation_ptr);    
        scopelock_cont_params& params = ptr->first;
        dgraph.receive_external_update(data);
        
        if (synchronize_data == false || 
            params.num_syncs_pending.dec() == 0) {
          
          params.handler(params.globalvid);
          scopelock_lock.lock();
          scopelock_continuation.erase(ptr);
          scopelock_lock.unlock();
        }

    }

  
    void continue_scope_lock(typename lazy_deque<scopelock_cont_params>::value_type* ptr,
                             uint32_t num_unsent_syncs) {
      // for convenience, lets take a reference to the params
      scopelock_cont_params& params = ptr->first;
      // check if I need to actually lock my replicas
      // do not need to if the adjacent lock type is no lock
      if (adjacent_vertex_lock_type(params.scopetype) != scope_range::NO_LOCK) {
        // the complicated case. I need to lock on my neighbors

        const fixed_dense_bitset<MAX_N_PROCS>& procs = dgraph.localvid_to_replicas(params.localvid);
        uint32_t curidx = params.nextowneridx;
        if (curidx < rmi.numprocs()) {
          // I am done after this
          params.nextowneridx = rmi.numprocs();
          // process procs[curidx]
          // send a lock request, setting myself as the continuation
          partial_lock_request(rmi.procid(),
                               curidx,
                               params.globalvid,
                               params.neighbor_dirty_bloom,
                               num_unsent_syncs,
                               procs,
                               params.scopetype,
                               params.priority,
                               (size_t)(ptr));
        }
        else {
          // finish the continuation by erasing the lazy_deque entry
          // if synchronize data is set, issue one more continuation which
          // goes to a global fuctnction
          if (synchronize_data == false || 
              params.num_syncs_pending.dec(num_unsent_syncs) == 0) {
            
            params.handler(params.globalvid);
            scopelock_lock.lock();
            scopelock_continuation.erase(ptr);
            scopelock_lock.unlock();
          }
        }
      }
      else {
        // this is the easy case. I only need to lock on myself.
        // first check if I am actually in fact. done.
        if (params.nextowneridx == 0) {
          // no I am not done
          ++params.nextowneridx;
          // issue a partial lock request to to the current machine
          fixed_dense_bitset<MAX_N_PROCS> db;
          db.set_bit_unsync(rmi.procid());
          partial_lock_request(rmi.procid(),
                               rmi.procid(),
                               params.globalvid,
                               params.neighbor_dirty_bloom,
                               num_unsent_syncs,
                               db,
                               params.scopetype,
                               params.priority,
                               (size_t)(ptr));
        }
        else {
          // finish the continuation by erasing the lazy_deque entry
          // if synchronize data is set, issue one more continuation which
          // goes to a global fuctnction
          if (synchronize_data == false || 
              params.num_syncs_pending.dec(num_unsent_syncs) == 0) {
            
            params.handler(params.globalvid);
            scopelock_lock.lock();
            scopelock_continuation.erase(ptr);
            scopelock_lock.unlock();
          }
        }
      }
    }



    struct plock_req_impl_arguments{
      vertex_id_type globalvid;
      uint64_t neighbor_dirty_bloom;
      procid_t srcproc;
      uint32_t num_unsent_syncs;
      
      plock_req_impl_arguments() { }
      plock_req_impl_arguments(vertex_id_type globalvid,
                               uint64_t neighbor_dirty_bloom,
                               procid_t srcproc,
                               uint32_t num_unsent_syncs):
                                  globalvid(globalvid),
                                  neighbor_dirty_bloom(neighbor_dirty_bloom),
                                  srcproc(srcproc),
                                  num_unsent_syncs(num_unsent_syncs) { };

      
      inline void save(oarchive &oarc) const {      
        serialize(oarc, this, sizeof(plock_req_impl_arguments));
      }

      inline void load(iarchive &iarc) {
        deserialize(iarc, this, sizeof(plock_req_impl_arguments));
      }
    };
    

    /**
       lock request implementation on the receiving processor
    */
    void partial_lock_request_impl(plock_req_impl_arguments args,
                                   const fixed_dense_bitset<MAX_N_PROCS> &procset,
                                   size_t scopetype,
                                   bool priority,
                                   size_t src_tag) {
      procid_t srcproc = args.srcproc;
      vertex_id_type globalvid = args.globalvid;
      uint64_t neighbor_dirty_bloom = args.neighbor_dirty_bloom;
      uint32_t num_unsent_syncs = args.num_unsent_syncs;
                                   
      // construct a partiallock_continuation
#ifdef DISTRIBUTED_LOCK_DEBUG
      //  logstream(LOG_DEBUG) << rmi.procid() << ": p-lock request from "<< srcproc << " : " << globalvid << std::endl;
#endif
      partiallock_cont_params plockparams;
      plockparams.srcproc = srcproc;
      plockparams.inidx = 0;
      plockparams.outidx = 0;
      plockparams.priority = priority;
      plockparams.procset = procset;
      plockparams.procset.clear_bit_unsync(rmi.procid());
      plockparams.neighbor_dirty_bloom = neighbor_dirty_bloom;
      plockparams.num_unsent_syncs = num_unsent_syncs;
      plockparams.src_tag = src_tag;
      plockparams.curlocked = false;
      plockparams.req.next = NULL;
      plockparams.scopetype = (scope_range::scope_range_enum)(scopetype);
      // if no lock needed on adjacent vertices
      // set inidx and outidx to infty
      if (adjacent_vertex_lock_type(plockparams.scopetype) == scope_range::NO_LOCK) {
        plockparams.inidx = (vertex_id_type)(-1);
        plockparams.outidx = (vertex_id_type)(-1);
      }
    
      plockparams.localvid = dgraph.globalvid_to_localvid(globalvid);
      dgraph.get_local_store().set_vertex_dirty(plockparams.localvid, true);
      // put it inside the partiallock continuation
      partiallock_lock.lock();
      typename lazy_deque<partiallock_cont_params>::value_type* 
        ptr = partiallock_continuation.push_anywhere(plockparams);
      partiallock_lock.unlock();
      //sqeeze the pointer into WORD_SIZE - 2 bits.
      // essentially assuming a minimum of 4 byte alignment
      assert(((size_t)(ptr) & 3) == 0);
      ptr->first.req.id = (size_t)(ptr) >> 2;
      continue_partial_lock(ptr);
    }


    void continue_partial_lock(typename lazy_deque<partiallock_cont_params>::value_type* ptr) {
      partiallock_cont_params& params = ptr->first;
      // perform a parallel iteration across in and out edges

      vertex_id_type curv = params.localvid;

      local_edge_list_type inedges = dgraph.get_local_store().in_edge_ids(curv);
      local_edge_list_type outedges = dgraph.get_local_store().out_edge_ids(curv);
    
      vertex_id_type inv  = 
        (inedges.size() > params.inidx) ? 
        dgraph.get_local_store().source(inedges[params.inidx]) : vertex_id_type(-1);
      vertex_id_type outv  = 
        (outedges.size() > params.outidx) ? 
        dgraph.get_local_store().target(outedges[params.outidx]) : vertex_id_type(-1);
      // iterate both in order and lock
      // include the current vertex in the iteration
      while (params.inidx < inedges.size() || 
             params.outidx < outedges.size()) {
        if (!params.curlocked && curv < inv  && curv < outv) {
          // if curv is smaller than inv and outv
          params.curlocked = true;
          COMPILER_WRITE_BARRIER;
          // acquire current vertex
          if (dgraph.localvid_is_ghost(curv) == false &&
              issue_deferred_lock(curv, 
                                  params.req, 
                                  params.priority,
                                  central_vertex_lock_type(params.scopetype)) == false) {
            return;
          }
        } 
        else if (inv < outv) {
          ++params.inidx;
          COMPILER_WRITE_BARRIER;
          if (dgraph.localvid_is_ghost(inv) == false &&
              issue_deferred_lock(inv, 
                                  params.req, 
                                  params.priority,
                                  adjacent_vertex_lock_type(params.scopetype)) == false) {
            return;
          }
          inv  = (inedges.size() > params.inidx) ? dgraph.get_local_store().source(inedges[params.inidx]) : (vertex_id_type)(-1);
        } 
        else if (outv < inv) {
          ++params.outidx;
          COMPILER_WRITE_BARRIER;
          if (dgraph.localvid_is_ghost(outv) == false &&
              issue_deferred_lock(outv, 
                                  params.req, 
                                  params.priority,
                                  adjacent_vertex_lock_type(params.scopetype)) == false) {
            return;
          }
          outv  = (outedges.size() > params.outidx) ? dgraph.get_local_store().target(outedges[params.outidx]) : (vertex_id_type)(-1);
        } 
        else if (inv == outv){
          ++params.inidx; ++params.outidx;
          COMPILER_WRITE_BARRIER;
          if (dgraph.localvid_is_ghost(outv) == false &&
              issue_deferred_lock(outv, 
                                  params.req, 
                                  params.priority,
                                  adjacent_vertex_lock_type(params.scopetype)) == false) {
            return;
          }
          inv  = (inedges.size() > params.inidx) ? dgraph.get_local_store().source(inedges[params.inidx]) : (vertex_id_type)(-1);
          outv  = (outedges.size() > params.outidx) ? dgraph.get_local_store().target(outedges[params.outidx]) : (vertex_id_type)(-1);
        }
      }
      // just in case we never got around to locking it
      if (!params.curlocked) {
        params.curlocked = true;
        COMPILER_WRITE_BARRIER;
        // acquire current vertex
        if (dgraph.localvid_is_ghost(curv) == false &&
            issue_deferred_lock(curv, 
                                params.req, 
                                params.priority,
                                central_vertex_lock_type(params.scopetype)) == false) {
          return;
        }
      }
    
      // resolve the synchronization
      if (params.srcproc != rmi.procid() && synchronize_data) {
        std::string ret = dgraph.external_push_ghost_scope_to_owner(dgraph.localvid_to_globalvid(params.localvid),
                                                                    params.neighbor_dirty_bloom);
        if (ret.size() > 0) {
          rmi.remote_call(params.srcproc,
                         &distributed_mutex_lock<GraphType>::push_synchronization,
                         params.src_tag,
                         ret);
        }
        else {
          params.num_unsent_syncs++;
        }
      }
      
      // if we get here, the lock is complete on this machine
      // do we hand over to the next machine?
      uint32_t nextmachine;
      if (params.procset.first_bit(nextmachine)) {
        // ok pass it on to the next machine
        partial_lock_request(params.srcproc,
                             nextmachine,
                             dgraph.localvid_to_globalvid(params.localvid),
                             params.neighbor_dirty_bloom,
                             params.num_unsent_syncs,
                             params.procset,
                             params.scopetype,
                             params.priority,
                             params.src_tag);        
      }
      else {
        if (params.srcproc == rmi.procid()) {
          partial_lock_completion(params.src_tag, params.num_unsent_syncs);
        }
        else {
  #ifdef DISTRIBUTED_LOCK_DEBUG
          logstream(LOG_DEBUG) << "Replying to successful remote lock of " << dgraph.local2globalvid[curv] << std::endl;
  #endif
          rmi.remote_call(params.srcproc,
                          &distributed_mutex_lock<GraphType>::partial_lock_completion,
                          (size_t)params.src_tag,
                          params.num_unsent_syncs);
        }
      }
      partiallock_lock.lock();
      partiallock_continuation.erase(ptr);
      partiallock_lock.unlock();
    }
  
    /**
       Issue a deferred lock of type locktype 
       on lock[id] using req as the request handler.
       returns true if lock completes immediately.
       returns false otherwise.
       Calling this function requires great care as the  continuation params
       must be complete and valid at this point. 
       If this function return false, the caller must assume that the 
       continuation params may be invalid or even no longer exist.
    */
    bool issue_deferred_lock(size_t id, deferred_rwlock::request &req,
                             bool priority, scope_range::lock_type_enum locktype) {
      ASSERT_LT(id, locks.size());
      deferred_rwlock::request* released = NULL;
      size_t numreleased = 0;
      switch(locktype) {
      case scope_range::READ_LOCK:
#ifdef DISTRIBUTED_LOCK_DEBUG
        //       logstream(LOG_DEBUG) << "read lock on " << dgraph.local2globalvid[id] << std::endl;
#endif
        if (priority) {
          numreleased = locks[id].readlock_priority(&req, released);
        }
        else {
          numreleased = locks[id].readlock(&req, released);
        }

        return complete_release(released, numreleased, &req);
      case scope_range::WRITE_LOCK:
#ifdef DISTRIBUTED_LOCK_DEBUG
        //        logstream(LOG_DEBUG) << "write lock on " << dgraph.local2globalvid[id] << std::endl;
#endif
        if (priority) {
          return locks[id].writelock(&req);
        }
        else {
          return locks[id].writelock_priority(&req);
        }
      default:
        return true;
      }
    }

    void issue_deferred_unlock(size_t id,
                               scope_range::lock_type_enum locktype) {  
                           
      ASSERT_LT(id, locks.size());
      deferred_rwlock::request* released = NULL;
      size_t numreleased = 0;
      switch(locktype) {
      case scope_range::READ_LOCK:
#ifdef DISTRIBUTED_LOCK_DEBUG
        //        logstream(LOG_DEBUG) << "read unlock on " << dgraph.local2globalvid[id] << std::endl;
#endif
        numreleased = locks[id].rdunlock(released);
        complete_release(released, numreleased, NULL);
        break;
      case scope_range::WRITE_LOCK:
#ifdef DISTRIBUTED_LOCK_DEBUG
        //       logstream(LOG_DEBUG) << "write unlock on " << dgraph.local2globalvid[id] << std::endl;
#endif
        numreleased = locks[id].wrunlock(released);
        complete_release(released, numreleased, NULL);
        break;
      default:
	break;
      }
    }

    void partial_unlock(vertex_id_type globalvid, size_t scopetypeint) {
#ifdef DISTRIBUTED_LOCK_DEBUG
      //    logstream(LOG_DEBUG) << rmi.procid() << ": p-unlock request for " << globalvid << std::endl;
#endif

      scope_range::scope_range_enum scopetype = (scope_range::scope_range_enum)(scopetypeint);


      vertex_id_type localvid = dgraph.globalvid_to_localvid(globalvid);
    
      local_edge_list_type inedges = dgraph.get_local_store().in_edge_ids(localvid);
      local_edge_list_type outedges = dgraph.get_local_store().out_edge_ids(localvid);

      size_t inidx = 0;
      size_t outidx = 0;
      // quick check for fast path
      if (adjacent_vertex_lock_type(scopetype) == scope_range::NO_LOCK) {
        inidx = (vertex_id_type)(-1);
        outidx = (vertex_id_type)(-1);
      }
      bool curunlocked = false;
      vertex_id_type curv = localvid;
    
      vertex_id_type inv  =  (inedges.size() > inidx) ? 
        dgraph.get_local_store().source(inedges[inidx]) : vertex_id_type(-1);
      vertex_id_type outv  = (outedges.size() > outidx) ? 
        dgraph.get_local_store().target(outedges[outidx]) : vertex_id_type(-1);

      // iterate both in order and unlock
      // include the current vertex in the iteration
      while (inidx < inedges.size() || outidx < outedges.size()) {
        if (!curunlocked && curv < inv  && curv < outv) {
          curunlocked = true;
          if (dgraph.localvid_is_ghost(curv) == false) {
            issue_deferred_unlock(curv, central_vertex_lock_type(scopetype));
          }
        } 
        else if (inv < outv) {
          ++inidx;
          if (dgraph.localvid_is_ghost(inv) == false) {
            issue_deferred_unlock(inv, adjacent_vertex_lock_type(scopetype));
          }
          inv = (inedges.size() > inidx) ? dgraph.get_local_store().source(inedges[inidx]) : (vertex_id_type)(-1);
        } 
        else if (outv < inv) {
          ++outidx;
          if (dgraph.localvid_is_ghost(outv) == false) {
            issue_deferred_unlock(outv, adjacent_vertex_lock_type(scopetype));
          }
          outv = (outedges.size() > outidx) ? dgraph.get_local_store().target(outedges[outidx]) : (vertex_id_type)(-1);
        } 
        else if (inv == outv){
          ++inidx; ++outidx;
          if (dgraph.localvid_is_ghost(outv) == false) {
            issue_deferred_unlock(outv, adjacent_vertex_lock_type(scopetype));
          }
          inv  = (inedges.size() > inidx) ? dgraph.get_local_store().source(inedges[inidx]) : (vertex_id_type)(-1);
          outv  = (outedges.size() > outidx) ? dgraph.get_local_store().target(outedges[outidx]) : (vertex_id_type)(-1);
        }
      }
      // just in case we never got around to locking it
      if (!curunlocked && dgraph.localvid_is_ghost(curv) == false) {
        issue_deferred_unlock(curv, central_vertex_lock_type(scopetype));
      }
    }

    /**
       Starts the continuation on a collection of released requests
       beginning at the link list pointed to by "released" and for 
       "numreleased" entries.
       If watch is one of the released requests, this function
       will not call the continuation on the "watch" request, and will
       return true.
       Returns false otherwise.
    */
    bool complete_release(deferred_rwlock::request* released,
                          size_t numreleased,
                          deferred_rwlock::request* watch) {
      bool ret = false;
      for (size_t i = 0;i < numreleased; ++i) {
        deferred_rwlock::request* nextptr = (deferred_rwlock::request*)(released->next);
        if (released == watch) {
          ret = true;
        }
        else {
          // decompress the pointer
          size_t ptr = released->id;
          ptr = ptr << 2;
          typename lazy_deque<partiallock_cont_params>::value_type*
            realptr = (typename lazy_deque<partiallock_cont_params>::value_type*)(ptr);
          continue_partial_lock(realptr);
        }
        released = nextptr;
      }
      return ret;
    }
  };
#ifdef DISTRIBUTED_LOCK_DEBUG
#undef DISTRIBUTED_LOCK_DEBUG
#endif
#undef COMPILER_WRITE_BARRIER
} // namespace graphlab
#endif

