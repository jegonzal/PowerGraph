#ifndef GRAPHLAB_DISTRIBUTED_SHARED_DATA_HPP
#define GRAPHLAB_DISTRIBUTED_SHARED_DATA_HPP

#include <unistd.h>
#include <map>
#include <limits>
#include <boost/iostreams/stream.hpp>

#include <graphlab/util/timer.hpp>
#include <graphlab/util/stl_util.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/shared_data/ishared_data.hpp>
#include <graphlab/shared_data/ishared_data_manager.hpp>
#include <graphlab/distributed/distributed_scope.hpp>
#include <graphlab/distributed/distributed_hash_table.hpp>
#include <graphlab/distributed/graph_lock_manager.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/util/synchronized_unordered_map.hpp>


#include <graphlab/macros_def.hpp>

namespace graphlab {


  
  template<typename Graph>
  class distributed_shared_data : 
    public ishared_data_manager<Graph> {
    // Typedefs
    // ==============================================================>
  public:

    typedef Graph graph_type;
    typedef distributed_shared_data<Graph> base;

    typedef typename base::iscope_type iscope_type;
    typedef typename base::iscope_factory_type iscope_factory_type;
    typedef typename base::sync_function_type sync_function_type;
    typedef typename base::apply_function_type apply_function_type;
    typedef typename base::merge_function_type merge_function_type;


  private:
    static distributed_shared_data<Graph>* receive_target;
    
    struct sync_params{
      sync_function_type sync_fun;
      apply_function_type apply_fun;
      merge_function_type merge_fun;
      size_t rangelow;
      size_t rangehigh;
      scope_range::scope_range_enum scopetype;
      void save(graphlab::oarchive &oarc) const{
        serialize(oarc, this, sizeof(sync_params));  
      }
      void load(graphlab::iarchive &iarc) {
        deserialize(iarc, this, sizeof(sync_params));
      }    
    };
    struct sync_task {
      sync_params params;
      size_t sync_interval;
      float next_time;
      bool sync_in_progress;
      any zero;
      mutex lock; // not really used at all
      sync_task() :
        sync_interval(-1),
        next_time(0) { }
    };
   
    typedef synchronized_unordered_map<any> constant_map_type;
    typedef std::map<size_t, sync_task> sync_map_type;
    
    atomic<size_t> pendinglocks;
    atomic<size_t> syncs_in_progress;
  public:

    distributed_shared_data(distributed_control &dc) : 
                    dc(dc), dht(dc),constant_map(11), glm(NULL), graph(NULL) { 
      receive_target = this;
      pendinglocks.value = 0;
      syncs_in_progress.value = 0;
      dht.set_pushed_updates(true);
    }
    


    /** Set an immutable constant */
    void set_constant(size_t index, const any& new_value) {
      setconstreplies.value = 0;
      for (procid_t i = 0; i < dc.numprocs(); ++i) {
        dc.remote_callxs(i,
                       distributed_shared_data<Graph>::set_constant_handler,
                       NULL,
                       0,
                       index,
                       new_value);
      }
      while (setconstreplies.value != dc.numprocs()) {
        sched_yield();
      }
    } // end of 

    //! \todo RETURN BY REFERENCE
    const any& get_constant(size_t index) const {
      typename constant_map_type::const_datapointer iter = constant_map.find(index);
      ASSERT_TRUE(iter.first);
      return *(iter.second);
    } // end of get constant
    
    /** register a sync.  Requies that the data be already set. */
    void set_sync(size_t index,
                  sync_function_type sync_fun,
                  apply_function_type apply_fun,
                  const any& zero,
                  size_t sync_interval = -1,
                  size_t rangelow = 0,
                  size_t rangehigh = -1) {
      sync_task& sync = sync_map[index];
      sync.params.sync_fun = sync_fun;
      sync.params.apply_fun = apply_fun;
      sync.params.merge_fun = NULL;
      sync.params.rangelow = rangelow;
      sync.params.rangehigh = rangehigh;
      sync.params.scopetype = scope_range::READ_CONSISTENCY;
      sync.zero = zero;
      sync.sync_interval = sync_interval;
      sync.next_time = sync_interval;
      sync.sync_in_progress = false;
      ASSERT_LE(sync.params.rangelow, sync.params.rangehigh);
      create_atomic(index, zero);
      
      // create the sync state on all machines
      std::vector<size_t> syncids = keys_as_vector(sync_map);
      if (dc.procid() == 0) {
        for (procid_t i = 0; i < dc.numprocs(); ++i) {
          volatile size_t trigger = 0;
          dc.remote_call(i, 
                         distributed_shared_data<Graph>::create_lock_blocks_handler,
                         &(index),
                         sizeof(size_t),
                         (handlerarg_t)(&trigger)
                         );
          while(trigger == 0) sched_yield();
        }
      }
    }

    void set_tree_sync(size_t index,
                  sync_function_type sync_fun,
                  apply_function_type apply_fun,
                  merge_function_type merge_fun,
                  const any& zero,
                  size_t sync_interval = -1,
                  scope_range::scope_range_enum scopetype = scope_range::READ_CONSISTENCY,
                  size_t rangelow = 0,
                  size_t rangehigh = -1) {
      sync_task& sync = sync_map[index];
      sync.params.sync_fun = sync_fun;
      sync.params.apply_fun = apply_fun;
      sync.params.merge_fun = merge_fun;
      sync.params.rangelow = rangelow;
      sync.params.rangehigh = rangehigh;
      sync.params.scopetype = scopetype;
      sync.zero = zero;
      sync.sync_interval = sync_interval;
      sync.sync_in_progress = false;
      create_atomic(index, zero);

      // create the sync state on all machines
      std::vector<size_t> syncids = keys_as_vector(sync_map);
      if (dc.procid() == 0) {
        for (procid_t i = 0; i < dc.numprocs(); ++i) {
          volatile size_t trigger = 0;
          dc.remote_call(i,
                         distributed_shared_data<Graph>::create_lock_blocks_handler,
                         &(index),
                         sizeof(size_t),
                         (handlerarg_t)(&trigger)
                         );
          while(trigger == 0) sched_yield();
        }
      }
    }
    void create_atomic(size_t index, const any& initial_value) {
      dht.set(index, initial_value);
    }



    /**
     * Notify a change in the index sync dependencies.  If it is time
     * the calling thread is used to complete a sync.
     */
    void signal(size_t index) {
      typedef typename sync_map_type::iterator iterator_type;
      iterator_type iter = sync_map.find(index);
      if(iter != sync_map.end()) {
        sync_task& stask = iter->second;
        if(stask.next_time < lowres_time_millis()) {
          sync(index);
        }
      }
    }
    
    void signal_all() {
      if (signallock.try_lock()) {
        typedef typename sync_map_type::value_type value_type;
        foreach(value_type& pair, sync_map) {
          size_t index = pair.first;
          sync_task& stask = pair.second;
          if(stask.next_time < lowres_time_millis()) {
            sync(index);
          }
        }
        signallock.unlock();
      }
    }
    
    void trigger_sync(size_t index) {
      if (dc.procid() == 0) {
        typedef typename sync_map_type::iterator iterator_type;
        iterator_type iter = sync_map.find(index);
        if(iter != sync_map.end()) {
          sync_task& stask = iter->second;
          stask.next_time = lowres_time_millis();
        }
      }
      else {
        dc.remote_call(0,
                       trigger_sync_handler,
                       NULL,
                       0,
                       (handlerarg_t)index);
      }
    }
      
    void trigger_sync_all() {
      if (dc.procid() == 0) {
        typedef typename sync_map_type::value_type value_type;
        foreach(const value_type& pair, sync_map) trigger_sync(pair.first);
      }
      else {
        dc.remote_call(0,
                       trigger_sync_all_handler,
                       NULL,
                       0);
      }
    }


    /**
     * \brief run sync on everything (should be called by graphlab
     * engine)
     *
     * User sync requests should use the sync which takes a graph
     * argument.
     *
     * \todo this should be optimized
     */
    void sync_all() {
      typedef typename sync_map_type::value_type value_type;
      foreach(const value_type& pair, sync_map) sync(pair.first);
    }
    
    /**
     * \brief run sync on a particular field
     *
     * User sync requests should use the sync which takes a graph
     * argument. This can only be called by node 0. 
     * All nodes must call "progress" ocassionally.
     * The linear sync follows the procid order. proc 0 will complete the
     * local accumulation first. This will then be passed onto proc 1, and
     * so on. The final proc will perform the apply and write
     */
    void sync(size_t index) {
      ASSERT_EQ(dc.procid(), 0);
      ASSERT_TRUE(sync_map.find(index) != sync_map.end());
      sync_map_lock.lock();
      if (sync_map[index].sync_in_progress) {
        sync_map_lock.unlock();
        return;
      }
      sync_map[index].sync_in_progress = true;
      syncs_in_progress.inc();
      sync_map_lock.unlock();
      // issue the sync
      // if we do not support merge, do the linear sync
      if (sync_map[index].params.merge_fun == NULL) {
        issue_local_sync(index,
                        sync_map[index].params.sync_fun,
                        sync_map[index].params.apply_fun,
                        sync_map[index].params.merge_fun,
                        sync_map[index].zero,
                        sync_map[index].params.rangelow,
                        sync_map[index].params.rangehigh,
                        sync_map[index].params.scopetype);
      }
      else {
        // other wise broad cast sync to everyone
        for (size_t i = 0;i < dc.numprocs(); ++i) {
          dc.remote_callxs(i,
                  distributed_shared_data<Graph>::begin_local_sync_handler,
                  NULL,
                  0,
                  index,
                  sync_map[index].params,
                  sync_map[index].zero);
        }
      }
    } // end of sync

    /**
     * \brief Run all sync tasks using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine. 
     * TODO: This is not-parallel. (can we make it parallel? seems like will
     * be very complicated
     */
    void sync_all(Graph& graph) {
      typedef typename sync_map_type::value_type value_type;
      foreach(const value_type& pair, sync_map) sync(graph, pair.first);      
    }


    /**
     * \brief Run the particular sync task using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine. 
     * TODO: This is not-parallel. (can we make it parallel? seems like will
     * be very complicated
     */
    void sync(Graph& graph, size_t index) {
      // issue a sequence of lock requests for all vertices in this partition
      //foreach(v, graph.my_vertices()) {
      //}
      // Requires a scope factory
      assert(scope_factory != NULL);
      // Get the table entry
      typedef typename sync_map_type::iterator iterator_type;
      iterator_type iter = sync_map.find(index);
      assert(iter != sync_map.end());
      sync_task& sync = iter->second;
      if (sync.params.sync_fun == NULL && sync.params.apply_fun == NULL) return;
      // Ensure that the sync is well defined
/*      assert(sync.sync_fun != NULL);
      assert(sync.apply_fun != NULL);*/
      any accumulator = sync.zero;
      if(sync.lock.try_lock()) {
        // Copy the accumulator
        if (sync.params.sync_fun != NULL) {
          // Try and grab the lock for the sync
          size_t cpuid = thread::thread_id();
          size_t h = std::min(scope_factory->num_vertices(), sync.params.rangehigh);
          for(size_t v = sync.params.rangelow; v < h; ++v) {
            // get the scope for the vertex
            iscope_type* scope = scope_factory->get_scope(cpuid, v);
            assert(scope != NULL);
            // Apply the sync function
            sync.params.sync_fun(index, *this, *scope, accumulator);
            // Commit and free the scope
            scope->commit();
            scope_factory->release_scope(scope);
          }
        }
        if (sync.params.apply_fun != NULL) {
          // Apply the result
          atomic_apply(index, sync.params.apply_fun, accumulator);        
        }
        // release the lock
        sync.lock.unlock();
        // set the next sync time
        size_t lrtime = lowres_time_millis();
        sync.next_time = lrtime + sync.sync_interval;
        if (sync.next_time < lrtime) { // overflow 
          sync.next_time = -1;
        }
      }      
    } // end of sync
    
  void sync_from_local(size_t index) {
      if (dc.procid() == 0) {
        ASSERT_TRUE(sync_map.find(index) != sync_map.end());
        sync_map_lock.lock();
        if (sync_map[index].sync_in_progress) {
          sync_map_lock.unlock();
          return;
        }
        sync_map[index].sync_in_progress = true;
        syncs_in_progress.inc();
        sync_map_lock.unlock();
        // issue the sync
        // if we do not support merge, do the linear sync
        if (sync_map[index].params.merge_fun == NULL) {
          issue_local_sync(index,
                          sync_map[index].params.sync_fun,
                          sync_map[index].params.apply_fun,
                          sync_map[index].params.merge_fun,
                          sync_map[index].zero,
                          sync_map[index].params.rangelow,
                          sync_map[index].params.rangehigh,
                          sync_map[index].params.scopetype);
        }
        else {
          // other wise broad cast sync to everyone
          for (size_t i = 0;i < dc.numprocs(); ++i) {
            dc.remote_callxs(i,
                    distributed_shared_data<Graph>::begin_local_sync_handler,
                    NULL,
                    0,
                    index,
                    sync_map[index].params,
                    sync_map[index].zero);
          }
        }
      }
      while(!progress(0)) {
        sched_yield();
      }
      while(progress(0)) {
        sched_yield();
      }
      while(syncs_in_progress.value > 0) {
        progress(0);
        sched_yield();
      }
    } // end of sync

    bool has_pending_tasks() {
      return (pendinglocks.value > 0) || (syncs_in_progress.value > 0);
    }
    // this must be called ocassionally
    bool progress(size_t cpuid) {
      // pop an element of the sync progression and run it
      std::pair<size_t, bool> elem = active_sync_progressions.try_dequeue();
      if (elem.second == false) {
        return false; // nothing to do
      }
      // ok, we have a sync to progress to run
      // grab a reference to the state we are progressing
      sync_state& curstate = sync_progression[elem.first];
      progress_state(cpuid, curstate);
      if (curstate.active || curstate.mergewaiting) active_sync_progressions.enqueue(elem.first);
      return true;
    }

    any atomic_get(size_t index) const {
      any value;
      ASSERT_MSG(dht.get(index, value), "DHT Failed to get entry %d", index);
      return value;
    }

    any get(size_t index) const {
      any value;
      ASSERT_MSG(dht.get_cached(index, value), "DHT Failed to get entry %d", index);
      return value;
    }

    void atomic_set(size_t index, const any& data) {
      dht.set(index, data);
    }


    any atomic_exchange(size_t index, const any& data) {
      any oldvalue;
      dht.exchange(index, data, oldvalue);
      return oldvalue;
    }


    void atomic_apply(size_t index,
                      apply_function_type fun,
                      const any& data) {
      std::stringstream strm;
      oarchive arc(strm);
      arc << data;
      volatile size_t trigger = 0;
      dc.remote_call(dht.key_node_hash(index), 
                     distributed_shared_data<Graph>::distributed_shared_data_apply,
                     (void*)(strm.str().c_str()),
                     strm.str().length(),
                     size_t(index),
                     size_t(fun),
                     (handlerarg_t)(&trigger));
      while(trigger == 0) sched_yield();
      
      dht.invalidate(index);
    }


  public:
    //=========================== MESSAGE HANDLERS =============================
  
    static void set_constant_handler(distributed_control& dc, 
                             procid_t source,  
                             void* ptr,    //unused
                             size_t len,
                             handlerarg_t key,
                             any &data) {
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      dsdm.constant_map.insert(key, data);
      // issue reply
      dc.remote_call(source,
                     distributed_shared_data<Graph>::set_constant_reply_handler,
                     NULL, 0);
    }
    
    static void set_constant_reply_handler(distributed_control& dc, 
                                    procid_t source,  
                                    void* ptr,    //unused
                                    size_t len) {
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      dsdm.setconstreplies.inc();
    }
    
    static void trigger_sync_handler(distributed_control& dc, 
                                    procid_t source,  
                                    void* ptr,    //unused
                                    size_t len,
                                    handlerarg_t index) {
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      dsdm.trigger_sync(index);
    }

    static void trigger_sync_all_handler(distributed_control& dc, 
                                    procid_t source,  
                                    void* ptr,    //unused
                                    size_t len) {
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      dsdm.trigger_sync_all();
    }
    
    static void create_lock_blocks_handler(distributed_control& dc, 
                             procid_t source,  
                             void* ptr,    // sync ids
                             size_t len,
                             handlerarg_t replyptr) { // numsyncids * size_t
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      size_t numsyncids = len / sizeof(size_t);
      size_t* syncids = (size_t*)ptr;
      
      std::vector<dist_scope_request> emptyreqs;
      for (size_t i = 0;i < numsyncids; ++i) {
//        logstream(LOG_INFO) << "Creating Sync Progression state for " << syncids[i] << std::endl;
        sync_state &sstate = dsdm.sync_progression[syncids[i]];
        sstate.index = syncids[i];
        sstate.lockblockid = dsdm.glm->block_deferred_lock(emptyreqs, 300);  // throttle_rate
        sstate.params.sync_fun = NULL;
        sstate.params.apply_fun = NULL;
        sstate.params.merge_fun = NULL;
        sstate.params.rangehigh = -1;
        sstate.params.rangelow = 0;
        sstate.params.scopetype = scope_range::READ_CONSISTENCY;
        sstate.active = false;
        sstate.mergewaiting = false;
        sstate.acc = any();
        sstate.mergeacc = any();
        sstate.mergectr = 0;
      }
      dc.remote_call(source, set_ptr_to_value_1_handler, NULL, 0,replyptr);
    }
/*
This handler can only be called on the machine holding on to the DHT
entry of the key
*/
    static void distributed_shared_data_apply(distributed_control& dc, 
                                              procid_t source,  
                                              void* ptr,    //serialized any
                                              size_t len,
                                              handlerarg_t key,
                                              handlerarg_t applyfn,
                                              handlerarg_t applyptr) {
      // ptr/len has the stream for the new data
      boost::iostreams::stream<boost::iostreams::array_source>  streamin((char*)(ptr), len);
      iarchive iarc(streamin);
      
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);

      distributed_hash_table::map_type::datapointer iter = dsdm.dht.data.find(key);
      if (iter.first == false) {
        logstream(LOG_ERROR) <<  "Attempted apply to empty DSDM key " << key << std::endl;
      } else {
        any newdata;
        iarc >> newdata;
        iter.second->first.writelock();
        apply_function_type fun = apply_function_type(applyfn);
        fun(key, dsdm, iter.second->second, newdata); 
        iter.second->first.unlock();
        dsdm.dht.modified(key);
        dc.remote_call(source, set_ptr_to_value_1_handler, NULL, 0,applyptr);
      }
    }

    static void begin_local_sync_handler(distributed_control &dc,
                                         procid_t source,
                                         void* ptr,
                                         size_t len,
                                         size_t index, 
                                         sync_params params,
                                         any& acc) {
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      dsdm.issue_local_sync(index, 
                            params.sync_fun, 
                            params.apply_fun,
                            params.merge_fun,
                            acc,
                            params.rangelow,
                            params.rangehigh,
                            params.scopetype);
    }
    
    // only callable at node 0
    static void sync_complete_handler(distributed_control &dc,
                                         procid_t source,
                                         void* ptr,
                                         size_t len,
                                         size_t index) {
//      logstream(LOG_INFO) << "Sync " << index << " complete "<< std::endl;
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      dsdm.sync_map_lock.lock();
      ASSERT_TRUE(dsdm.sync_map.find(index) != dsdm.sync_map.end());
      sync_task& stask = dsdm.sync_map[index];
      stask.sync_in_progress = false;
      dsdm.syncs_in_progress.dec();
      // update the next time
      
      // set the next sync time
      size_t lrtime = lowres_time_millis();
      stask.next_time = lrtime + stask.sync_interval;
      if (stask.next_time < lrtime) { // overflow 
        stask.next_time = -1;
      }
      dsdm.sync_map_lock.unlock();
      std::cout << "Completion of " << index << ". Syncs in progress: " << dsdm.syncs_in_progress.value << std::endl;
    }

    static void tree_merge_handler(distributed_control &dc,
                                         procid_t source,
                                         void* ptr,
                                         size_t len,
                                         size_t index,
                                         any &acc) {
      logstream(LOG_INFO) << "Tree Merge on " << index << " receieved"<< std::endl;
      distributed_shared_data<Graph> &dsdm = *(distributed_shared_data<Graph>::receive_target);
      ASSERT_TRUE(dsdm.sync_progression.find(index) != dsdm.sync_progression.end());
      dsdm.perform_tree_merge(dsdm.sync_progression[index], acc, true);
    }
// ======= Stuff called only by the engine =====================

    // must be called by engine
    void set_scope_factory(iscope_factory_type* factory) {
      scope_factory = factory;
    }

    // this must be called simultenously by all processors
    // this will call all other MPI nodes to create the lock blocks they need
    void set_lock_manager(graph_type* _graph, graph_lock_manager<Graph>* _glm) {
      graph = _graph;
      glm = _glm;
    }

    // Data Members
    // ==============================================================>
  private:
    // for communication
    distributed_control &dc;
    
    // DHT stores the atomics
    mutable distributed_hash_table dht;
    
    // this manages the sync instructions. This is only set on the root node
    sync_map_type sync_map;
    mutex sync_map_lock;
    
    struct sync_state {
      sync_params params;
      size_t index;
      size_t lockblockid;
      any acc;
      any mergeacc; // only used if merge_fun is set
      size_t mergectr;
      mutex lock;
      bool active;
      bool mergewaiting;
    };
    blocking_queue<size_t> active_sync_progressions;
    std::map<size_t, sync_state> sync_progression;
    
    // constants are stored here
    constant_map_type constant_map;
    atomic<size_t> setconstreplies;
    
    graph_lock_manager<Graph>* glm;
    graph_type* graph;
    iscope_factory_type* scope_factory;
    // sync progression 
    
    mutex signallock;
    
    /// Issue the sync to the local set of vertices
    void issue_local_sync(size_t index, 
                          sync_function_type sync_fun,
                          apply_function_type apply_fun,
                          merge_function_type merge_fun,
                          any& acc,
                          size_t rangelow,
                          size_t rangehigh,
                          scope_range::scope_range_enum scopetype) {
//      logstream(LOG_INFO) << "Issuing local sync on : " << index << std::endl;
      ASSERT_TRUE(sync_progression.find(index) != sync_progression.end());
      // acquire a lock on all of my vertices
      std::vector<dist_scope_request> reqs;
      size_t numissued = 0;
      if (rangelow == 0 && rangehigh >= graph->num_vertices()) {
        // if it is the complete range. do the resize once
        reqs.resize(graph->my_vertices().size());
        for (size_t i = 0;i < graph->my_vertices().size(); ++i) {
          reqs[i].vertex = graph->my_vertices()[i];
          reqs[i].scoperange = scopetype; //TODO: what shuold this be??
          ++numissued;
        }
      }
      else {
         for (size_t i = 0;i < graph->my_vertices().size(); ++i) {
          vertex_id_t v = graph->my_vertices()[i];
          if (v >= rangelow && v < rangehigh) {
            reqs.push_back(dist_scope_request(v, scopetype));
            ++numissued;
          }
        }
      }
      //std::random_shuffle(reqs.begin(), reqs.end());
      pendinglocks.inc(numissued);
      // update the sync state
      sync_state &sstate = sync_progression[index];
      sstate.acc = acc;
      sstate.mergeacc = any();
      sstate.params.sync_fun = sync_fun;
      sstate.params.apply_fun = apply_fun;
      sstate.params.merge_fun = merge_fun;
      sstate.params.rangelow = rangelow;
      sstate.params.rangehigh = rangehigh;
      sstate.params.scopetype = scopetype;
      sstate.active = true;
      sstate.mergewaiting = true;
      logstream(LOG_INFO) << "local sync issued on : " << index 
                          << " with " << numissued << " vertices" << " " << scopetype << std::endl;
      glm->block_add_deferred_lock(sstate.lockblockid, reqs);

      active_sync_progressions.enqueue(index);
      // issue the lock request
    }
    
    void progress_state(size_t cpuid, sync_state &state) {
      // while there are locks we can grab
      dist_scope_request scopereqdone;
      size_t remaininglocks;
      bool firstcompletion = state.active;
      if (state.active) {
        remaininglocks = glm->block_status(state.lockblockid, scopereqdone);
        while (scopereqdone.vertex != vertex_id_t(-1)) {
  //        logstream(LOG_INFO) << "Progress on: " << state.index << std::endl;
          // get the scope
          iscope_type* scope = scope_factory->get_scope(cpuid, scopereqdone.vertex);
          // call the function
          state.params.sync_fun(state.index, *this, *scope, state.acc);
          // release the scope
          // as again don't use commit. we are managing it from the blocks
          glm->block_release_partial(state.lockblockid, scopereqdone);
          scope_factory->release_scope(scope);
          remaininglocks = glm->block_status(state.lockblockid, scopereqdone);
          pendinglocks.dec();
        }
        //logstream(LOG_INFO) << "local sync progress on : " << state.index 
  //                << " with " << remaininglocks << " vertices remaining" << std::endl;
        state.active = remaininglocks > 0;
      }
      // ok. If I am totally done
      if (state.active == false && state.mergewaiting) {
        // do we do linear sync?
        if (state.params.merge_fun == NULL) {
          if (dc.procid() == dc.numprocs() - 1) {
            // if I am the last processor, I need to apply
            atomic_apply(state.index, state.params.apply_fun, state.acc);
            // call back proc 0 to tell it I am done
            dc.remote_call(0,
                          distributed_shared_data<Graph>::sync_complete_handler,
                          NULL,
                          0,
                          state.index);
          }
          else {  // if I am not the last processor, I need to push this state
          //TODO: push the state to the next processor
          dc.remote_callxs(dc.procid() + 1,
                            distributed_shared_data<Graph>::begin_local_sync_handler,
                            NULL,
                            0,
                            state.index,
                            state.params,
                            state.acc);
          }
        }
        else {
          //we can do the tree merge.
          // if mergewaiting is false then this is a fresh call to tree merge
          perform_tree_merge(state, state.acc, firstcompletion);
        }
      }
    }

    void perform_tree_merge(sync_state &state, any &acc, bool mergeacc) {
      state.lock.lock();
      if (state.mergewaiting) {
        // merge in if this is the first time
        if (mergeacc) {
          std::cout << dc.procid() << "Merging in index " << state.index << std::endl;
          if (state.mergeacc.empty()) state.mergeacc.swap(acc);
          else state.params.merge_fun(state.index, *this, state.mergeacc, acc);
          state.mergectr++;
        }

        // see if I got everyone
        // this uses tree reduction
        size_t numchildren = 0;
        numchildren = (dc.procid() * 2 + 1 < dc.numprocs()) +
                      (dc.procid() * 2 + 2 < dc.numprocs());
        size_t parent = (dc.procid() - 1) / 2;
        
        // alternative rule. Let node 0 do everything .This appears to 
        // be quite a bit faster
        parent = 0;
        if (dc.procid() == 0) numchildren = dc.numprocs() - 1;
        else numchildren = 0;
        if (state.mergectr == numchildren + 1) {
          state.mergewaiting = false;
          if (dc.procid() != 0) {
            // clear the state before I issue remote
            // lest we get caught in the next sync
            // though it shouldn't happen cos of the lock anyway
            any transmitacc;
            transmitacc.swap(state.mergeacc);
            state.mergectr = 0;
            state.params.merge_fun = NULL;
            dc.remote_callxs(parent,
                            distributed_shared_data<Graph>::tree_merge_handler,
                            NULL,
                            0,
                            state.index,
                            transmitacc);
          }
          else {
            // I am the root
            // I finish here
            atomic_apply(state.index, state.params.apply_fun, state.mergeacc);
            state.mergeacc = any();
            state.mergectr = 0;
            dc.remote_call(0,
                            distributed_shared_data<Graph>::sync_complete_handler,
                            NULL,
                            0,
                            state.index);
          }
        }
      }
      state.lock.unlock();
    }
    
  }; // end of class thread_share_data

template <typename Graph>
distributed_shared_data<Graph>* distributed_shared_data<Graph>::receive_target;

}; // end graphlab namespace
#include <graphlab/macros_undef.hpp>



#endif
