#ifndef GRAPHLAB_DISTRIBUTED_FULLSWEEP_SDM_HPP
#define GRAPHLAB_DISTRIBUTED_FULLSWEEP_SDM_HPP

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
  class distributed_fullsweep_sdm : 
    public ishared_data_manager<Graph> {
    // Typedefs
    // ==============================================================>
  public:

    typedef Graph graph_type;
    typedef distributed_fullsweep_sdm<Graph> base;

    typedef typename base::iscope_type iscope_type;
    typedef typename base::iscope_factory_type iscope_factory_type;
    typedef typename base::sync_function_type sync_function_type;
    typedef typename base::apply_function_type apply_function_type;
    typedef typename base::merge_function_type merge_function_type;


  private:
    static distributed_fullsweep_sdm<Graph>* receive_target;
    
    struct sync_params{
      sync_function_type sync_fun;
      apply_function_type apply_fun;
      merge_function_type merge_fun;
      size_t rangelow;
      size_t rangehigh;
      size_t nmyverts_in_range;
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
     
  public:

    distributed_fullsweep_sdm(distributed_control &dc, unsigned int _numcpus, Graph &_g) : 
      dc(dc), dht(dc),constant_map(11), graph(&_g), numcpus(_numcpus), fast_data_vector(256, cached_any()) { 
      receive_target = this;
    }
    


    /** Set an immutable constant */
    void set_constant(size_t index, const any& new_value) {
      setconstreplies.value = 0;
      for (procid_t i = 0; i < dc.numprocs(); ++i) {
        dc.remote_callxs(i,
                         distributed_fullsweep_sdm<Graph>::set_constant_handler,
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
    
    
    // Has to be called on all procs.
    void set_fullsweep_sync(size_t index,
                            sync_function_type sync_fun,
                            apply_function_type apply_fun,
                            merge_function_type merge_fun,
                            const any& zero,
                            size_t sync_interval = -1,
                            scope_range::scope_range_enum scopetype = scope_range::READ_CONSISTENCY,
                            size_t rangelow = 0,
                            size_t rangehigh = -1
                            ) {
      if (rangehigh < 0) rangehigh = graph->num_vertices();
      
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
      
      // Compute how many vertices in the range 
      sync.params.nmyverts_in_range = 0;
      foreach(vertex_id_t vid, graph->my_vertices()) {
        if (vid>=rangelow && vid < rangehigh)  sync.params.nmyverts_in_range++;
      }
      std::cout << "### " <<  "Vertices in range: " <<  sync.params.nmyverts_in_range << std::endl;
        
      if (dc.procid() == 0) {
        create_atomic(index, zero);
      }
      
      /* Initialize sync state */
      sync_state &sstate = sync_progression[index];
      sstate.index = index;
    
      sstate.acc = std::vector<any>(numcpus, zero);
      sstate.mergeacc = zero;
      sstate.global_merge_counter = 0;
       
      dc.barrier();
    }
    
    
    void create_atomic(size_t index, const any& initial_value) {
      dht.set(index, initial_value);
    }

    // Execute sync on local vertices
    void sync_from_local(size_t index) {
      std::cout << "### " <<  "Doing local sync " << index << std::endl;
      sync_task& stask = sync_map[index];

      foreach(vertex_id_t v, graph->my_vertices()) {
        if (v >= stask.params.rangelow && v < stask.params.rangehigh) {
          iscope_type* scope = scope_factory->get_scope(0, v);
          exec_sync_on_vertex(0, scope, index);
          scope->commit();
          scope_factory->release_scope(scope); 
        }     
      }
      std::cout << "### " <<  "Waiting for local sync " << index << std::endl;
      dc.barrier();
    }
   
   
    /**
     * \brief run sync on a particular field
     *
     * User code should not call this function and instead use the
     * sync which takes a graph argument.
     *
     * \todo this should be optimized
     */
    virtual void sync(size_t index) {}
    
    /**
     * \brief Run the particular sync task using the graph data.
     *
     * This function should not be called by user code outside of the
     * graphlab engine.
     *
     * \todo this should be optimized
     *
     */
    virtual void sync_all() {}

    
    /**
     * \brief Run the particular sync task using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine.
     *
     */
    virtual void sync(Graph& graph, size_t index) {}


    /**
     * \brief Run all sync tasks using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine.
     *
     */

    virtual void sync_all(Graph& graph) {}
    
    

    /** Notify that a sync may be needed on index */
    virtual void signal(size_t index) {}
    virtual void signal_all() {}

    
    /** syncs as soon as possible */
    virtual void trigger_sync(size_t index) {}
    virtual void trigger_sync_all() {}
   
    virtual void set_sync(size_t index,
                          sync_function_type sync_fun,
                          apply_function_type apply_fun,
                          const any& zero,
                          size_t sync_interval,
                          size_t rangelow = 0,
                          size_t rangehigh = -1) {
      assert(false);
    }
    
    // Merge accumulator removed from a remote node.
    void global_merge(size_t index, any acc) {
      sync_state& sstate = sync_progression[index];
      sync_task& stask = sync_map[index];

      sstate.mergelock.lock();
        
      stask.params.merge_fun(index, *this, sstate.mergeacc, acc);
      sstate.global_merge_counter++;
        
      std::cout << "### " <<  "Merge counter for " << index << " now " <<  sstate.global_merge_counter << std::endl;
        
      // Was this last merge?
      if ( sstate.global_merge_counter == dc.numprocs()) {
        std::cout << "### " <<  "Received last merge for " << index << " ..." << std::endl;
        atomic_apply(index, stask.params.apply_fun, sstate.mergeacc);
        sstate.global_merge_counter = 0;
        sstate.mergeacc = stask.zero;
      }
        
      sstate.mergelock.unlock();
    }
   
    void merge_and_apply(size_t index) {
      sync_task& stask = sync_map[index];
      any mergeacc = stask.zero;
      sync_state& sstate = sync_progression[index];

      // Merge over cpus
      if (stask.params.merge_fun != NULL) {
        //   std::cout << "### " <<  "Merging " << index << " over cpus. " << std::endl;
        for(size_t icpu = 0; icpu < numcpus; icpu++) {
          stask.params.merge_fun(index, *this, mergeacc, sstate.acc[icpu]);
          sstate.acc[icpu] = stask.zero;
        }
            
        // Send to proc id 0
        dc.remote_callxs(0, 
                         distributed_fullsweep_sdm<Graph>::global_merge_handler,
                         NULL,
                         0,
                         index,
                         mergeacc);
                            
        // Barrier
        ASSERT_EQ(sstate.counter.dec((size_t) stask.params.nmyverts_in_range), 0);
        // std::cout << "### " <<  "Sync iteration ok : " << index << std::endl;
      }
    }
   
    void fillcache(size_t index, any value) {
      fast_data_vector[index].valid = true;
      fast_data_vector[index].payload = value;
    }


   
    void exec_sync_on_vertex(size_t cpuid, iscope_type * scope, size_t index) {
      sync_state& sstate = sync_progression[index];
      sync_task& stask = sync_map[index];
      
        
      unsigned int nvert = stask.params.nmyverts_in_range;
      timer t;
      t.start();
      bool waslast = false;
      sstate.lock.lock();
   
      size_t c = sstate.counter.inc();
      if (c == nvert) {
        waslast = true;
        sstate.wait_for_update = true;
      } else {
        // If was last, then we release lock much later!
        sstate.lock.unlock();
      }	
      

      sstate.justsyncing.inc();
      stask.params.sync_fun(index, *this, *scope, sstate.acc[cpuid]);
      sstate.justsyncing.dec();

      if (waslast) {            
         while(sstate.justsyncing.value > 0) {
          usleep(50);
        }
            

        merge_and_apply(index);
        ASSERT_EQ(sstate.counter.value, 0);
            
        while(sstate.wait_for_update) {
          sched_yield();
        }
        sstate.lock.unlock();
            
        // Fill the 1st level cache with the new data
        fillcache(index, atomic_get(index));

        printf("Last sync on %d took %lf secs\n", (int) index, t.current_time());
      }
      ASSERT_FALSE(c > nvert);
        
        
    }
   
    // this must be called ocassionally
    void progress(size_t cpuid, iscope_type * scope) {
      typename sync_map_type::iterator keyiter;
      vertex_id_t v = scope->vertex();
      for(keyiter = sync_map.begin(); keyiter != sync_map.end(); keyiter++) {
        sync_task stask = keyiter->second;
        if (v >= stask.params.rangelow && v < stask.params.rangehigh) {
          exec_sync_on_vertex(cpuid, scope, keyiter->first);
        }
      }
        
    }

    any atomic_get(size_t index) const {
      any value;
      ASSERT_MSG(dht.get(index, value), "DHT Failed to get entry %d", index);
      return value;
    }


    any get(size_t index) const {
      
      if (fast_data_vector[index].valid) {
      	return fast_data_vector[index].payload; 
      }
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
                     distributed_fullsweep_sdm<Graph>::distributed_fullsweep_sdm_apply,
                     (void*)(strm.str().c_str()),
                     strm.str().length(),
                     size_t(index),
                     size_t(fun),
                     (handlerarg_t)(&trigger));
      while(trigger == 0) sched_yield();
      dht.invalidate(index);
      fast_data_vector[index].valid = false;
      // Send ready message to everyone
      for(size_t proc=0; proc<dc.numprocs(); proc++) {
      	dc.remote_call(proc, distributed_fullsweep_sdm<Graph>::global_sync_ready_handler,
                       NULL, 0, size_t(index));
      }
    }


  public:
    //=========================== MESSAGE HANDLERS =============================
  
    static void set_constant_handler(distributed_control& dc, 
                                     procid_t source,  
                                     void* ptr,    //unused
                                     size_t len,
                                     handlerarg_t key,
                                     any &data) {
      distributed_fullsweep_sdm<Graph> &dsdm = *(distributed_fullsweep_sdm<Graph>::receive_target);
      dsdm.constant_map.insert(key, data);
      // issue reply
      dc.remote_call(source,
                     distributed_fullsweep_sdm<Graph>::set_constant_reply_handler,
                     NULL, 0);
    }
    
    static void set_constant_reply_handler(distributed_control& dc, 
                                           procid_t source,  
                                           void* ptr,    //unused
                                           size_t len) {
      distributed_fullsweep_sdm<Graph> &dsdm = *(distributed_fullsweep_sdm<Graph>::receive_target);
      dsdm.setconstreplies.inc();
    }
    
    static void distributed_fullsweep_sdm_apply(distributed_control& dc, 
                                                procid_t source,  
                                                void* ptr,    //serialized any
                                                size_t len,
                                                handlerarg_t key,
                                                handlerarg_t applyfn,
                                                handlerarg_t applyptr) {
      // ptr/len has the stream for the new data
      boost::iostreams::stream<boost::iostreams::array_source>  streamin((char*)(ptr), len);
      iarchive iarc(streamin);
      
      distributed_fullsweep_sdm<Graph> &dsdm = *(distributed_fullsweep_sdm<Graph>::receive_target);

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
        dc.remote_call(source, set_ptr_to_value_1_handler, NULL, 0,applyptr);
      }
    }

    static void global_sync_ready_handler(distributed_control &dc,
                                          procid_t source,
                                          void* ptr,
                                          size_t len,
                                          size_t index) {
      distributed_fullsweep_sdm<Graph> &dsdm = *(distributed_fullsweep_sdm<Graph>::receive_target);
      
      
  //    std::cout << "### " <<  "global_sync_ready:" << dc.procid() << ": index: " << std::endl;
      
      // Invalidate from cache
      if (dsdm.dht.key_node_hash(index) != dc.procid()) {
      	dsdm.dht.invalidate(index);
      	dsdm.fast_data_vector[index].valid = false;	
      }
      dsdm.sync_progression[index].wait_for_update = false;
      //     std::cout << "### " <<  "Ready :: " << index << std::endl;
    }

    static void global_merge_handler(distributed_control &dc,
                                     procid_t source,
                                     void* ptr,
                                     size_t len,
                                     size_t index,
                                     any &acc) {
      ASSERT_EQ(dc.procid(), 0);
      logstream(LOG_INFO) << "Global Merge on " << index << " received from " << source << std::endl;
      distributed_fullsweep_sdm<Graph> &dsdm = *(distributed_fullsweep_sdm<Graph>::receive_target);
      ASSERT_TRUE(dsdm.sync_progression.find(index) != dsdm.sync_progression.end());
      dsdm.global_merge(index, acc);
    }

   
    // ======= Stuff called only by the engine =====================

    // must be called by engine
    void set_scope_factory(iscope_factory_type* factory) {
      scope_factory = factory;
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
      size_t index;
      mutex lock;
      mutex mergelock;
      std::vector< any > acc;
      bool wait_for_update;
      atomic< size_t > counter;
      int global_merge_counter;
      any mergeacc; // only used if merge_fun is set
      atomic< int > justsyncing;
    };    

    std::map<size_t, sync_state> sync_progression;

    
    // constants are stored here
    constant_map_type constant_map;
    atomic<size_t> setconstreplies;
    
    graph_type* graph;
    iscope_factory_type* scope_factory;
    // sync progression 
    
    unsigned int  numcpus;

    
    struct cached_any {
      bool valid;
      any payload;
      cached_any() {valid=false;}
    };
    
    // Fast cache for shared data
    std::vector< cached_any > fast_data_vector;
    
  }; // end of class thread_share_data

  template <typename Graph>
  distributed_fullsweep_sdm<Graph>* distributed_fullsweep_sdm<Graph>::receive_target;

}; // end graphlab namespace
#include <graphlab/macros_undef.hpp>



#endif
