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


#ifndef GRAPHLAB_DISTRIBUTED_AGGREGATOR
#define GRAPHLAB_DISTRIBUTED_AGGREGATOR

#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/pthread_tools.hpp>

#include <graphlab/rpc/dc_dist_object.hpp>



#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename Engine>
  class distributed_aggregator {
  public:
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
    typedef typename graph_type::lvid_type lvid_type;
    typedef typename engine_type::update_functor_type update_functor_type;
    typedef typename engine_type::context_type context_type;
    typedef typename graph_type::local_edge_list_type local_edge_list_type;
    typedef typename graph_type::edge_type edge_type;

    typedef dc_dist_object< distributed_aggregator > rmi_type;


    // Global Aggregates ------------------------------------------------------
    /**
     * The base class for registered aggregation operation
     */
    struct isync {
      mutex lock;
      graphlab::barrier* barrier_ptr;
      atomic<size_t> counter;
      const size_t interval;
      isync(size_t interval) : barrier_ptr(NULL), counter(0), interval(interval) { }
      virtual ~isync() {
        if(barrier_ptr != NULL) { delete barrier_ptr; barrier_ptr = NULL;}
      }
      virtual void run_aggregator(std::string key,
                                  size_t cpuid,
                                  distributed_aggregator* dist_aggregator,
                                  rmi_type* rmi_ptr,
                                  engine_type* engine_ptr,                          
                                  graph_type*  graph_ptr) = 0;
    }; // end of isync

    /**
     * The derived class used to store the various aggregation
     * operations
     */
    template<typename Aggregator >
    struct sync : public isync {
      typedef Aggregator       aggregator_type;
      using isync::lock;
      using isync::barrier_ptr;
      using isync::counter;
      using isync::interval;
      const aggregator_type zero;
      aggregator_type shared_aggregator;
      sync(const size_t interval, const aggregator_type& zero) : 
        isync(interval), zero(zero), shared_aggregator(zero) { }
      void run_aggregator(std::string key,
                          size_t cpuid,
                          distributed_aggregator* dist_aggregator,
                          rmi_type* rmi_ptr, 
                          engine_type* engine_ptr,
                          graph_type*  graph_ptr) { 
        // Thread zero must initialize the the final shared accumulator
        const size_t nverts = graph_ptr->num_local_vertices();
        
        // construct the local (to this thread) accumulator and
        // context
        aggregator_type local_accum(zero);
        context_type context(engine_ptr, graph_ptr);

        // Apply the update to the local vertices
        for(lvid_type lvid = counter++; lvid < nverts; lvid = counter++) {
          if(local_accum.is_factorizable()) {       
            const local_edge_list_type in_edges = graph_ptr->l_in_edges(lvid);
            const local_edge_list_type out_edges = graph_ptr->l_out_edges(lvid);
            context.init(graph_ptr->global_vid(lvid), EDGE_CONSISTENCY);
            if(local_accum.gather_edges() == IN_EDGES || 
               local_accum.gather_edges() == ALL_EDGES) {
              foreach(const edge_type& edge, in_edges) 
                local_accum.gather(context, edge);
            }
            if(local_accum.gather_edges() == OUT_EDGES || 
               local_accum.gather_edges() == ALL_EDGES) {
              foreach(const edge_type& edge, out_edges) 
                local_accum.gather(context, edge);
            }
          }

          if(graph_ptr->l_is_master(lvid)) {
            context.init(graph_ptr->global_vid(lvid), VERTEX_CONSISTENCY);
            local_accum(context);
          }
        }
        // std::cout << "Finished local sync: " 
        //           << lowres_time_millis() / 1000 << std::endl;
        // Merge with master
        lock.lock(); 
        shared_aggregator += local_accum; 
        lock.unlock();
        barrier_ptr->wait();  // Wait until all merges are complete

        if(cpuid == 0) {
          // std::cout << "Merging sync: " 
          //           << lowres_time_millis() / 1000 << std::endl;
          std::vector<aggregator_type> result(rmi_ptr->numprocs());
          result[rmi_ptr->procid()] = shared_aggregator;
          const size_t ROOT(0);
          rmi_ptr->gather(result, ROOT);
          if(rmi_ptr->procid() == ROOT) {
            // Sum up all the results
            shared_aggregator = zero;
            for(size_t i = 0; i < result.size(); ++i) 
              shared_aggregator += result[i];
            // Finalize with the global context
            iglobal_context& global_context = context;
            shared_aggregator.finalize(global_context);
            const size_t time_in_seconds = lowres_time_millis() / 1000;
            std::cout << "Sync Finished: " << time_in_seconds << std::endl;

          }
          // Zero out the shared accumulator for the next run
          shared_aggregator = zero;
          // update the sync queue
          dist_aggregator->master_lock.lock();
          dist_aggregator->schedule_prelocked(key, interval);
          dist_aggregator->sync_in_progress = false;
          rmi_ptr->barrier();
          dist_aggregator->master_lock.unlock();          
        }
      } // end of run aggregator
    }; // end of sync



  private:
    friend class isync;
    
    template <typename Aggregator>
    friend class sync;
    
    //! The base communication object
    rmi_type rmi;

    engine_type& engine;
    graph_type& graph;


    //! A lock used to manage access to internal data-structures
    mutex master_lock;
    bool sync_in_progress;
    
    //! The container of all registered sync operations
    typedef std::map<std::string, isync*> sync_map_type;
    sync_map_type sync_map;
   
    //! The priority queue over sync operations
    typedef mutable_queue<std::string, long> sync_queue_type;   
    sync_queue_type sync_queue;

    //! The pool of threads used to run the sync
    thread_pool threads;


  public:

    distributed_aggregator(distributed_control& dc, engine_type& engine,
                           graph_type& graph) :
      rmi(dc, this), engine(engine), graph(graph), 
      sync_in_progress(false) { rmi.barrier(); }


    /** Destroy members of the map */
    ~distributed_aggregator() { 
      threads.join();
      master_lock.lock();
      sync_queue.clear();
      typedef typename sync_map_type::value_type pair_type;
      foreach(pair_type& pair, sync_map) {        
        ASSERT_TRUE(pair.second != NULL);
        delete pair.second; pair.second = NULL;
      }     
      sync_map.clear();
      master_lock.unlock();
    }


    /**
     * Get a reference to the internal threads in the aggregator.
     * This is used by the engine to join threads and to resize the
     * thread pool.
     */
    thread_pool& get_threads() { return threads; }

    /**
     * This is used by the engine at start to initialize the
     * aggregation task queue.
     */
    void initialize_queue() { 
      // Only maintian a queue on processor 0
      if(rmi.procid() != 0) return;
      master_lock.lock();
      sync_queue.clear();
      typedef typename sync_map_type::value_type pair_type;
      foreach(const pair_type& pair, sync_map) {        
        ASSERT_TRUE(pair.second != NULL);
        // If their is a sync associated with the global record and
        // the sync has a non-zero interval than schedule it.
        if(pair.second->interval > 0) 
          schedule_prelocked(pair.first, pair.second->interval);
      }
      master_lock.unlock();    
    } // end of initialize queue

  

    //! \brief Registers an aggregator with the engine
    template<typename Aggregator>
    void add_aggregator(const std::string& key,            
                        const Aggregator& zero,                 
                        size_t interval) {
      isync*& sync_ptr = sync_map[key];
      // Clear the old sync and remove from scheduling queue
      if(sync_ptr != NULL) { 
        sync_queue.remove(key);
        delete sync_ptr; sync_ptr = NULL; 
      }
      ASSERT_TRUE(sync_ptr == NULL);
      // Attach a new sync type
      sync_ptr = new sync<Aggregator>(interval, zero);
    }// end of add_sync
    

    /**
     * Run an aggregator.  This must be called on all machines
     */
    void aggregate_now(const std::string& key) {
      // BUG: check that the engine is ready to sync?
      // (engine.initialize_members()?)
      typename sync_map_type::iterator iter = sync_map.find(key);
      if(iter == sync_map.end()) {
        logstream(LOG_FATAL) 
          << "Key \"" << key << "\" is not in sync map!"
          << std::endl;
        return;
      }
      // The current implementation will lead to a deadlock if called
      // from within an update function
      master_lock.lock();
      rpc_aggregate(key);
      rmi.barrier();
      threads.join();
      master_lock.unlock();    
    } // end of aggregate_now


    void evaluate_queue() {
      if(rmi.procid() != 0) return;
      // if the engine is no longer running or there is nothing in the
      // sync queue then we terminate early
      if(sync_queue.empty()) return;
      // if there is a sync in progress
      if(sync_in_progress) return;
      // Try to grab the lock if we fail just return
      if(!master_lock.try_lock()) return;
      if(sync_in_progress) { master_lock.unlock(); return; }
      // Now we are ready to evaluate the update count
      const long negated_next_ucount = sync_queue.top().second;
      ASSERT_LE(negated_next_ucount, 0);
      const size_t next_ucount = size_t(-negated_next_ucount);
      const size_t time_in_seconds = lowres_time_millis() / 1000;
      // if it is time to run the sync then spin off the threads
      if(next_ucount < time_in_seconds) { // Run the actual sync
        std::cout << "Time requested:    \t" << next_ucount << std::endl;
        std::cout << "Sync initiated at: \t" << time_in_seconds << std::endl;
        const std::string key = sync_queue.top().first;
        sync_queue.pop();
        initiate_aggregate(key);
      }    
      master_lock.unlock();    
    } // end of evaluate_sync_queue
    

    void initiate_aggregate(const std::string& key) {
      if(rmi.procid() != 0) return;
      ASSERT_FALSE(sync_in_progress);
      sync_in_progress = true;
      for(procid_t proc = 0; proc < rmi.numprocs(); ++proc) {
        rmi.remote_call(proc, &distributed_aggregator::rpc_aggregate, key);
      }
    } // end of innitiate aggregate


    void rpc_aggregate(const std::string& key) {
      master_lock.lock();    
      if(rmi.procid() != 0) {
        ASSERT_FALSE(sync_in_progress);
        sync_in_progress = true;
      } else {
        ASSERT_TRUE(sync_in_progress);
      }
      master_lock.unlock();     
      // Get the sync
      isync* sync_ptr = sync_map[key];
      ASSERT_FALSE(sync_ptr == NULL);
      // Initialize the counter
      sync_ptr->barrier_ptr = new graphlab::barrier(threads.size());
      sync_ptr->counter = 0;
      // Launch the threads
      for(size_t i = 0; i < threads.size(); ++i) {
        const boost::function<void (void)> sync_function = 
          boost::bind(&isync::run_aggregator, 
                      sync_ptr, key,  i, this, &rmi, &engine, &graph);
        threads.launch(sync_function);
      }
    } // end of rpc aggregate

  

    

    
    void schedule_prelocked(const std::string& key, size_t sync_interval) {
      const size_t time_in_seconds = lowres_time_millis() / 1000;
      const long negated_next = -long(time_in_seconds + sync_interval); 
      sync_queue.push(key, negated_next);
    } // end of schedule_prelocked
    
    
    
  }; // end of class shared memory aggregator


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif
