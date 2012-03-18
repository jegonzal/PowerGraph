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


#ifndef GRAPHLAB_SHARED_MEMORY_AGGREGATOR
#define GRAPHLAB_SHARED_MEMORY_AGGREGATOR

#include <graphlab/logger/assertions.hpp>
#include <graphlab/parallel/pthread_tools.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename Engine>
  class shared_memory_aggregator {
  public:
    typedef Engine engine_type;
    typedef typename engine_type::graph_type graph_type;
    typedef typename engine_type::update_functor_type update_functor_type;
    typedef typename engine_type::context_type context_type;

    // Global Aggregates ------------------------------------------------------
    /**
     * The base class for registered aggregation operation
     */
    struct isync {
      size_t interval;
      vertex_id_type begin_vid, end_vid;
      virtual ~isync() { }
      virtual void run_aggregator(engine_type* engine_ptr,                          
                                  graphlab::barrier* barrier_ptr,
                                  size_t ncpus, size_t cpuid) = 0;
    }; // end of isync

    /**
     * The derived class used to store the various aggregation
     * operations
     */
    template<typename Aggregator >
    struct sync : public isync {
      typedef Aggregator       aggregator_type;
      using isync::begin_vid;
      using isync::end_vid;     
      const aggregator_type zero;
      aggregator_type shared_aggregator;
      mutex lock;
      sync(const aggregator_type& zero) : zero(zero), 
                                          shared_aggregator(zero) { }
      void run_aggregator(engine_type* engine_ptr,
                          graphlab::barrier* barrier_ptr,
                          size_t ncpus, size_t cpuid) { 
        // Thread zero must initialize the the final shared accumulator
        const size_t nverts = engine_ptr->graph.num_vertices();
        // Compute the true begin and end.
        const size_t global_begin = std::min(nverts, size_t(begin_vid));;
        const size_t global_end = std::min(nverts, size_t(end_vid));
        ASSERT_LE(global_begin, global_end);
        ASSERT_LE(global_end, nverts);
        // Compute the span of each subtask.  The span should not be
        // less than some minimal span.
        const vertex_id_type span = 
          std::max((global_end - global_begin)/ncpus, size_t(1));
        const size_t true_begin_vid = std::min(cpuid*span, nverts);
        const size_t true_end_vid   = std::min((cpuid+1)*span, nverts);
    
        // construct the local (to this thread) accumulator and
        // context
        aggregator_type local_accum(zero);
        // Do map computation;
        for(vertex_id_type vid = true_begin_vid; vid < true_end_vid; ++vid) {
          engine_ptr->evaluate_update(vid, local_accum);
        }
        // Merge with master
        lock.lock(); 
        shared_aggregator += local_accum; 
        lock.unlock();
        barrier_ptr->wait();  // Wait until all merges are complete
        if(cpuid == 0) {
          // Recast the context as a global context.  This ensures
          // that the user implements finalize correctly;
          context_type context = engine_ptr->get_context(0);
          iglobal_context& global_context = context;
          shared_aggregator.finalize(global_context);
          // Zero out the shared accumulator for the next run
          shared_aggregator = zero;
        }
        barrier_ptr->wait();
      } // end of run aggregator
    }; // end of sync


  private:
    //! The aggregator this engine is operating on. 
    engine_type& engine;

    //! A lock used to manage access to internal data-structures
    mutex sync_master_lock;
    
    //! The container of all registered sync operations
    typedef std::map<std::string, isync*> sync_map_type;
    sync_map_type sync_map;
   
    //! The priority queue over sync operations
    typedef mutable_queue<std::string, long> sync_queue_type;   
    sync_queue_type sync_queue;

    //! The pool of threads used to run the sync
    thread_pool threads;


  public:

    shared_memory_aggregator(engine_type& engine) :
      engine(engine) { }

    /** Destroy members of the map */
    ~shared_memory_aggregator() { 
      threads.join();
      sync_master_lock.lock();
      sync_queue.clear();
      typedef typename sync_map_type::value_type pair_type;
      foreach(pair_type& pair, sync_map) {        
        ASSERT_TRUE(pair.second != NULL);
        delete pair.second; pair.second = NULL;
      }     
      sync_map.clear();
      sync_master_lock.unlock();
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
      sync_master_lock.lock();
      sync_queue.clear();
      typedef typename sync_map_type::value_type pair_type;
      foreach(const pair_type& pair, sync_map) {        
        ASSERT_TRUE(pair.second != NULL);
        // If their is a sync associated with the global record and
        // the sync has a non-zero interval than schedule it.
        if(pair.second->interval > 0) 
          schedule_sync_prelocked(pair.first, pair.second->interval);
      }
      sync_master_lock.unlock();    
    } // end of initialize queue


  

    //! \brief Registers an aggregator with the engine
    template<typename Aggregator>
    void add_aggregator(const std::string& key,            
                        const Aggregator& zero,                 
                        size_t interval,
                        vertex_id_type begin_vid = 0,
                        vertex_id_type end_vid = 
                        std::numeric_limits<vertex_id_type>::max()) {
      isync*& sync_ptr = sync_map[key];
      // Clear the old sync and remove from scheduling queue
      if(sync_ptr != NULL) { delete sync_ptr; sync_ptr = NULL; }
      sync_queue.remove(key);
      ASSERT_TRUE(sync_ptr == NULL);
      // Attach a new sync type
      typedef sync<Aggregator> sync_type;
      sync_ptr = new sync_type(zero);
      sync_ptr->interval    = interval;
      sync_ptr->begin_vid   = begin_vid;
      sync_ptr->end_vid     = end_vid;
    }// end of add_sync
    

    //! Performs a sync immediately.
    void aggregate_now(const std::string& key) {
      engine.initialize_members();    
      typename sync_map_type::iterator iter = sync_map.find(key);
      if(iter == sync_map.end()) {
        logstream(LOG_FATAL) 
          << "Key \"" << key << "\" is not in sync map!"
          << std::endl;
        return;
      }
      isync* sync = iter->second;
      ASSERT_NE(sync, NULL);
      // The current implementation will lead to a deadlock if called
      // from within an update function
      sync_master_lock.lock();
      aggregate_prelocked(key, sync);
      sync_master_lock.unlock();
    } // end of sync_now


    void evaluate_queue(size_t last_update_count) {
      // if the engine is no longer running or there is nothing in the
      // sync queue then we terminate early
      if(sync_queue.empty()) return;
      // Try to grab the lock if we fail just return
      if(!sync_master_lock.try_lock()) return;
      // ASSERT: the lock has been aquired. Test for a task at the top
      // of the queue
      const long negated_next_ucount = sync_queue.top().second;
      ASSERT_LE(negated_next_ucount, 0);
      const size_t next_ucount = size_t(-negated_next_ucount);
      // if we have more updates than the next update count for this
      // task then run it
      if(next_ucount < last_update_count) { // Run the actual sync
        const std::string key = sync_queue.top().first;
        sync_queue.pop();
        isync* sync = sync_map[key];
        ASSERT_NE(sync, NULL);
        aggregate_prelocked(key, sync);
        // Reschedule the sync record
        schedule_sync_prelocked(key, sync->interval);
      }    
      sync_master_lock.unlock();    
    } // end of evaluate_sync_queue
    

    void aggregate_prelocked(const std::string& key, isync* sync) { 
      ASSERT_NE(sync, NULL);
      graphlab::barrier barrier(threads.size());
      for(size_t i = 0; i < threads.size(); ++i) {
        const boost::function<void (void)> sync_function = 
          boost::bind(&(isync::run_aggregator), 
                      sync, &engine, &barrier, threads.size(), i);
        threads.launch(sync_function);
      }
      engine.join_threads(threads);
    } // end of launch sync prelocked
  

    

    
    void schedule_sync_prelocked(const std::string& key, size_t sync_interval) {
      const size_t ucount = engine.last_update_count();
      const long negated_next_ucount = -long(ucount + sync_interval); 
      sync_queue.push(key, negated_next_ucount);
    } // end of schedule_sync
    
    
    
  }; // end of class shared memory aggregator


}; // end of graphlab namespace
#include <graphlab/macros_undef.hpp>

#endif
