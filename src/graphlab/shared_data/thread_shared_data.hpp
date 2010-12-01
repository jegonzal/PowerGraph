#ifndef GRAPHLAB_THREAD_SHARED_DATA_HPP
#define GRAPHLAB_THREAD_SHARED_DATA_HPP

#include <map>
#include <limits>

#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/timer.hpp>

#include <graphlab/shared_data/ishared_data.hpp>
#include <graphlab/scope/general_scope.hpp>

#include <unistd.h>

#include <graphlab/macros_def.hpp>
namespace graphlab {

  template<typename Graph>
  class thread_shared_data : 
    public ishared_data_manager<Graph> {
    // Typedefs
    // ==============================================================>
  public:

    typedef Graph graph_type;
    typedef thread_shared_data<Graph> base;

    typedef typename base::iscope_type iscope_type;
    typedef typename base::iscope_factory_type iscope_factory_type;
    
    typedef typename base::sync_function_type sync_function_type;
    typedef typename base::apply_function_type apply_function_type;


  private:
   
    struct sync_task {
      sync_function_type sync_fun;
      apply_function_type apply_fun;
      size_t sync_interval;
      size_t next_time;
      size_t min_range;
      size_t max_range;
      any zero;
      mutex lock;
      size_t rangelow;
      size_t rangehigh;
      sync_task() :
        sync_fun(NULL), apply_fun(NULL),
        sync_interval(-1),
        next_time(0) { }
    };

    struct atomic_entry {
      any value;
      rwlock lock;
    };
    
    typedef std::vector<any> constant_map_type;
    typedef std::map<size_t, atomic_entry> atomic_map_type;
    typedef std::map<size_t, sync_task> sync_map_type;
    
  public:

    thread_shared_data() : scope_factory(NULL) {
      logger(LOG_WARNING, 
            "The use of the shared_data table has been deprecated. "
            "Please use glshared");
    }


    /** Set an immutable constant */
    void set_constant(size_t index, const any& new_value) {
      if(constants.size() <= index) constants.resize(index+1);
      constants[index] = new_value;
    } // end of 

    const any& get_constant(size_t index) const {
      assert(index < constants.size());
      assert(!constants[index].empty());
      return constants[index];
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
      sync.sync_fun = sync_fun;
      sync.apply_fun = apply_fun;
      sync.zero = zero;
      sync.sync_interval = sync_interval;
      if(sync_interval != size_t(-1) ) {
        // Force a sync on the first time around
        sync.next_time = 0;
      } else {
        sync.next_time = -1;
      }
      
      sync.rangelow = rangelow;
      sync.rangehigh = rangehigh;
      ASSERT_LE(sync.rangelow, sync.rangehigh);
      create_atomic(index, zero);
    }

    void create_atomic(size_t index, const any& initial_value) {
      atomic_entry& entry = atomic_map[index];
      entry.value = initial_value;
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
        if(stask.next_time < lowres_time_millis() &&
           stask.sync_interval != size_t(-1) ) {
          sync(index);
        }
      }
    }
    
    void signal_all() {
      typedef typename sync_map_type::value_type value_type;
      foreach(value_type& pair, sync_map) {
        size_t index = pair.first;
        sync_task& stask = pair.second;
        if(stask.next_time < lowres_time_millis() &&
           stask.sync_interval != size_t(-1) ) {
          sync(index);
        }
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
    
    
    void trigger_sync_all() {
      typedef typename sync_map_type::value_type value_type;
      foreach(const value_type& pair, sync_map) trigger_sync(pair.first);
    }
      
    void trigger_sync(size_t index) {
      typedef typename sync_map_type::iterator iterator_type;
      iterator_type iter = sync_map.find(index);
      if(iter != sync_map.end()) {
        sync_task& stask = iter->second;
        stask.next_time = lowres_time_millis();
      }
    }
    /**
     * \brief run sync on a particular field
     *
     * User sync requests should use the sync which takes a graph
     * argument.
     *
     * \todo this should be optimized
     */
    void sync(size_t index) {

      // Requires a scope factory
      assert(scope_factory != NULL);
      // Get the table entry
      typedef typename sync_map_type::iterator iterator_type;
      
      iterator_type iter = sync_map.find(index);
      assert(iter != sync_map.end());
      sync_task& sync = iter->second;
      if (sync.sync_fun == NULL && sync.apply_fun == NULL) return;
      
      if(sync.lock.try_lock()) {
        sync.next_time = -1;
        // std::cout << "Sync: " << index << std::endl;
        // Copy the accumulator
        any accumulator = sync.zero;
        if (sync.sync_fun != NULL) {
          // Try and grab the lock for the sync
          size_t cpuid = thread::thread_id();
          size_t vmax = std::min(scope_factory->num_vertices(), sync.rangehigh);
          size_t vmin = std::max(size_t(0), sync.rangelow);
          ASSERT_LE(vmin, vmax);
          for(size_t v = vmin; v < vmax; ++v) {
            // get the scope for the vertex
            iscope_type* scope = scope_factory->get_scope(cpuid, v);
            assert(scope != NULL);
            // Apply the sync function
            sync.sync_fun(index, *this, *scope, accumulator);
            // Commit and free the scope
            scope->commit();
            scope_factory->release_scope(scope);
          }
        }
        
        if (sync.apply_fun != NULL) {
          // Apply the result
          atomic_apply(index, sync.apply_fun, accumulator);        
        }
        
        // set the next sync time
        if(sync.sync_interval != size_t(-1) ) {
          size_t lrtime = lowres_time_millis();
          sync.next_time = lrtime + sync.sync_interval;
          assert(sync.next_time >= lrtime);
        }
        // release the lock
        sync.lock.unlock();  
      }
    } // end of sync

    /**
     * \brief Run all sync tasks using the graph data.
     *
     * This function can be called by user code outside of the
     * graphlab engine.
     *
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
     *
     */
    void sync(Graph& graph, size_t index) {
      // Get the table entry
      typedef typename sync_map_type::iterator iterator_type;
      iterator_type iter = sync_map.find(index);
      assert(iter != sync_map.end());
      sync_task& sync = iter->second;
      // Ensure that the sync is well defined
      assert(sync.sync_fun != NULL);
      assert(sync.apply_fun != NULL);

      if(sync.lock.try_lock()) {
        // Copy the accumulator
        any accumulator = sync.zero;
        typedef general_scope<Graph> scope_type;              
        size_t vmax = std::min(graph.num_vertices(), sync.rangehigh);
        size_t vmin = std::max(size_t(0), sync.rangelow);
        ASSERT_LE(vmin, vmax);
        for(size_t v = vmin; v < vmax; ++v) {
          // get the scope for the vertex
          scope_type scope(&graph, v,NULL);
          // Apply the sync function
          sync.sync_fun(index, *this, scope, accumulator);
          // Commit and free the scope
          scope.commit();
        }
        // Apply the result
        atomic_apply(index, sync.apply_fun, accumulator);        
        // set the next sync time
        float lrtime = lowres_time_millis();
        sync.next_time = lrtime + sync.sync_interval;
        if (sync.next_time < lrtime) { // overflow 
          sync.next_time = -1;
        }
        // release the lock
        sync.lock.unlock();

        // sync.next_time = sync.next_time + sync.sync_interval;
      } 
    } // end of sync


    any get(size_t index) const {
      // Get the field
      typedef typename atomic_map_type::const_iterator iterator_type;
      iterator_type iter = atomic_map.find(index);
      assert(iter != atomic_map.end());
      const atomic_entry& entry = iter->second;
      entry.lock.readlock();
      any res = entry.value;
      entry.lock.unlock();
      return res;
    }
    // in the shared memmory setting, get and atomic get are the same
    any atomic_get(size_t index) const {
      return get(index);
    }


    void atomic_set(size_t index, const any& data) {
      typedef typename atomic_map_type::iterator iterator_type;
      iterator_type iter = atomic_map.find(index);
      // The shared data item must already have been created
      assert(iter != atomic_map.end());
      atomic_entry& entry = iter->second;
      entry.lock.writelock();
      entry.value = data;
      entry.lock.unlock();
    }


    any atomic_exchange(size_t index, const any& data) {
      typedef typename atomic_map_type::iterator iterator_type;
      iterator_type iter = atomic_map.find(index);
      // The shared data item must already have been created
      assert(iter != atomic_map.end());
      atomic_entry& entry = iter->second;
      entry.lock.writelock();
      any old_value = entry.value;
      entry.value = data;
      entry.lock.unlock();
      return old_value;
    }


    void atomic_apply(size_t index,
                      apply_function_type fun,
                      const any& data) {
      typedef typename atomic_map_type::iterator iterator_type;
      iterator_type iter = atomic_map.find(index);
      // The shared data item must already have been created
      assert(iter != atomic_map.end());
      atomic_entry& entry = iter->second;
      entry.lock.writelock();
      fun(index, *this, entry.value, data);
      entry.lock.unlock();
    }



    /** Set the scpe factory, called by the engine */
    void set_scope_factory(iscope_factory_type* factory) {
      scope_factory = factory;
    }


    
    
    // Data Members
    // ==============================================================>
  private:
    constant_map_type constants;
    atomic_map_type atomic_map;
    sync_map_type sync_map;
    iscope_factory_type* scope_factory;

  }; // end of class thread_share_data



}; // end graphlab namespace
#include <graphlab/macros_undef.hpp>



#endif
