#ifndef DISTRIBUTED_SHARED_DATA_HPP
#define DISTRIBUTED_SHARED_DATA_HPP

#include <map>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/generics/any_vector.hpp>
namespace graphlab {

template <typename ContextType> 
class distributed_shared_data {
 public:
  struct distributed_global_record { 
    std::vector<spinlock> locks;
    graphlab::any_vector values;
    bool is_const;
  };
  
  typedef std::map<std::string, distributed_global_record> global_map_type;
  global_map_type global_records;
  
  size_t numthreads;

  // Global Aggregates ------------------------------------------------------
  struct isync {
    std::string key;
    isync(std::string key): key(key) { }
    virtual ~isync() { }
    virtual void reset() = 0;
    virtual void accumulate(ContextType* context,
                            size_t threadid) = 0;
    virtual void finalize(ContextType* context,
                dc_dist_object<distributed_shared_data<ContextType> >& rmi) = 0;
  }; // end of isync
  
  template<typename Accum >
  struct sync : public isync {
    typedef Accum       accumulator_type;
    using isync::key;
    const accumulator_type zero;
    accumulator_type shared_accumulator;
    std::vector<accumulator_type> local_accumulator;
    mutex lock;
    size_t numthreads;
    sync(const accumulator_type& zero,
         std::string key,
         size_t numthreads) : isync(key), zero(zero), 
                            shared_accumulator(zero),
                            local_accumulator(numthreads, zero),
                            numthreads(numthreads){ }
    void reset() {
      shared_accumulator = zero;
      local_accumulator = std::vector<accumulator_type>(numthreads, zero);
    }
    
    void accumulate(ContextType* context,
                    size_t threadid) {
      ASSERT_LT(threadid, local_accumulator.size());
      local_accumulator[threadid](*context);
    }
    
    void finalize(ContextType* context,
                dc_dist_object<distributed_shared_data<ContextType> >& rmi) {
      for (size_t i = 0;i < local_accumulator.size(); ++i) {
        shared_accumulator += local_accumulator[i];
      }
      // gather all accumulators to 0
      if (rmi.procid() == 0) {
        for (procid_t i = 1;i < rmi.numprocs(); ++i) {
          Accum acc;
          rmi.recv_from(i, acc);
          shared_accumulator += acc;
        }
        shared_accumulator.finalize(*context);
      }
      else {
        rmi.send_to(0, shared_accumulator);
      }
      
    }
  }; // end of sync

  std::map<std::string, isync*> sync_map;

  mutable dc_dist_object<distributed_shared_data<ContextType> > rmi;
  
  distributed_shared_data(distributed_control &dc,
                          size_t numthreads): numthreads(numthreads),rmi(dc, this) { }
  
  /**
    * Define a global mutable variable (or vector of variables).
    * Must be called on all machines simultaneously
    *
    * \param key the name of the variable (vector)
    * \param value the initial value for the variable (vector)
    * \param size the initial size of the global vector (default = 1)
    * 
    */
  template< typename T >
  void add_global(const std::string& key, const T& value, 
                  size_t size = 1) {
    distributed_global_record& record = global_records[key];
    // Set the initial value (this can change the type)
    typedef std::vector<T> vector_type;
    record.values = vector_type(size, value);
    record.is_const = false;
    record.locks.resize(size);
  }

  /**
    * Define a global constant.
    * Must be called on all machines simultaneously
    */
  template< typename T >
  void add_global_const(const std::string& key, const T& value, 
                        size_t size = 1) {
    distributed_global_record& record = global_records[key];
    // Set the initial value (this can change the type)
    typedef std::vector<T> vector_type;
    record.values = vector_type(size, value);
    record.is_const = true;
    record.locks.resize(size);
  }


  //! Change the value of a global entry
  template< typename T >
  void set_global(const std::string& key, const T& value, 
                  size_t index = 0) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;
      return;
    }
    distributed_global_record& record = iter->second;
    std::vector<T>& values = record.values.template as<T>();
    ASSERT_EQ(values.size(), record.locks.size());
    ASSERT_LT(index, values.size());
    record.locks[index].lock();
    values[index] = value; 
    record.locks[index].unlock();
  }
  
  

  //! Change the value of a global entry
  void set_global_any(const std::string& key, const any& value, 
                      size_t index = 0) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;
      return;
    }
    distributed_global_record& record = iter->second;
    ASSERT_LT(index, record.values.size());
    record.locks[index].lock();
    record.values.set(index, value); 
    record.locks[index].unlock();
  }
  
  void synchronize_global(const std::string key, size_t index) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;
      return;
    }
    distributed_global_record& record = iter->second;
    
    for (procid_t i = 0;i < rmi.numprocs(); ++i) {
      if (i != rmi.procid()) {
        rmi.remote_call(i, &distributed_shared_data::set_global_any, 
                        key, record.values.get(index), index);
      }
    }
  }

  //! Get a copy of the value of a global entry
  template< typename T >
  T get_global(const std::string& key, size_t index = 0) const {
    typename global_map_type::const_iterator iter = global_records.find(key);
    if(iter == global_records.end()) {
      logstream(LOG_FATAL) 
        << "Key \"" << key << "\" is not in global map!"
        << std::endl;      
    }
    const distributed_global_record& record = iter->second;
    typedef std::vector<T> vector_type;
    const vector_type& values = record.values.template as<T>();
    ASSERT_EQ(values.size(), record.locks.size());
    ASSERT_LT(index, values.size());
    record.locks[index].lock();
    T ret_value = values[index];
    record.locks[index].unlock();
    return ret_value;
  }

  /// Must be called on all machines simultaneously
  template<typename Accum>
  void add_sync(const std::string& key,            
                const Accum& zero) {
    isync*& sync_ptr = sync_map[key];
    // Clear the old sync and remove from scheduling queue
    if(sync_ptr != NULL) { delete sync_ptr; sync_ptr = NULL; }
    ASSERT_TRUE(sync_ptr == NULL);
    // Attach a new sync type
    sync_ptr = new sync<Accum>(zero, key, numthreads);
  }

  void reset_sync(std::string key) {
    sync_map[key]->reset();
  }

  void reset_all_syncs() {
    typename std::map<std::string, isync*>::iterator iter = sync_map.begin();
    while (iter != sync_map.end()) {
      iter->second->reset();
      ++iter;
    }
  }

  void wait_for_all_communication() {
    rmi.full_barrier();
  }

  void accumulate(std::string key, 
                  ContextType* context,
                  size_t threadid) {
    typename std::map<std::string, isync*>::iterator iter = sync_map.find(key);
    ASSERT_TRUE(iter != sync_map.end());
    iter->second->accumulate(context, threadid);
  }
  
  /// Must be called on all machines simultaneously
  void finalize(std::string key, 
                ContextType* context) {
    typename std::map<std::string, isync*>::iterator iter = sync_map.find(key);
    ASSERT_TRUE(iter != sync_map.end());
    iter->second->finalize(context, rmi);
  }
  
  void acquire_lock(const std::string& key, size_t index = 0) {
    typename global_map_type::iterator iter = global_records.find(key);
    ASSERT_TRUE(iter != global_records.end());
    // Get the global record
    distributed_global_record& rec = iter->second;
    ASSERT_LT(index, rec.locks.size());
    rec.locks[index].lock();
  }

  void get_global(const std::string& key, 
                          any_vector*& ret_vec_ptr,
                          bool& ret_is_const) {
    typename global_map_type::iterator iter = global_records.find(key);
    if(iter == global_records.end()) ret_vec_ptr = NULL;
    else {
      // Get the global record
      distributed_global_record& rec = iter->second;
      ret_is_const = rec.is_const;      
      ret_vec_ptr = &rec.values;
    }
  }

  void commit_change(const std::string& key, size_t index = 0) {
    synchronize_global(key, index);
  }

  void release_lock(const std::string& key, size_t index = 0) {
    typename global_map_type::iterator iter = global_records.find(key);
    ASSERT_TRUE(iter != global_records.end());
    // Get the global record
    distributed_global_record& rec = iter->second;
    ASSERT_LT(index, rec.locks.size());
    rec.locks[index].unlock();
  }

};

}
#endif