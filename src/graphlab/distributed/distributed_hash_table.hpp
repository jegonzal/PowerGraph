/*
  \author Yucheng Low (ylow)
  An implementation of a distributed integer -> integer map with caching
  capabilities. 

*/

#ifndef DISTRIBUTED_HASH_TABLE_HPP
#define DISTRIBUTED_HASH_TABLE_HPP
#include <boost/unordered_map.hpp>
#include <boost/intrusive/list.hpp>

#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/distributed/distributed_control.hpp>
#include <graphlab/util/generics/any.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/synchronized_unordered_map.hpp>
#include <graphlab/util/dense_bitset.hpp>

namespace graphlab {

/**
  A cache entry. Boost intrusive is used to provide the LRU capabilities here
*/
class any_lru_list{
 public:

  size_t key; /// the key assiciated with this cache entry
  any value; /// the value assiciated with this cache entry
  typedef boost::intrusive::list_member_hook<
            boost::intrusive::link_mode<boost::intrusive::auto_unlink> >
                                                          lru_member_hook_type;

  lru_member_hook_type member_hook_;
  ~any_lru_list() { }
  explicit any_lru_list(size_t k = 0, const any &i = any()) : key(k) {
    value = i;
  }
};

/**
This implements a distributed size_t => size_t map with caching capabilities
*/
class distributed_hash_table{
 public:

  /// datatype of the data map
  typedef synchronized_unordered_map<std::pair<rwlock, any> > map_type;
  /// datatype of the local cache map
  typedef boost::unordered_map<size_t, any_lru_list*> cache_type;


  typedef boost::intrusive::member_hook<any_lru_list,
                                        any_lru_list::lru_member_hook_type,
                                        &any_lru_list::member_hook_> MemberOption;
  /// datatype of the intrusive LRU list embedded in the cache map
  typedef boost::intrusive::list<any_lru_list, MemberOption, boost::intrusive::constant_time_size<false> > lru_list_type;

  /// Constructor. Creates the integer map.
  distributed_hash_table(distributed_control &dc, size_t max_cache_size = 1024);

  void set_pushed_updates(bool _pushed_updates) {
    pushed_updates = _pushed_updates;
  }

  ~distributed_hash_table();
  /// Sets the key to the value
  void set(size_t key, const any &value);

    /// Sets the key to the value
  void exchange(size_t key, const any &newvalue, any &oldvalue);

  void modified(size_t key);

  /** Gets the value associated with the key. returns true on success.
      The any returned is a copy. */
  bool get(size_t key, any &value);

  /** Gets the value associated with the key, reading from cache if available
      Note that the cache may be out of date.
      The any returned is a copy. */
  bool get_cached(size_t key, any &value);

  template <typename T>
  bool get(size_t key, T &value) {
    any t;
    bool ret = get(key, t);
    if (ret == false) return ret;
    value = t.as<T>();
    return true;
  }


  template <typename T>
  bool get_cached(size_t key, T &value) {
    any t;
    bool ret = get_cached(key, t);
    if (ret == false) return ret;
    value = t.as<T>();
    return true;
  }

  template <typename T>
  void exchange(size_t key, const T &newvalue, T &oldvalue) {
    any newval = newvalue;
    any oldval = oldvalue;
    exchange(key, newval, oldval);
  }

  /// Invalidates the cache entry associated with this key
  void invalidate(size_t key);

  double cache_miss_rate();

  size_t num_gets() const {
    return reqs;
  }
  size_t num_misses() const {
    return misses;
  }

  size_t cache_size() const {
    return cache.size();
  }

  struct mapreplydata{
    bool reply_set;
    bool reply_found;
    any reply_value;
  };

 private:

  distributed_control &dc;
  map_type data;  /// The actual table data that is distributed

  
  spinlock cachelock; /// lock for the cache datastructures
  cache_type cache;   /// The cache table
  lru_list_type lruage; /// THe LRU linked list associated with the cache


  procid_t numprocs;   /// NUmber of processors
  size_t maxcache;     /// Maximum cache size allowed

  size_t mapid;

  bool pushed_updates;

  /// Updates the cache with this new value
  void update_cache(size_t key, const any &val);

  /// Removes the least recently used element from the cache
  void remove_lru();

  size_t reqs;
  size_t misses;

  /// The key to nodeid hash function
  inline procid_t key_node_hash(size_t keyval) {
    return keyval % numprocs;
  }

  friend void int_any_map_get_handler(distributed_control& dc, 
                                      procid_t source,  
                                      void* ptr,    //unused
                                      size_t len,   //unused
                                      handlerarg_t mapid,
                                      handlerarg_t key,
                                      handlerarg_t reqid);
  
  friend void int_any_map_set_handler(distributed_control& dc, 
                                      procid_t source,  
                                      void* ptr,    //serialized any
                                      size_t len,   
                                      handlerarg_t mapid,
                                      handlerarg_t key,
                                      handlerarg_t exchangereqid,
                                      handlerarg_t setreplyreqid);
  
  friend void int_any_map_get_reply_handler(distributed_control& dc, 
                                            procid_t source,  
                                            void* ptr,    //serialized any
                                            size_t len,   
                                            handlerarg_t mapid,
                                            handlerarg_t key,
                                            handlerarg_t reqid);

  friend void cache_write_handler(distributed_control& dc,
                                procid_t source,
                                void* ptr,     //serialized any
                                size_t len,  
                                handlerarg_t mapid,
                                handlerarg_t key);


  template <typename Graph>
  friend class distributed_shared_data;
  template <typename Graph>
  friend class distributed_fullsweep_sdm;
// 
};

}
#endif
