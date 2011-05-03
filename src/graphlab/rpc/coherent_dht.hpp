/*
This file is part of GraphLab.

GraphLab is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

GraphLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with GraphLab.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef COHERENT_DHT_HPP
#define COHERENT_DHT_HPP
#include <boost/unordered_map.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/function.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/util/timer.hpp>
namespace graphlab {

namespace dc_impl {
  
#define COHERENT_DHT_COMPRESSED_HASH 32768
#define COHERENT_DHT_SUBSCRIBE_IF_ACCESSES_PER_INVALIDATE 10
/**
 * \ingroup rpc_internal
  A cache entry for the coherent_dht. 
  Boost intrusive is used to provide the LRU capabilities here
*/
template<typename KeyType, typename ValueType>
class coherent_lru_list{
 public:

  KeyType key; /// the key assiciated with this cache entry
  ValueType value; /// the value associated with this cache entry
  uint32_t accesses;
  typedef boost::intrusive::list_member_hook<
            boost::intrusive::link_mode<boost::intrusive::auto_unlink> >
                                                          lru_member_hook_type;

  lru_member_hook_type member_hook_;
  ~coherent_lru_list() { }
  explicit coherent_lru_list(const KeyType& k = KeyType(), const ValueType &v = ValueType()) 
          : key(k), value(v),accesses(0) {  }
};

} // namespace dc_impl

/**
 * \ingroup rpc
This implements a processor consistent cache coherent distributed hash table.
Each machine has a part of the hash table as well as a cache. The system
implements automatic cache invalidation as well as automatic cache subscription
(currently through a rather poor heuristic). 
This class also implements some functionality which provide more fine-grained 
control over the invalidation policy, allowing stronger consistency levels to
be implemented on top of this class.
*/
template<typename KeyType, typename ValueType>
class coherent_dht{
 public:

  typedef dc_impl::coherent_lru_list<KeyType, ValueType> lru_entry_type;
  /** datatype of the data map. maps from key to the value */
  typedef boost::unordered_map<KeyType, ValueType> map_type;
  /// datatype of the local cache map
  typedef boost::unordered_map<KeyType, lru_entry_type* > cache_type;


  typedef boost::intrusive::member_hook<lru_entry_type,
                                        typename lru_entry_type::lru_member_hook_type,
                                        &lru_entry_type::member_hook_> MemberOption;
  /// datatype of the intrusive LRU list embedded in the cache map
  typedef boost::intrusive::list<lru_entry_type, 
                                 MemberOption, 
                                 boost::intrusive::constant_time_size<false> > lru_list_type;

  /// Constructor. Creates the integer map.
  coherent_dht(distributed_control &dc, 
                         size_t max_cache_size = 1024):rpc(dc, this),data(11) {
    has_mod_trigger = false;
                           
    cache.rehash(max_cache_size);
    maxcache = max_cache_size;
    logger(LOG_INFO, "%d Creating distributed_hash_table. Cache Limit = %d", 
           dc.procid(), maxcache);
    reqs = 0;
    misses = 0;
    
    for (size_t i = 0;i < COHERENT_DHT_COMPRESSED_HASH; ++i) {
      subscription[i].resize(dc.numprocs());
      subscription[i].clear();
    }
    dc.barrier();
  }


  ~coherent_dht() {
    data.clear();
    typename cache_type::iterator i = cache.begin();
    while (i != cache.end()) {
      delete i->second;
      ++i;
    }
    cache.clear();
  }
  
  /**
  Attaches a modification trigger which is signalled whenever
  the local cache/storage of any index is updated or invalidated.
  Only one modification trigger can be attached.
  The trigger is only signalled after the update/invalidation 
  is complete.
  The trigger is called with three parameters, the key, the value as well
  as a boolean flag "is_in_cache". 
  
  The "key" is the key of the entry which was just updated/invalidated.
  The "value" is the current value of the entry. This is only set if is_in_cache is true
  Note that the call to the trigger is not locked
  and the "is_in_cache" flag could very well be outdated when the call is issued.
  The flag should therefore not be treated as "truth" but simply as a hint
  about the state of the internal cache.

  */
  void attach_modification_trigger(boost::function<void(const KeyType&, const ValueType&, bool)> trigger) {
    ASSERT_FALSE(has_mod_trigger);
    has_mod_trigger = true;
    modification_trigger = trigger;
  }
  
  /**
  Detaches the modification trigger
  */
  void detach_modification_trigger() {
    has_mod_trigger = false;
  }
  
  /** acquire a lock on the key.
  The lock should be released using end_critical_section() as 
  soon as possible. There is no guarantee that the lock on this key
  is fine-grained.
  */
  ValueType& begin_critical_section(const KeyType &key) {
    ASSERT_EQ(owning_machine(key), rpc.procid());
    // get the compressed hash
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    // acquire the fine grained lock
    finegrained_lock[compressedhash].lock();
    typename map_type::iterator iter = data.find(key);
    assert(iter != data.end());
    return iter->second;
  }
  
  /**
    Releases a lock on the key as a acquired by begin_critical_section()
  */
  void end_critical_section(const KeyType &key) {
    ASSERT_EQ(owning_machine(key), rpc.procid());
    // get the compressed hash
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    // release the fine grained lock
    finegrained_lock[compressedhash].unlock();
  }  
  
  
  void set(const KeyType& key, const ValueType &newval)  {
    set_impl(key, newval, rpc.procid());
  }

   /** Sets the key to the value
    * if the key belongs to a remote machine.
    * It is guaranteed that if the current machine sets a key to a new value
    * subsequent reads will never return the previous value. (i.e. it will
    * return the new value or later values set by other processors).
    */
  void set_impl(const KeyType& key, const ValueType &newval, 
                 procid_t source)  {
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    size_t owningmachine = hashvalue % rpc.dc().numprocs();
    if (owningmachine == rpc.dc().procid()) {
      // use the full lock to get the iterator
      datalock.lock();
      typename map_type::iterator iter = data.find(key);
      // if the key does not exist, create the key.
      if (iter == data.end()) {
        finegrained_lock[compressedhash].lock();
        data[key] = newval;
        finegrained_lock[compressedhash].unlock();
        datalock.unlock();
      }
      else {
        // the key exists! switch to the fine grained lock and set the 
        // value
        datalock.unlock();
        finegrained_lock[compressedhash].lock();
        iter->second = newval;
        finegrained_lock[compressedhash].unlock();
      }
      push_changes(key, true, source);
    }
    else {
      rpc.remote_call(owningmachine, 
                      &coherent_dht<KeyType,ValueType>::set_impl, 
                      key,
                      newval,
                      source);
      update_cache(key, newval, false);
    }
  }


  void set_synchronous(const KeyType& key, const ValueType &newval) {
    set_synchronous_impl(key, newval, rpc.procid());
  }
  
  /**
  Forces synchronization of this key
  This operation is synchronous. When this function returns
  all machines are guarnateed to have the updated value
  */
  void set_synchronous_impl(const KeyType& key, const ValueType &newval,
                       procid_t source) {
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    size_t owningmachine = hashvalue % rpc.dc().numprocs();
    if (owningmachine == rpc.dc().procid()) {
      // use the full lock to get the iterator
      datalock.lock();
      typename map_type::iterator iter = data.find(key);
      // if the key does not exist, create the key.
      if (iter == data.end()) {
        finegrained_lock[compressedhash].lock();
        data[key] = newval;
        finegrained_lock[compressedhash].unlock();
        datalock.unlock();
      }
      else {
        // the key exists! switch to the fine grained lock and set the 
        // value
        datalock.unlock();
        finegrained_lock[compressedhash].lock();
        iter->second = newval;
        finegrained_lock[compressedhash].unlock();
      }
      push_changes(key, false, source);
    }
    else {
      rpc.remote_request(owningmachine, 
                        &coherent_dht<KeyType,ValueType>::set_synchronous_impl, 
                        key,
                        newval,
                        source);
      update_cache(key, newval, false);
    }
  }
  
  /**
  Push the current value of the key to all machines.
  If async=true, when this call returns, all machines are guaranteed to have
  the most up to date value of the key.
  */
  void push_changes(const KeyType& key, bool async, procid_t ignoreproc) {
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    size_t owningmachine = hashvalue % rpc.dc().numprocs();

    if (owningmachine == rpc.procid()) {
      // get a copy of the data to call the modification trigger with
      // local modification trigger,
      if (has_mod_trigger) {
        finegrained_lock[compressedhash].lock();
        typename map_type::const_iterator iter = data.find(key);
        assert(iter != data.end());
        ValueType v = iter->second;
        finegrained_lock[compressedhash].unlock();

        modification_trigger(key,v, true);
      }
      // switch to finegrained lock
      finegrained_lock[compressedhash].lock();
      typename map_type::iterator iter = data.find(key);
      finegrained_lock[compressedhash].unlock();
      assert(iter != data.end());
      if (async) {
        update_cache_coherency_set(compressedhash, key, iter->second, ignoreproc);  
      }
      else {
        update_cache_coherency_set_synchronous(compressedhash, key, iter->second, ignoreproc);
      }
    }    
    else {
      // key is not on this machine. Get the owning machine to do it
      if (async) {
        rpc.remote_call(owningmachine,
                        &coherent_dht<KeyType, ValueType>::push_changes,
                        key,
                        async,
                        ignoreproc);
      }
      else {
        rpc.remote_request(owningmachine,
                          &coherent_dht<KeyType, ValueType>::push_changes,
                          key,
                          async,
                          ignoreproc);
      }
    }
  }

  /** Gets the value associated with the key. returns true on success.
   *  get will read from the cache if data is already available in the cache.
   * If not, get will obtain the data from across the network
   */
  std::pair<bool, ValueType> get(const KeyType &key) const {
    // if this is to my current machine, just get it and don't go to cache
    procid_t owningmachine = owning_machine(key);
    if (owningmachine == rpc.dc().procid()) return get_non_cached(key);
      
    reqs++;
    cachelock.lock();
    // check if it is in the cache
    typename cache_type::iterator i = cache.find(key);
    if (i == cache.end()) {
      // nope. not in cache. Call the regular get
      cachelock.unlock();
      misses++;
      return get_non_cached(key);
    }
    else {
      // yup. in cache. return the value
      std::pair<bool, ValueType> ret;
      ret.first = true;
      ret.second = i->second->value;
      i->second->accesses++;
      // shift the cache entry to the head of the LRU list
      lruage.erase(lru_list_type::s_iterator_to(*(i->second)));
      lruage.push_front(*(i->second));
      cachelock.unlock();
      return ret;
    }
  }

  /**
    Returns the machine responsible for storing the key
  */
  procid_t owning_machine(const KeyType &key) const {
    size_t hashvalue = hasher(key);
    size_t owningmachine = hashvalue % rpc.dc().numprocs();
    return owningmachine;
  }
   
  /**
  Returns true of the key is current in the cache
  */
  bool in_cache(const KeyType &key) const {
    // if this is to my current machine, just get it and don't go to cache
    procid_t owningmachine = owning_machine(key);
    if (owningmachine == rpc.dc().procid()) return true;
    cachelock.lock();
    // check if it is in the cache
    typename cache_type::iterator i = cache.find(key);
    if (i != cache.end()) {
      cachelock.unlock();
      return true;
    }
    cachelock.unlock();
    return false;
  }
   
  /**
    Puts out a prefetch request for this key.
  */
  bool asynchronous_get(const KeyType &key) const {
    if (in_cache(key)) return true;

    procid_t owningmachine = owning_machine(key);
    rpc.remote_call(owningmachine,
                   &coherent_dht<KeyType,ValueType>::asychronous_get_handler,
                   key,
                   rpc.procid());
    return false;
  }

  /// Returns the number of misses divided by the number of requests
  double cache_miss_rate() {
    return double(misses) / double(reqs);
  }
  /// Returns the number of requests 
  size_t num_gets() const {
    return reqs;
  }
  /// Returns the number of cache misses
  size_t num_misses() const {
    return misses;
  }
  /// Returns the current size of the cache
  size_t cache_size() const {
    return cache.size();
  }
  
  /**
  Subscribes to this key. This key will be a permanent entry
  in the cache and can not be invalidated. Key modifications
  are automatically sent to this machine.
  */
  void subscribe(const KeyType &key, bool async = false) const{
    procid_t owningmachine = owning_machine(key);
    // do not subscribe if this is my machine
    if (owningmachine == rpc.dc().procid()) return;
    if (async) {
      rpc.remote_call(owningmachine,
                  &coherent_dht<KeyType,ValueType>::register_subscription, 
                  key,
                  rpc.dc().procid());
    }
    else {
      rpc.remote_request(owningmachine,
                  &coherent_dht<KeyType,ValueType>::register_subscription_synchronous, 
                  key,
                  rpc.dc().procid());
    }
  }

  void full_barrier() {
    rpc.full_barrier();
  }
  
  /// Invalidates the cache entry associated with this key
  void invalidate(const KeyType &key) const{
    bool haschanges = false;
    bool isincache = false;
    cachelock.lock();
    // is the key I am invalidating in the cache?
    typename cache_type::iterator i = cache.find(key);
    if (i != cache.end()) {
      // drop it from the lru list
      // if it is frequently accessed, don't invalidate it but subscribe.
      if (i->second->accesses >= COHERENT_DHT_SUBSCRIBE_IF_ACCESSES_PER_INVALIDATE) {
        subscribe(key, false);
        
        isincache = true;
        haschanges = true;
      }
      else {
        delete i->second;
        haschanges = true;
        cache.erase(i);
      }
    }
    cachelock.unlock();
    if (haschanges && has_mod_trigger) {
      if (isincache == false) {
        modification_trigger(key, ValueType(), isincache);
      }
      else {
      
        size_t hashvalue = hasher(key);
        size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;

        finegrained_lock[compressedhash].lock();
        typename map_type::const_iterator iter = data.find(key);
        assert(iter != data.end());
        ValueType v = iter->second;
        finegrained_lock[compressedhash].unlock();
        
        modification_trigger(key, v, isincache);
      }
    }
  }
 private:

  mutable dc_dist_object<coherent_dht<KeyType, ValueType> > rpc;
  
  map_type data;  /// The actual table data that is distributed
  mutex datalock;
  


  mutex finegrained_lock[COHERENT_DHT_COMPRESSED_HASH];  
  dense_bitset subscription[COHERENT_DHT_COMPRESSED_HASH];
                                

  
  mutex cachelock; /// lock for the cache datastructures
  mutable cache_type cache;   /// The cache table
  mutable lru_list_type lruage; /// THe LRU linked list associated with the cache


  size_t maxcache;     /// Maximum cache size allowed

  mutable size_t reqs;
  mutable size_t misses;
  


  boost::hash<KeyType> hasher;
  
  bool has_mod_trigger;
  boost::function<void(const KeyType&, const ValueType&,  bool)> modification_trigger;

  void update_cache_from_remote(const KeyType &key, const ValueType &val) const {
    return update_cache(key, val);
  }
  /** Updates the internal cache with this new value. The cache
   * entry is also moved to the head of the LRU list
  */
  void update_cache(const KeyType &key, const ValueType &val, bool calltrigger = true) const{

    cachelock.lock();
    typename cache_type::iterator i = cache.find(key);
    // create a new entry
    if (i == cache.end()) {
      cachelock.unlock();
      // if we are out of room, remove the lru entry
      if (cache.size() >= maxcache) remove_lru();
      cachelock.lock();
      // insert the element, remember the iterator so we can push it
      // straight to the LRU list
      std::pair<typename cache_type::iterator, bool> ret = cache.insert(std::make_pair(key, new lru_entry_type(key, val)));
      if (ret.second)  lruage.push_front(*(ret.first->second));
    }
    else {
        // modify entry in place
        i->second->value = val;
        // swap to front of list
        //boost::swap_nodes(lru_list_type::s_iterator_to(i->second), lruage.begin());
        lruage.erase(lru_list_type::s_iterator_to(*(i->second)));
        lruage.push_front(*(i->second));
    }
    cachelock.unlock();
     if (has_mod_trigger && calltrigger) modification_trigger(key,val, true);
  }

  
  /// Removes the least recently used element from the cache
  void remove_lru() const{
    cachelock.lock();
    KeyType keytoerase = lruage.back().key;
    // is the key I am invalidating in the cache?
    typename cache_type::iterator i = cache.find(keytoerase);
    if (i != cache.end()) {
      // drop it from the lru list
      delete i->second;
      cache.erase(i);
    }
    cachelock.unlock();
  }





  
  /**
   * Gets the true value of this key
   */
  std::pair<bool, ValueType> get_non_cached(const KeyType &key) const {
    // figure out who owns the key
    procid_t owningmachine = owning_machine(key);    
    std::pair<bool, ValueType> ret;
    // if I own the key, get it from the map table
    if (owningmachine == rpc.dc().procid()) {
      datalock.lock();
      typename map_type::const_iterator iter = data.find(key);    
      datalock.unlock();
      if (iter == data.end()) {
        ret.first = false;
      }
      else {
        ret.first = true;
        ret.second = iter->second;
      }
    }
    else {
      ret = rpc.remote_request(owningmachine, 
                               &coherent_dht<KeyType,ValueType>::get_non_cached, 
                               key);
      if (ret.first) update_cache(key, ret.second);
      else invalidate(key);
    }
    return ret;
  }


  /**
   * Called when the value changes. This updates the cache of all machines.
   * All machines subscribed to this key get updated, while all
   * machines not subscribed to this key get invalidated
   */
  void update_cache_coherency_set(size_t compressedhash, 
                                  const KeyType &key,
                                  const ValueType &value,
                                  procid_t except = procid_t(-1)) {
    // broadcast invalidate
    for (procid_t i = 0;i < rpc.dc().numprocs(); ++i) {
      if (i != rpc.dc().procid() && i != except) {
        if (subscription[compressedhash].get(i)) {
          rpc.remote_call(i, 
                          &coherent_dht<KeyType,ValueType>::update_cache_from_remote, 
                          key, value);
        }
        else {
          rpc.remote_call(i, &coherent_dht<KeyType,ValueType>::invalidate, key);
        }
      }
    }
  }

  void invalidate_reply(const KeyType &key, 
                        procid_t source, size_t reply) const {
    invalidate(key);
    if (source != procid_t(-1)) {
      rpc.dc().remote_call(source, reply_increment_counter,
                            reply, dc_impl::blob());
    }
  }

  void update_cache_reply(const KeyType &key, const ValueType &value, 
                          procid_t source, size_t reply) const {
    update_cache(key, value);
    if (source != procid_t(-1)) {
      rpc.dc().remote_call(source, reply_increment_counter,
                            reply, dc_impl::blob());
    }
  }

  void update_cache_coherency_set_synchronous(size_t compressedhash, 
                                              const KeyType &key,
                                              const ValueType &value,
                                              procid_t except = procid_t(-1)) {
    // broadcast invalidate
    dc_impl::reply_ret_type repret(true, rpc.numprocs() - 1);
    if (except < rpc.numprocs() && except != rpc.procid()) repret.flag.dec();
    
    size_t r = reinterpret_cast<size_t>(&repret); 
    for (procid_t i = 0;i < rpc.numprocs(); ++i) {
      if (i != rpc.procid() && i != except) {
        if (subscription[compressedhash].get(i)) {
          rpc.remote_call(i,
                          &coherent_dht<KeyType,ValueType>::update_cache_reply, 
                          key, value,                           
                          rpc.procid(), r);
        }
        else {
          rpc.remote_call(i, &coherent_dht<KeyType,ValueType>::invalidate_reply, 
                          key, rpc.procid(), r);
        }
      }
    }
    repret.wait();
  }
  
  void asychronous_get_handler(const KeyType &key, procid_t source) {
    std::pair<bool, ValueType> ret = get_non_cached(key);
    if (ret.first) {
      rpc.remote_call(source,
                      &coherent_dht<KeyType, ValueType>::update_cache_reply,
                      key, ret.second,
                      procid_t(-1),
                      0);
    }
  }
  
  /**
   * Subscribes the source machine to this key.
   * Naturally this key must exist.
   * We send an update to the source machine upon subscription
   */
  void register_subscription(const KeyType &key, procid_t source) {
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    subscription[compressedhash].set_bit(source);
    
    std::pair<bool, ValueType> val = get_non_cached(key);
    ASSERT_TRUE(val.first);
    rpc.remote_call(source, 
                    &coherent_dht<KeyType,ValueType>::update_cache_from_remote, 
                    key, val.second);
  }
  
    /**
   * Subscribes the source machine to this key.
   * Naturally this key must exist.
   * We send an update to the source machine upon subscription
   */
  void register_subscription_synchronous(const KeyType &key, procid_t source) {
    size_t hashvalue = hasher(key);
    size_t compressedhash = hashvalue % COHERENT_DHT_COMPRESSED_HASH;
    subscription[compressedhash].set_bit(source);
    
    std::pair<bool, ValueType> val = get_non_cached(key);
    ASSERT_TRUE(val.first);
    rpc.remote_request(source, 
                        &coherent_dht<KeyType,ValueType>::update_cache_from_remote, 
                        key, val.second);
  }
};

}
#endif
