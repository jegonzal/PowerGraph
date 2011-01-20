#ifndef DHT_HPP
#define DHT_HPP
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
namespace graphlab {

template <typename KeyType, typename ValueType>
class dht { 
 public:
  typedef boost::unordered_map<size_t, ValueType> storage_type;
  
  dht(distributed_control &dc): rpc(dc, this) { }
  
  std::pair<bool, ValueType> get(const KeyType &key) const {
    // who owns the data?
    size_t hashvalue = hasher(key);
    size_t owningmachine = hashvalue % rpc.dc().numprocs();

    std::pair<bool, ValueType> retval;
    // if it is me, we can return it
    if (owningmachine == rpc.dc().procid()) {
      lock.lock();
      typename storage_type::const_iterator iter = storage.find(hashvalue);
      retval.first = iter != storage.end();
      if (retval.first) retval.second = iter->second;
      lock.unlock();
    }
    else {
      retval = rpc.fast_remote_request(owningmachine, 
                                       &dht<KeyType,ValueType>::get, 
                                       key);
    }
    return retval;
  }
  
  
  void set(const KeyType &key, const ValueType &newval) {  
    // who owns the data?
    size_t hashvalue = hasher(key);
    size_t owningmachine = hashvalue % rpc.dc().numprocs();
    // if it is me, set it
    if (owningmachine == rpc.dc().procid()) {
      lock.lock();
      storage[hashvalue] = newval;
      lock.unlock();
    }
    else {
      rpc.fast_remote_call(owningmachine, 
                           &dht<KeyType,ValueType>::set, 
                           key, newval);
    }
  }
  
 private:
  mutable dc_dist_object<dht<KeyType, ValueType> > rpc;
  
  boost::hash<KeyType> hasher;
  mutex lock;
  storage_type storage;


};

};
#endif
