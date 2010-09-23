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
  
  std::pair<bool, ValueType> get_from_hashvalue(size_t hashvalue) const {
    std::pair<bool, ValueType> retval;
    // who owns the data?
    size_t owningmachine = hashvalue % rpc.dc().numprocs();
    // if it is me, we can return it
    if (owningmachine == rpc.dc().procid()) {
      lock.lock();
      typename storage_type::const_iterator iter = storage.find(hashvalue);
      retval.first = iter != storage.end();
      if (retval.first) retval.second = iter->second;
      lock.unlock();
      return retval;
    }
    else {
      retval = rpc.fast_remote_request(owningmachine, &dht<KeyType,ValueType>::get_from_hashvalue, hashvalue);
    }
    return retval;
  }
  
  std::pair<bool, ValueType> get(const KeyType &key) const {
    size_t hashvalue = hasher(key);
    return get_from_hashvalue(hashvalue);
  }
  
  void set_from_hashvalue(size_t hashvalue, const ValueType &newval)  {
    // who owns the data?
    size_t owningmachine = hashvalue % rpc.dc().numprocs();
    // if it is me, we can return it
    if (owningmachine == rpc.dc().procid()) {
      //std::cerr << "local set" << std::endl;
      lock.lock();
      storage[hashvalue] = newval;
      lock.unlock();
    }
    else {
      //std::cerr << "remote set to " << owningmachine << std::endl;
      rpc.fast_remote_call(owningmachine, &dht<KeyType,ValueType>::set_from_hashvalue, hashvalue, newval);
    }
  }
  
  void set(const KeyType &key, const ValueType &newval) {  
    size_t hashvalue = hasher(key);
    set_from_hashvalue(hashvalue, newval);
  }
  
 private:
  mutable dc_dist_object<dht<KeyType, ValueType> > rpc;
  
  boost::hash<KeyType> hasher;
  mutex lock;
  storage_type storage;


};

};
#endif
