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

#ifndef DHT_HPP
#define DHT_HPP
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
namespace graphlab {

/**
 * \ingroup rpc
 * Implements a very rudimentary distributed key value store.
  */
template <typename KeyType, typename ValueType>
class dht { 
 public:
  typedef boost::unordered_map<size_t, ValueType> storage_type;
  
  dht(distributed_control &dc): rpc(dc, this) { }
  
  /**
   * gets the value associated with a key.
   * Returns (true, Value) if the entry is available.
   * Returns (false, undefined) otherwise.
   */
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
  
  /**
   * Sets the newval to be the value associated with the key
   */
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
  
  void print_stats() const {
    std::cerr << rpc.calls_sent() << " calls sent\n";
    std::cerr << rpc.calls_received() << " calls received\n";
  }
 private:
  mutable dc_dist_object<dht<KeyType, ValueType> > rpc;
  
  boost::hash<KeyType> hasher;
  mutex lock;
  storage_type storage;


};

};
#endif
