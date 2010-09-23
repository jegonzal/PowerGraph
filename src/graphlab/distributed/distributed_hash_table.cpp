#include <vector>
#include <boost/iostreams/stream.hpp>
#include <graphlab/distributed/distributed_hash_table.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/logger/logger.hpp>
#include <graphlab/util/timer.hpp>
namespace graphlab {

static std::vector<distributed_hash_table*> allmaps;


/// Response to a get request
void int_any_map_get_reply_handler(distributed_control& dc, 
                                          procid_t source,  
                                          void* ptr,    //serialized any
                                          size_t len,   
                                          handlerarg_t mapid,
                                          handlerarg_t key,
                                          handlerarg_t reqid){
  distributed_hash_table::mapreplydata *mapreply = 
                  reinterpret_cast<distributed_hash_table::mapreplydata*>(reqid);
  mapreply->reply_found = (ptr != NULL);

  if (len > 0) {
    boost::iostreams::stream<boost::iostreams::array_source>  streamin((char*)(ptr), len);
    iarchive iarc(streamin);
    iarc >> mapreply->reply_value;
  }
  mapreply->reply_set = true;
}


/// Response to a get request
void cache_write_handler(distributed_control& dc,
                        procid_t source,
                        void* ptr,    //serialized any
                        size_t len,
                        handlerarg_t mapid,
                        handlerarg_t key){
  ASSERT_GE(len, 0);
  any value;
  boost::iostreams::stream<boost::iostreams::array_source>  streamin((char*)(ptr), len);
  iarchive iarc(streamin);
  iarc >> value;
  allmaps[mapid]->update_cache(key, value);
}



/// Set handler. Sets the local table to contain a key->value pair
void int_any_map_set_handler(distributed_control& dc, 
                             procid_t source,  
                             void* ptr,    //serialized any
                             size_t len,   
                             handlerarg_t mapid,
                             handlerarg_t key,
                             handlerarg_t exchangereqid,
                             handlerarg_t setreplyreqid) {
  boost::iostreams::stream<boost::iostreams::array_source>  streamin((char*)(ptr), len);
  iarchive iarc(streamin);
  // if this is an exchange, serialize and send the reply.
  if (exchangereqid) {
    
    distributed_hash_table::map_type::datapointer iter = allmaps[mapid]->data.find(key);
    
    // issue a get_reply
    if (iter.first == false) {
      // nope! we do not have this key. Create it, unlock and return 0
      // this is a new key. We do not need to trigger any cache updates
      // since no one will have this key anyway.
      std::pair<rwlock, any> dataelem;
      iarc >> dataelem.second;
      allmaps[mapid]->data.insert(key, dataelem);

      dc.remote_call(source,
                    int_any_map_get_reply_handler,
                    NULL,
                    0,
                    mapid, key , exchangereqid);
    }
    else {
      // nope! we have this key. Switch to the fine grain lock.
      std::pair<rwlock, any> &i = *(iter.second);
      i.first.writelock();

      // save it
      std::stringstream str;
      oarchive oarc(str);
      oarc << i.second;
      iarc >> i.second;
      i.first.unlock();

      // trigger cache updates to everyone except the source.
      // assume that the source should know to either do a cache invalidate
      // or update itself
      if (allmaps[mapid]->pushed_updates) {
        for (procid_t i = 0;i < dc.numprocs(); ++i) {
          if (i != source) {
            dc.remote_call(i, cache_write_handler,
                          ptr,len, mapid, key);
          }
        }
      }
      // issue the get reply to the source
      dc.remote_call(source, int_any_map_get_reply_handler, 
                    (void*)(str.str().c_str()), str.str().length(),mapid, key,
                    exchangereqid);
    }   
  }
  else {   
    distributed_hash_table::map_type::datapointer iter = allmaps[mapid]->data.find(key);
    // issue a get_reply
    if (iter.first == false) {
      // nope! we do not have this key. Create it, unlock and return 0  
      std::pair<rwlock, any> dataelem;
      iarc >> dataelem.second;
      allmaps[mapid]->data.insert(key, dataelem);
    }
    else {
      std::pair<rwlock, any> &i = *(iter.second);
      // nope! we have this key. Switch to the fine grain lock.
      i.first.writelock();
      iarc >> i.second;
      i.first.unlock();
      // trigger cache updates to everyone except the source.
      // assume that the source should know to either do a cache invalidate
      // or update itself
      if (allmaps[mapid]->pushed_updates) {
        for (procid_t i = 0;i < dc.numprocs(); ++i) {
          if (i != source) {
            dc.remote_call(i, cache_write_handler,
                          ptr,len, mapid, key);
          }
        }
      }

    }   
    dc.remote_call(source, set_ptr_to_value_1_handler, NULL, 0,setreplyreqid);
    
  }
  // then set it
  
}



/// Get handler. Remote request for the value to a key.
void int_any_map_get_handler(distributed_control& dc, 
                                          procid_t source,  
                                          void* ptr,    //unused
                                          size_t len,   //unused
                                          handlerarg_t mapid,
                                          handlerarg_t key,
                                          handlerarg_t reqid) {
  // see if we have it
  distributed_hash_table::map_type::datapointer iter = allmaps[mapid]->data.find(key);

  if (iter.first == false) {
    // nope! we do not have this key. Return 0
    dc.remote_call(source, int_any_map_get_reply_handler, NULL, 0,mapid, key , reqid);
  }
  else {
    std::pair<rwlock, any> &i = *(iter.second);
    std::stringstream str;
    oarchive oarc(str);
    i.first.readlock();
    oarc << i.second; 
    i.first.unlock();
    dc.remote_call(source, int_any_map_get_reply_handler, 
                   (void*)(str.str().c_str()), str.str().length(),mapid, key , reqid);
  }
}

distributed_hash_table::distributed_hash_table(distributed_control &dc,
                                   size_t max_cache_size):dc(dc),data(11) {
  numprocs = dc.numprocs();
  pushed_updates = true;
  allmaps.push_back(this);

  mapid = allmaps.size() - 1;
  // create the cache. Make sure it has room so we don't need to resize
  // too much
  cache.rehash(max_cache_size);
  maxcache = max_cache_size;
  //if (dc.procid() == 0) {
    logger(LOG_INFO, "%d Creating distributed_hash_table %d. Cache Limit = %d", dc.procid(), mapid, maxcache);
  //}
  reqs = 0;
  misses = 0;
}


void distributed_hash_table::set(size_t key, const any &value) {
  // pack the key/value pair
  std::stringstream str;
  oarchive oarc(str);
  oarc << value; 
  // figure out which node to store it
  procid_t nodeid = key_node_hash(key);
  // Run the request
  volatile size_t trigger = 0;
  dc.remote_call(nodeid, int_any_map_set_handler, 
                 (void*)(str.str().c_str()), str.str().length(), mapid, key, 0,
                 reinterpret_cast<handlerarg_t>(&trigger));
  while(trigger == 0) {
    sched_yield();
  }
  // update the cache entry locally
  update_cache(key,value);
}

void distributed_hash_table::exchange(size_t key, 
                                      const any &newvalue, 
                                      any &oldvalue) {
  // reply goes here
  mapreplydata mapreply;
  mapreply.reply_set = false;

  
  // pack the key/value pair
  std::stringstream str;
  oarchive oarc(str);
  oarc << newvalue; 
  // figure out which node to store it
  procid_t nodeid = key_node_hash(key);
  // Run the request
  dc.remote_call(nodeid, int_any_map_set_handler, 
                 (void*)(str.str().c_str()), str.str().length(), mapid, key, (size_t)(&mapreply), 0);
  // update the cache entry locally
  update_cache(key,newvalue);
  // wait for the reply
  volatile bool* b = &(mapreply.reply_set);
  while(!(*b)) {
    sched_yield();
  }
  oldvalue.swap(mapreply.reply_value);  
}

bool distributed_hash_table::get(size_t key, any &value) {
  mapreplydata mapreply;
  mapreply.reply_set = false;
  // figure out which node is supposed to own it
  procid_t nodeid = key_node_hash(key);
  // send the request
  dc.remote_call(nodeid, int_any_map_get_handler, NULL, 0, mapid, key, (size_t)(&mapreply));
  volatile bool* b = &(mapreply.reply_set);
  while(!(*b)) {
    sched_yield();
  }

  if (mapreply.reply_found) {
    // we got a value
    // copy out and release the lock
    value.swap(mapreply.reply_value);
    // set the cache
    update_cache(key,value);
    return true;
  }
  else {
    // we did not get a value. Invalidate cache and return
    invalidate(key);
    return false;
  }
}

bool distributed_hash_table::get_cached(size_t key, any &value) {
  reqs++;
  cachelock.lock();
  // check if it is in the cache
  cache_type::iterator i = cache.find(key);
  if (i == cache.end()) {
    // nope. not in cache. Call the regular get
    cachelock.unlock();
    misses++;
    return get(key, value);
  }
  else {
    // yup. in cache. return the value
    value = i->second->value;
    cachelock.unlock();
    return true;
  }
}
void distributed_hash_table::invalidate(size_t key) {
  cachelock.lock();
  // is the key I am invalidating in the cache?
  cache_type::iterator i = cache.find(key);
  if (i != cache.end()) {
    // drop it from the lru list
    delete i->second;
    cache.erase(i);
  }
  cachelock.unlock();
}

void distributed_hash_table::update_cache(size_t key, const any& val) {
  cachelock.lock();
  cache_type::iterator i = cache.find(key);
  // create a new entry
  if (i == cache.end()) {
    cachelock.unlock();
    // if we are out of room, remove the lru entry
    if (cache.size() >= maxcache) remove_lru();
    cachelock.lock();
    // insert the element, remember the iterator so we can push it
    // straight to the LRU list
    std::pair<cache_type::iterator, bool> ret = cache.insert(std::make_pair(key, new any_lru_list(key, val)));
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
}

void distributed_hash_table::remove_lru() {
  cachelock.lock();
  size_t keytoerase = lruage.back().key;
  // is the key I am invalidating in the cache?
  cache_type::iterator i = cache.find(keytoerase);
  if (i != cache.end()) {
    // drop it from the lru list
    delete i->second;
    cache.erase(i);
  }
  cachelock.unlock();
}

double distributed_hash_table::cache_miss_rate() {
  if (reqs == 0.0) return 0;
  return double(misses) / double(reqs);
}



distributed_hash_table::~distributed_hash_table() {
  data.clear();
  cache_type::iterator i = cache.begin();
  while (i != cache.end()) {
    delete i->second;
    ++i;
  }
  cache.clear();
}
}
