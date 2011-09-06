#include <pthread.h>
#include <map>
#include <graphlab/util/cache.hpp>
#include <graphlab/shared_data/sharedsum.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {
  namespace sharedsum_impl {


    typedef cache::lru<isharedsum*, icache_entry*> cache_type;


    size_t max_cache_size = 1024;
    
    size_t& cache_size() { return max_cache_size; }


    void destroy_tls_data(void* ptr) {
      cache_type* cache_ptr = static_cast<cache_type*>(ptr);
      if(cache_ptr != NULL) { 
        typedef cache_type::cache_map_type::iterator iterator;
        cache_type& cache = *cache_ptr;
        for(iterator iter = cache.begin(), iend = cache.end(); iter != iend; ++iter) {
          if(iter->get_right() != NULL) {
            delete iter->get_right();
            iter->get_right() = NULL;
          }
        }
        delete cache_ptr;
      }
    }
    struct tls_key_creator {
      pthread_key_t TLS_KEY;
      tls_key_creator() : TLS_KEY(0) {
        pthread_key_create(&TLS_KEY, destroy_tls_data);
      }
    }; 
    const tls_key_creator key;     
    cache_type& get_cache() {
      cache_type* cache_ptr = static_cast<cache_type*> 
        (pthread_getspecific(key.TLS_KEY));
      if(cache_ptr == NULL) {
        cache_ptr = new cache_type();
        pthread_setspecific(key.TLS_KEY, cache_ptr);
      }
      ASSERT_NE(cache_ptr, NULL);
      return *cache_ptr;
    }



    icache_entry* get_cache_entry(isharedsum* key) {
      icache_entry* result = NULL;
      const bool successful = get_cache().get(key, result);
      return successful? result : NULL;
    }    


    void add_cache_entry(isharedsum* id, icache_entry* ptr) {
      while(size() > max_cache_size) evict();
      cache_type& cache = get_cache();
      cache[id] = ptr;
    }

    bool evict(isharedsum* key) {
      icache_entry* entry = NULL;
      const bool success = get_cache().evict(key, entry);
      if(success) {
        ASSERT_NE(entry, NULL);
        ASSERT_NE(key, NULL);
        key->flush(entry);
      }
      return success;
    }


    void evict() {
      typedef std::pair<isharedsum*, icache_entry*> pair_type;
      pair_type evicted_pair = get_cache().evict();
      ASSERT_NE(evicted_pair.first, NULL);
      ASSERT_NE(evicted_pair.second, NULL);
      evicted_pair.first->flush(evicted_pair.second);
    }


    size_t size() { return get_cache().size(); }



  }; // end of sharedsum_impl;




}; // end of graphlab namespace


