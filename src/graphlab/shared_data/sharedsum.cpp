#include <pthread.h>
#include <map>

#include <graphlab/shared_data/sharedsum.hpp>

namespace graphlab {
  namespace sharedsum_impl {
    typedef std::map<void* ipartial_sum*> partial_sum_map_type;

    void destroy_tls_data(void* ptr) {
      partial_sum_map_type* ps_map_ptr = 
        static_cast<partial_sum_map_type*>(ptr);
      if(ps_map_ptr != NULL) { delete ps_map_ptr; }
    }
    struct tls_key_creator {
      pthread_key_t TLS_KEY;
      tls_key_creator() : TLS_KEY(0) {
        pthread_key_create(&TLS_KEY, destroy_tls_data);
      }
    }; 
    const tls_key_creator key;     


    partial_sum_map_type& get_partial_sum_map() {
      partial_sum_map_type* ps_map_ptr = 
        static_cast<partial_sum_map_type*> 
        (pthread_getspecific(key.TLS_KEY));
      if(ps_map_ptr == NULL) {
        ps_map_ptr = new partial_sum_map_type();
        pthread_setspecific(key.TLS_KEY, ps_map_ptr);
      }
      ASSERT_NE(ps_map_ptr, NULL);
      return *ps_map_ptr;
    }


    ipartial_sum* get_tl_partial_sum(void* sharedsum_ptr) {
      partial_sum_map_type& ps_map = get_partial_sum_map();
      return ps_map[sharedsum_ptr];
    }

    
    void set_tl_partial_sum(void* sharedsum_ptr, ipartial_sum* ps_ptr) {
      get_partial_sum_map()[sharedsum_ptr] = ps_ptr;
    }

  }; // end of sharedsum_impl;




}; // end of graphlab namespace


