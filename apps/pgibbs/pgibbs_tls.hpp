#ifndef PGIBBS_TLS_HPP
#define PGIBBS_TLS_HPP



/**
 *
 * This code is used to represent thread local storage needed in some
 * of the sampler code
 *
 *  \author Joseph Gonzalez
 */

// INCLUDES ===================================================================>

#include <pthread.h>


#include "factorized_model.hpp"


//! Key used to get the pgibbs tls
pthread_key_t pgibbs_tls_key;

//! Local state available to each thread
struct pgibbs_tls {
  factor_t cavity;
  factor_t conditional_factor;
  factor_t belief;
  factor_t tmp_belief;
};




pgibbs_tls* create_pgibbs_tls() {
  assert(pthread_getspecific(pgibbs_tls_key) == NULL);
  pgibbs_tls* data = new pgibbs_tls();
  assert(data != NULL);
  pthread_setspecific(pgibbs_tls_key, data);
  return data;
}

pgibbs_tls& get_pgibbs_tls() {
  pgibbs_tls* tls =
    reinterpret_cast<pgibbs_tls*>
    (pthread_getspecific(pgibbs_tls_key) );
  // If no tsd be has been associated, create one
  if(tls == NULL) tls = create_pgibbs_tls();
  assert(tls != NULL);
  return *tls;
}

void destroy_pgibbs_tls(void* ptr) {
  pgibbs_tls* tls = 
    reinterpret_cast<pgibbs_tls*>(ptr);
  if(tls != NULL) delete tls;

}


struct pgibbs_tls_key_creater {
  pgibbs_tls_key_creater( )  {
    pthread_key_create(&pgibbs_tls_key,
                       destroy_pgibbs_tls);
  }
};
static const pgibbs_tls_key_creater make_pgibbs_tls_key;




#endif
