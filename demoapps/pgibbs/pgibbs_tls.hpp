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


// //! Key used to get the pgibbs tls
// extern pthread_key_t pgibbs_tls_key;

//! Local state available to each thread
struct pgibbs_tls {
  factor_t cavity;
  factor_t conditional_factor;
  factor_t belief;
  factor_t tmp_belief;
};


pgibbs_tls& get_pgibbs_tls();




#endif
