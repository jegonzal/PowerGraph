/**  
 * Copyright (c) 2009 Carnegie Mellon University. 
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


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
