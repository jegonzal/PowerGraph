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


#ifndef  PGIBBS_GLOBAL_VARIABLES
#define  PGIBBS_GLOBAL_VARIABLES


#include <graphlab.hpp>

#include "factorized_model.hpp"

// Global Shared Varaibles ====================================================>
//extern graphlab::glshared_const<factorized_model::factor_map_t*> SHARED_FACTORS;
extern const factorized_model::factor_map_t* SHARED_FACTORS_PTR;
extern graphlab::glshared_const<size_t> MAX_NSAMPLES;
extern graphlab::glshared<size_t> n_samples;

#endif
