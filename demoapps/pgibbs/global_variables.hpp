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
