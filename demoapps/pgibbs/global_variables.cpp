#include "global_variables.hpp"

// Global Shared Varaibles ====================================================>
//graphlab::glshared_const<factorized_model::factor_map_t> SHARED_FACTORS;
const factorized_model::factor_map_t* SHARED_FACTORS_PTR = NULL;
graphlab::glshared_const<size_t> MAX_NSAMPLES;
graphlab::glshared<size_t> n_samples;
