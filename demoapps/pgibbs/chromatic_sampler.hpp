#ifndef CHROMATIC_SAMPLER_HPP
#define CHROMATIC_SAMPLER_HPP

#include "mrf.hpp"

void single_gibbs_update(mrf_gl::iscope& scope, 
                         mrf_gl::icallback& scheduler);


/** Get the update counts for a vertex */
inline size_t get_nsamples(const mrf_vertex_data& vdata) { 
  return vdata.nsamples; 
}


// bool nsamples_terminator(const mrf_gl::ishared_data* shared_data);


//! Run the chromatic sampler for a fixed ammount of time
void run_chromatic_sampler(mrf_gl::core& core, 
                           const std::string& chromatic_results_fn,
                           const std::vector<double>& runtime,
                           const bool draw_images);



#endif
