#ifndef PGIBBS_RUN_STATISTICS_HPP
#define PGIBBS_RUN_STATISTICS_HPP

#include "mrf.hpp"

//! Statistics associated with a run
struct run_statistics {
  size_t nsamples;
  size_t nchanges;
  double loglik;
  size_t min_samples;
  size_t max_samples;
  run_statistics() :
    nsamples(0), nchanges(0), loglik(0.0),
    min_samples(std::numeric_limits<size_t>::max()), max_samples(0) { }
  run_statistics(const mrf_gl::core& core) :
    nsamples(0), nchanges(0), loglik(0.0),
    min_samples(std::numeric_limits<size_t>::max()), max_samples(0){
    // Compute the unnormalized log likelihood
    loglik = unnormalized_loglikelihood(core);
    for(vertex_id_t vid = 0; vid < core.graph().num_vertices(); ++vid) {
      const mrf_vertex_data& vdata = core.graph().vertex_data(vid);
      nsamples += vdata.nsamples;
      nchanges += vdata.nchanges;
      min_samples = std::min(min_samples, vdata.nsamples);
      max_samples = std::max(max_samples, vdata.nsamples);
    } // end of for loop
  } // end of compute run statistics
  
  void print() const {
    std::cout << "nsamples:        " << nsamples << std::endl
              << "nchanges:        " << nchanges << std::endl
              << "loglik:          " << loglik   << std::endl
              << "min_samples:     " << min_samples << std::endl
              << "max_samples:     " << max_samples << std::endl;
  }
};

#endif
