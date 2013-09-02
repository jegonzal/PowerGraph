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
  run_statistics(const mrf_graph_type& mrf) :
    nsamples(0), nchanges(0), loglik(0.0),
    min_samples(std::numeric_limits<size_t>::max()), max_samples(0) {
    typedef mrf_graph_type::vertex_id_type vertex_id_type;
    // Compute the unnormalized log likelihood
    loglik = unnormalized_loglikelihood(mrf);
    for(vertex_id_type vid = 0; vid < mrf.num_vertices(); ++vid) {
      const mrf_vertex_data& vdata = mrf.vertex_data(vid);
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
