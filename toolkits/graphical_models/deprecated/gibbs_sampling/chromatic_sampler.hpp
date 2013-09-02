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


#ifndef CHROMATIC_SAMPLER_HPP
#define CHROMATIC_SAMPLER_HPP

#include <graphlab.hpp>

#include "mrf.hpp"

class gibbs_update :
  public graphlab::iupdate_functor<mrf_graph_type, gibbs_update> {
  typedef graphlab::iupdate_functor<mrf_graph_type, gibbs_update> base;
  void operator()(base::icontext_type& context);
}; // end of class gibbs update

/** Get the update counts for a vertex */
inline size_t get_nsamples(const mrf_vertex_data& vdata) { 
  return vdata.nsamples; 
}


//! Run the chromatic sampler for a fixed ammount of time
void run_chromatic_sampler(graphlab::core<mrf_graph_type, gibbs_update>& core, 
                           const std::string& chromatic_results_fn,
                           const std::vector<double>& runtime,
                           const bool draw_images);



#endif
