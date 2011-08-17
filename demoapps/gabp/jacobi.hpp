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


/**
 * Functionality: The code solves the linear system Ax = b using
 * The Jacobi algorithm. (A is a square matrix). 
 * A assumed to be full column rank.  Algorithm is described
 * http://en.wikipedia.org/wiki/Jacobi_method
 * Written by Danny Bickson 
*/

#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <cmath>
#include <cstdio>
#include "linear.h"
#include <graphlab/macros_def.hpp>


/***
 * JACOBI UPDATE FUNCTION
 * x_i = (b_i - \sum_j A_ij * x_j)/A_ii
 */
void jacobi_update_function(gl_types::iscope &scope,
                            gl_types::icallback &scheduler) {
  

  /* GET current vertex data */
  vertex_data& vdata = scope.vertex_data();
  gl_types::edge_list outedgeid = scope.out_edge_ids();

  const bool& debug = DEBUG_KEY.get_val();

  //store last round values
  vdata.prev_mean = vdata.cur_mean;

  //initialize accumlated values in x_i
  sdouble x_i = vdata.prior_mean;
  sdouble A_ii = vdata.prior_prec;
  assert(A_ii != 0);

  if (debug) 
    std::cout << "entering node " << scope.vertex() << " P=" << vdata.prior_prec << " u=" << vdata.prior_mean << std::endl;
  
    for(size_t i = 0; i < outedgeid.size(); ++i) {
      edge_data& out_edge = scope.edge_data(outedgeid[i]);
      const vertex_data & other = scope.neighbor_vertex_data(scope.target(outedgeid[i]));
      x_i -= out_edge.weight * other.cur_mean;
    }

    x_i /= A_ii;
    vdata.cur_mean = x_i;

    if (debug)
       std::cout<<scope.vertex()<<") x_i: "<<x_i<<std::endl;
}

#include <graphlab/macros_undef.hpp>
#endif
