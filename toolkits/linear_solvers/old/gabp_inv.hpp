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
 * GRAPHLAB implementation of Gaussiabn Belief Propagation Code See
 * algrithm description and explanation in: Danny Bickson, Gaussian
 * Belief Propagation: Theory and Application. Ph.D. Thesis. The
 * Hebrew University of Jerusalem. Submitted October 2008.
 * http://arxiv.org/abs/0811.2518 By Danny Bickson, CMU. Send any bug
 * fixes/reports to bickson@cs.cmu.edu 
 *
 * Functionality: The code solves the linear system Ax = b using
 * Gaussian Belief Propagation. (A is either square matrix or
 * skinny). A assumed to be full column rank.  Algorithm is described
 * in Algorithm 1, page 14 of the above Phd Thesis.
 *
 * If you are using this code, you should cite the above reference. Thanks!
 */

#ifndef GABP_INV_HPP
#define GABP_INV_HPP

#include <cmath>
#include <cstdio>
#include "linear.h"
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>



/***
 * UPDATE FUNCTION
 */
void gabp_update_inv_function(gl_types_inv::iscope &scope,
                          gl_types_inv::icallback &scheduler) {


  /* GET current vertex data */
  vertex_data_inv& vdata = scope.vertex_data();
  gl_types_inv::edge_list inedgeid = scope.in_edge_ids();
  gl_types_inv::edge_list outedgeid = scope.out_edge_ids();

  const bool& support_null_variance  = SUPPORT_NULL_VARIANCE_KEY.get_val();
  const bool& round_robin = ROUND_ROBIN_KEY.get_val();
  const bool& debug = DEBUG_KEY.get_val();
  const int n = MATRIX_WIDTH_KEY.get_val();
  assert(n>0);

  //store last round values
  memcpy(vdata.prev_mean, vdata.cur_mean, sizeof(sdouble)*n);
  //sdouble prev_prec = vdata.cur_prec;

  //initialize accumlated values
  sdouble * mu_i = new sdouble[n];
  memcpy(mu_i, vdata.prior_mean, n*sizeof(sdouble));

  /* CALCULATE new value */
  if (debug) {
    std::cout << "entering node " << scope.vertex()
              << " P=" << vdata.prior_prec
              << " u=" << vdata.prior_mean[0]
              << std::endl;
  }


  for (int j=0; j<n; j++){
  //accumlate all messages (the inner summation in section 4 of Algorithm 1)
  sdouble J_i = vdata.prior_prec;
  if (!support_null_variance) assert(J_i != 0);
   foreach(gl_types_inv::edge_id eid, inedgeid) {
    const edge_data_inv& edata = scope.edge_data(eid);
    mu_i[j] += edata.mean[j];
    J_i +=  edata.prec;
  }

  if (debug) {
    std::cout << scope.vertex() << ") summing up all messages in problem " << j 
              << " Mu: " << mu_i[j] << " " << J_i << std::endl;
  }

  // optional support for null variances
  if (support_null_variance && J_i == 0){
    vdata.cur_mean[j] = mu_i[j];
    vdata.cur_prec = 0;
  } else {
    assert(J_i != 0);
    vdata.cur_mean[j] = mu_i[j] / J_i;
    assert(vdata.cur_mean[j] != NAN);
    vdata.cur_prec = J_i;
  }
  assert(vdata.cur_mean[j] != NAN);

  /* SEND new value and schedule neighbors */
    for(size_t i = 0; i < inedgeid.size(); ++i) {
      assert(scope.source(inedgeid[i]) == scope.target(outedgeid[i]));
      edge_data_inv& in_edge = scope.edge_data(inedgeid[i]);
      edge_data_inv& out_edge = scope.edge_data(outedgeid[i]);
      graphlab::vertex_id_t target = scope.target(outedgeid[i]);

      //substruct the sum of message sent from node j
      sdouble mu_i_j = mu_i[j] - in_edge.mean[j];
      sdouble J_i_j  = J_i - in_edge.prec;

      if (!support_null_variance)  assert(J_i_j != 0);
      assert(out_edge.weight != 0);

      if (support_null_variance && J_i_j == 0){
        out_edge.mean[j] = 0;
        out_edge.prec = 0;
      } else {
        //compute the update rule (Section 4, Algorithm 1)
        out_edge.mean[j] = -(out_edge.weight * mu_i_j / J_i_j);
        out_edge.prec = -((out_edge.weight * out_edge.weight) / J_i_j);//matrix is assumed symmetric!
      }

      if (!round_robin) {
        gl_types_inv::update_task task(target, gabp_update_inv_function);
        double priority = fabs(vdata.cur_prec) + 1e-5;
        scheduler.add_task(task, priority);
      }

      if (debug) {
        std::cout << "Sending to " << target << " "
                  << out_edge.mean[j] << " "
                  << out_edge.prec << " wdge weight "
                  << out_edge.weight << std::endl;
      }
    }
   }

}

#include <graphlab/macros_undef.hpp>
#endif
