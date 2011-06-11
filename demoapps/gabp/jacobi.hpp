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
  graphlab::edge_list outedgeid = scope.out_edge_ids();

  const bool& debug = DEBUG_KEY.get();

  //store last round values
  vdata.prev_mean = vdata.cur_mean;

  //initialize accumlated values
  sdouble x_i = vdata.prior_mean;
  sdouble A_ii = vdata.prior_prec;
  assert(A_ii != 0);

  if (debug) 
    std::cout << "entering node " << scope.vertex() << " P=" << vdata.prior_prec << " u=" << vdata.prior_mean << std::endl;
  
  /* SEND new value and schedule neighbors */
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
