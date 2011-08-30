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
 *  Code written by Danny Bickson, CMU
 *  Any changes to the code must include this original license notice in full.
 */


#ifndef _STATS_HPP
#define _STATS_HPP

#include <graphlab/macros_def.hpp>

extern advanced_config ac;
extern problem_setup ps;

double min(sparse_vec & dvec);
double max(sparse_vec & dvec);
double sum(sparse_vec & dvec);

// calc statistics about matrix/tensor and exit  
void calc_stats(){
   graph_type * gr = ps.g;

  double avgval=0, minval=1e100, maxval=-1e100;
  int nz=0;
 
  for (int i=0; i< ps.M; i++){ 
    vertex_data * data = &gr->vertex_data(i);
      if (min(data->datapoint) < minval)
	 minval = min(data->datapoint);
      if (max(data->datapoint) > maxval)
	 maxval = max(data->datapoint);
      
      nz += ((sparse_vec)data->datapoint).nnz();
      avgval += sum(data->datapoint);
  }
 
 avgval /= (double)nz;
 printf("Avg matrix value %g min val %g max value %g\n", avgval, minval, maxval);
 assert(nz == ps.L);
}


#include <graphlab/macros_undef.hpp>
#endif //_STATS_HPP
