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

#include "clustering.h"
#include <graphlab/macros_def.hpp>

extern advanced_config ac;
extern problem_setup ps;

extern const char * testtypename[];

// calc statistics about matrix/tensor and exit  
void calc_stats(testtype type){
   graph_type * gr = ps.g<graph_type>(type);

  double avgval=0, minval=1e100, maxval=-1e100;
  int nz=0;
  int reported =0;
  int col_reported = 0; 

  for (int i=0; i< (type == TRAINING ? ps.M: ps.M_validation); i++){ 
    vertex_data * data = &gr->vertex_data(i);
      if (min(data->datapoint) < minval)
	 minval = min(data->datapoint);
      if (max(data->datapoint) > maxval)
	 maxval = max(data->datapoint);
      if (data->reported)
	 reported++;

      nz += nnz(data->datapoint);
      avgval += sum(data->datapoint);
  }
 
 for (uint i=(type ==TRAINING? ps.M: ps.M_validation); i < gr->num_vertices(); i++){
    vertex_data * data = &gr->vertex_data(i);
    if (data->reported)
        col_reported++;

  }
 
 avgval /= (double)nz;
 printf("%s Matrix size is %d rows %d cols %d nnz\n", testtypename[type], ps.M, ps.N, ps.L);
 printf("Avg matrix value %g min val %g max value %g\n", avgval, minval, maxval); 
 printf("Rows with non-zero entries: %d\n", reported);
 if (col_reported > 0)
   printf("Cols with non-zero entries: %d\n", col_reported);
 assert(nz == ps.L);
}


#include <graphlab/macros_undef.hpp>
#endif //_STATS_HPP
