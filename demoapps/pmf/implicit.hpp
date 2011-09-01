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

#ifndef _IMPLICIT_HPP
#define _IMPLICIT_HPP

#include <graphlab/macros_def.hpp>

extern advanced_config ac;
extern problem_setup ps;

void add_implicit_edges(graph_type * g){
#ifdef GL_NO_MULT_EDGES
   assert(ac.implicitratingpercentage>= 0 && ac.implicitratingpercentage<=1);
   assert(ac.implicitratingtype != "none");

   if (ac.implicitratingtype != "uniform" && ac.implicitratingtype != "user"){
      logstream(LOG_ERROR) << "item implicit ratings not implemented yet!" << std::endl;
      return;
   }

   edge_data data;
   data.weight = ac.implicitratingvalue;
   data.time = ac.implicitratingweight;

   int added = 0;
 
   bool *flag_edges = new bool[ps.N];
   for (int i=0; i< ps.M; i++){
      memset(flag_edges, 0, ps.N*sizeof(bool));
      foreach(gl_types::edge_id oedgeid, g->out_edge_ids(i)) {
          int to = g->target(oedgeid)-ps.M;
          assert(to >= 0 && to < ps.N);
          flag_edges[to]=true;
      }
      float toadd;
      if (ac.implicitratingtype == "uniform")
	toadd  = ac.implicitratingpercentage*ps.N;
      else if (ac.implicitratingtype == "user")
        toadd = ac.implicitratingpercentage*g->out_edge_ids(i).size();
      if (ac.implicitratingtype == "uniform" && i == 0 && toadd < 1){
         logstream(LOG_WARNING) << "implicitratingpercentage given is too low, resulting in " << toadd << " new edges per node. No edges will be added." << std::endl;
      }
      ivec newedges = randi(toadd,0,ps.N-1);
      assert(newedges.size() <= ps.N);
      for (int j=0; j< newedges.size(); j++){
	 if (!flag_edges[newedges[j]]){
		ps.g->add_edge(i,ps.M+ newedges[j], data);
 		flag_edges[newedges[j]] = true;
		added++;
         }
      } 
   }

   ps.L+=added; //update edge count including added edges
   logstream(LOG_INFO) << "added " << added << " implicit edges, rating=" <<ac.implicitratingvalue << " weight=" << ac.implicitratingweight << " type=" << ac.implicitratingtype << std::endl;
#else
   logstream(LOG_ERROR) << "implicit edges addition is not supported in multiple edges between user and movie mode. Uncomment the flag GL_NO_MULT_EDGES in pmf.h and recompile" << std::endl;
#endif
}
















#include <graphlab/macros_undef.hpp>

#endif //_IMPLICIT_HPP
