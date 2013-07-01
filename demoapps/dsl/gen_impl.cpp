/*  
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


#include <cmath>
#include <iostream>
#include "graph_typedefs.gen"

using namespace graphlab;

#define ALPHA 0.87

extern "C" void user_program(user_funs* f) {
    float prev = f->get_vertex_data();
    float neighbors = f->reduce_neighbors(IN_EDGES);
    //neighbors = neighbors/out_edges
    float curr = ALPHA*neighbors + (1-ALPHA);
    f->set_vertex_data(curr);
    float last_change = std::abs(curr - prev);
    std::cout << "last change was: " << last_change << std::endl;
    if (last_change > 0.01) {
	f->signal_neighbors(OUT_EDGES);
    }
    std::cout.flush();
}

extern "C" void vertex_reduce(vertex_data_type& a, const vertex_data_type& b) {
    a += b;
}

