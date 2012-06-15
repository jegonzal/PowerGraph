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

#include <vector>
#include <iostream>
#include <graphlab/rpc/dc.hpp>
using namespace graphlab;


void print(distributed_control &dc, procid_t caller, int val) {
  std::cout << dc.procid() << ": Receiving print with value : " << val << std::endl;
}



int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;

  if (dc.numprocs() != 4) {
    std::cout<< "RPC Example 8: Basic Broadcast Test\n";
    std::cout << "Run with exactly 4 MPI nodes.\n";
    return 0;
  }
  
  if (dc.procid() == 0) {
    std::cout << "First set of calls... Proc 1 and 3 should receive" << std::endl;
    std::vector<procid_t> s;
    s.push_back(1); s.push_back(3);
    dc.remote_call(s.begin(), s.end(), print, 1);
  }
  dc.full_barrier();
  
  if (dc.procid() == 0) {
    std::cout << "Second set of calls... Proc 0 and 2 should receive" << std::endl;
    std::vector<procid_t> s;
    s.push_back(2); s.push_back(0);
    dc.remote_call(s.begin(), s.end(), print, 1);
  }
  dc.full_barrier();
  // terminate MPI
  mpi_tools::finalize();
}
