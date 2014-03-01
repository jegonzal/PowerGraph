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


#include <iostream>
#include <graphlab/rpc/dc.hpp>
using namespace graphlab;


void print(int val) {
  std::cout << val << std::endl;
}

int add_one(int val) {
  return val + 1;
}


int main(int argc, char ** argv) {
  // init MPI
  global_logger().set_log_level(LOG_INFO);
  mpi_tools::init(argc, argv);
  distributed_control dc;
  
  if (dc.numprocs() != 2) {
    std::cout<< "RPC Example 1: Basic Synchronous RPC\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
  
  if (dc.procid() == 0) {
    dc.control_call(1, print, 10);
    std::cout << "5 plus 1 is : " << dc.remote_request(1, add_one, 5) << std::endl;
    std::cout << "11 plus 1 is : " << dc.remote_request(1, add_one, 11) << std::endl;
  }
  dc.barrier();
  // terminate MPI
  mpi_tools::finalize();
}

