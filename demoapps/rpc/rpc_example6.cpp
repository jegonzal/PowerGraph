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
#include <string>
#include <vector>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/generics/any.hpp>
using namespace graphlab;


void print(any val) {
  val.print(std::cout);
  std::cout << std::endl;
}


int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;

  if (dc.numprocs() != 2) {
    std::cout<< "RPC Example 6: Asynchronous RPC with any \n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }

  if (dc.procid() == 0) {
    dc.remote_call(1, print, any(15));
    dc.remote_call(1, print, any(10.5));
    dc.remote_call(1, print, any(std::string("hello world")));    
  }
    
  int i = dc.procid() == 0 ? 10 : 100;
  dc.broadcast(i, dc.procid() == 0);
  std::cout << i << std::endl;
  assert(i == 10);
  
  mpi_tools::finalize();
}
