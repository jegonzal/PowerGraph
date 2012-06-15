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
#include <sstream>
#include <vector>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
using namespace graphlab;


void print(std::string val) {
  std::cout << val << std::endl;
}

std::vector<int> add_one(std::vector<int> val) {
  val.push_back(1);
  return val;
}


int main(int argc, char ** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  distributed_control dc;
  
  if (dc.numprocs() != 2) {
    std::cout<< "RPC Example 2: Asynchronous RPC with Built-in Serialization\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
  
  dc.barrier();
  if (dc.procid() == 0) {
    dc.remote_call(1, print, "hello world!");
    // Create a vector with a few elements
    std::vector<int> vec;
    vec.push_back(1); vec.push_back(2);
    // Call the remote machine 
    vec = dc.remote_request(1, add_one, vec);
    
    std::stringstream strm;
    // Print the vector 
    for (size_t i = 0; i < vec.size(); ++i) {
      strm << vec[i] << ", ";
    }
    strm << std::endl;
    strm.flush();
    dc.remote_call(1, print, strm.str());
  }
  dc.barrier();

  mpi_tools::finalize();
}

