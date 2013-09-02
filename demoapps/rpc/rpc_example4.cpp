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
#include <map>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/serialization/serialization_includes.hpp>

using namespace graphlab;

struct teststruct{
  int a;
  double b;
  std::vector<std::string> name;
  
  void save(oarchive &arc) const {
    arc << a << b << name;
  };
  void load(iarchive &arc) {
    arc >> a >> b >> name;
  };
};


void print(std::map<int, teststruct> &data1,  
           std::string data2) {
  std::cout << "1.a = " << data1[1].a << std::endl;
  std::cout << "2.b = " << data1[2].b << std::endl;
  std::cout << "string = " << data2 << std::endl;
}



int main(int argc, char ** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;

  if (dc.numprocs() != 2) {
    std::cout<< "RPC Example 4: Asynchronous RPC with Manual Serialization\n";
    std::cout << "Run with exactly 2 MPI nodes.\n";
    return 0;
  }
    
  if (dc.procid() == 0) {
    std::map<int, teststruct> data;
    data[1].a = 10;
    data[2].b = 15.0;
    dc.remote_call(1, print, data, "hello world!");
  }
  dc.barrier();
  
  mpi_tools::finalize();
}
