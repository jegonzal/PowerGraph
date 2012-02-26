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
 */

#include <vector>
#include <iostream>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
using namespace graphlab;


struct test_struct {
  dc_dist_object<test_struct> rmi;
  test_struct(distributed_control &dc):rmi(dc, this) { 
    dc.barrier();
  }
  
  void print(int val) {
    std::cout << rmi.procid() << ": Receiving print with value : " << val << std::endl;
  }

  void test() {
    if (rmi.procid() == 0) {
      std::cout << "First set of calls... Proc 1 and 3 should receive" << std::endl;
      std::vector<procid_t> s;
      s.push_back(1); s.push_back(3);
      rmi.remote_call(s.begin(), s.end(), &test_struct::print, 1);
    }
    rmi.full_barrier();
    
    if (rmi.procid() == 0) {
      std::cout << "Second set of calls... Proc 0 and 2 should receive" << std::endl;
      std::vector<procid_t> s;
      s.push_back(2); s.push_back(0);
      rmi.remote_call(s.begin(), s.end(), &test_struct::print, 1);
    }
    rmi.full_barrier();
  }
};

int main(int argc, char ** argv) {
  // init MPI
  mpi_tools::init(argc, argv);
  
  if (mpi_tools::size() != 4) {
    std::cout<< "RPC Example 8: Basic Broadcast Test\n";
    std::cout << "Run with exactly 4 MPI nodes.\n";
    return 0;
  }
  // set up parameters
  dc_init_param param;
  ASSERT_TRUE(init_param_from_mpi(param));
  // not necessary. Set log level. Change to LOG_WARNING to get less output
  global_logger().set_log_level(LOG_INFO);
  // create distributed control
  distributed_control dc(param);
  test_struct ts(dc);
  ts.test();

  // terminate MPI
  mpi_tools::finalize();
}
