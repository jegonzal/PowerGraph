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
#include <graphlab/util/timer.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/parallel/fiber_group.hpp>
#include <graphlab/parallel/fiber_remote_request.hpp>
using namespace graphlab;

atomic<size_t> complete_count;

size_t some_remote_function(size_t a) {
  return a;
}

void test_fiber(size_t sequential_count) {
  for (size_t i = 0;i < sequential_count; ++i) {
    request_future<size_t> ret = fiber_remote_request(1, some_remote_function, 1);
    complete_count.inc(ret()); 
  }
}

int main(int argc, char** argv) {
  mpi_tools::init(argc, argv);
  distributed_control dc;
  timer ti;
  // with fibers
  if (dc.procid() == 0) {
    fiber_group group(8192);
    for (int i = 0;i < 1600000; ++i) {
      group.launch(boost::bind(test_fiber, 1));
      if (i % 100000 == 0) std::cout << i << "\n";
    }
    group.join();
    std::cout << "completed requests: " << complete_count.value << " in " << ti.current_time() << "\n";  
  }

  dc.barrier();
  mpi_tools::finalize();
}

