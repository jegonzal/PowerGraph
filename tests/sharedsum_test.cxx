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



// Test the graph class

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <cxxtest/TestSuite.h>


#include <graphlab.hpp>
#include <graphlab/shared_data/sharedsum.hpp>

#include <graphlab/macros_def.hpp>



class GraphTestSuite: public CxxTest::TestSuite {
public:
  

  
    
  // construct the graph at startup 
  GraphTestSuite() {
    global_logger().set_log_level(LOG_INFO);
    global_logger().set_log_to_console(true);
  }
 
  void test_sequential() {
    graphlab::sharedsum<int> x;
    for(size_t i = 0; i < 100; ++i) {
      x += 1;
    }
    std::cout << std::endl;
    std::cout << x.get() << std::endl;
  } // end of test_sequential 

  
  struct functor {
    graphlab::sharedsum<int>& x;
    size_t min, max;
    functor(graphlab::sharedsum<int>& x,
            size_t min, size_t max) :
      x(x), min(min), max(max) { };
    void operator()() {
      for(size_t i = min; i < max; ++i) {
        x += 1;
      }
      std::cout << "value: " << x << std::endl;
      x.flush();
    }
  }; // end of functor

  void test_parallel() {
    const size_t total = 100;
    const size_t step_size = 10000;
    graphlab::sharedsum<int> x(0, 1000);
    graphlab::thread_group threads;
    for(size_t i = 0; (i+1) < total * step_size; i += step_size) {
      threads.launch(functor(x, i, i + step_size)); 
    }
    threads.join();
    std::cout << std::endl;
    std::cout << x.get() << std::endl;

  } // end of test parallel

};



#include <graphlab/macros_undef.hpp>
