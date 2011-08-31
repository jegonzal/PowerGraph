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



/**
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
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
 
  void test_sequential_sum() {
    graphlab::sharedsum<int> x;
    for(size_t i = 0; i < 100; ++i) {
      x += 1;
    }
    std::cout << std::endl;
    std::cout << x << std::endl;
    TS_ASSERT_EQUALS(x.val(), 100);
  } // end of test_sequential 

  
  struct functor {
    graphlab::sharedsum<int>& x;
    int min, max;
    functor(graphlab::sharedsum<int>& x,
            int min, int max) :
      x(x), min(min), max(max) { };
    void operator()() {
      for(int i = min; i < max; ++i) {
        x += 1;
      }
      x.flush();
    }
  }; // end of functor

  void test_parallel_sum() {
    const int total = 10;
    const int step_size = 1000;
    graphlab::sharedsum<int> x(0, 100);
    graphlab::thread_group threads;
    for(int i = 0; (i+1) < total * step_size; i += step_size) {
      threads.launch(functor(x, i, i + step_size)); 
    }
    threads.join();
    std::cout << std::endl;
    std::cout << x << std::endl;
    TS_ASSERT_EQUALS(x.val(), total*step_size); 
  } // end of test parallel



  // struct functor {
  //   graphlab::sharedsum<int>& x;
  //   int min, max;
  //   functor(graphlab::sharedsum<int>& x,
  //           int min, int max) :
  //     x(x), min(min), max(max) { };
  //   void operator()() {
  //     for(int i = min; i < max; ++i) {
  //       x += 1;
  //     }
  //     x.flush();
  //   }
  // }; // end of functor



  // void test_parallel_density() {
  //   const int total = 100;
  //   const int step_size = 10000;
  //   std::vector<graphlab::sharedsum<int>> counts(10000, 
  //                                                graphlab::shared_sum<int>(0, 100));
  //   graphlab::thread_group threads;
  //   for(int i = 0; (i+1) < total * step_size; i += step_size) {
  //     threads.launch(functor(x, i, i + step_size)); 
  //   }
  //   threads.join();
  //   std::cout << std::endl;
  //   std::cout << x << std::endl;
  //   TS_ASSERT_EQUALS(x.val(), total*step_size); 
  // } // end of test parallel




};



#include <graphlab/macros_undef.hpp>
