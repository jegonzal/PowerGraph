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


// standard C++ headers
#include <iostream>

// includes the entire graphlab framework
#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

graphlab::glshared<size_t> testint;

void increment(graphlab::any &a, const graphlab::any& b) {
  a.as<size_t>() += b.as<size_t>();
}

void testfn_single() {
  graphlab::any b = size_t(1);
  for (size_t i = 0;i < 100; ++i) {
    testint.apply(increment, b);
    boost::shared_ptr<const size_t> ptr = testint.get_ptr();
    ASSERT_EQ((*ptr), i+1);
  }
}

void testfn() {
  graphlab::any b = size_t(1);
  for (size_t i = 0;i < 100; ++i) {
    testint.apply(increment, b);
    boost::shared_ptr<const size_t> ptr = testint.get_ptr();
    usleep(10);
    ASSERT_GE((*ptr), i);
  }
}



class GLSharedTestSuite: public CxxTest::TestSuite {
public:

  void test_glshared_single(void) {
    global_logger().set_log_level(LOG_WARNING);
    global_logger().set_log_to_console(true);
    testint.set(0);
    ASSERT_EQ(testint.get_val(), 0);
    ASSERT_EQ(*(testint.get_ptr()), 0);
    testfn_single();
    ASSERT_EQ(testint.get_val(), 100);
  }
  
  void test_glshared_multi(void) {
    global_logger().set_log_level(LOG_WARNING);
    global_logger().set_log_to_console(true);
    testint.set(0);
    ASSERT_EQ(testint.get_val(), 0);
    ASSERT_EQ(*(testint.get_ptr()), 0);
    graphlab::thread threads[4];
    for (size_t i = 0; i < 4; ++i) {
      threads[i] = graphlab::launch_in_new_thread(testfn);
    }
    for (size_t i = 0; i < 4; ++i) {
      threads[i].join();
    }

    ASSERT_EQ(testint.get_val(), 400);
  }
};

