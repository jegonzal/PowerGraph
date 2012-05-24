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




#include <cxxtest/TestSuite.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <graphlab.hpp>

void add_numbers(graphlab::atomic_add_vector<size_t>* vec_ptr, size_t count) {
  for(size_t i = 0; i < count; ++i) {
    for(size_t j = 0; j < vec_ptr->size(); ++j) {
      vec_ptr.add(j, 1);
    }
  }
}


class atomic_add_vector_tests : public CxxTest::TestSuite {
  graphlab::atomic_add_vector<size_t> vec;
  graphlab::thread_pool threads;

public:

  void test_many_adds() {
    const size_t num_threads = 32;
    const size_t count = 10000;
    const size_t vec_size = 10;
    threads.resize(num_threads);
    vec.resize(vec_size);
    for(size_t i = 0; i < threads.size(); ++i) {
      threads.launch(boost::bind(add_numbers, &vec, count));
    }
    threads.join();
    const size_t true_value = num_threads * count;
    for(size_t i = 0; i < vec.size(); ++i) {
      size_t value(-1);
      const bool success = test_and_get(i, value);
      TS_ASSERT_TRUE(success);
      TS_ASSERT_EQUALS(value, true_value);
    }
  }

}; // end of test suite
