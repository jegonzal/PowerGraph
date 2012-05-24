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
      vec_ptr->add(j, 1);
    }
  }
}


void add_and_get_numbers(graphlab::atomic_add_vector<size_t>* vec_ptr,
                         graphlab::atomic<size_t>* shared_count_ptr,
                         size_t count) {
  size_t local_count = 0;
  for(size_t i = 0; i < count; ++i) {
    for(size_t j = 0; j < vec_ptr->size(); ++j) {
      vec_ptr->add(j, 1);
      if(i % 5 == 0) {
        size_t value;
        if(vec_ptr->test_and_get(j, value)) local_count += value;
      }
    }
  }
  for(size_t j = 0; j < vec_ptr->size(); ++j) {
    size_t value;
    if(vec_ptr->test_and_get(j, value)) local_count += value;
  }
  shared_count_ptr->inc(local_count);
} // end of add and get numbers


class atomic_add_vector_tests : public CxxTest::TestSuite {
public:
  void test_many_adds() {
    graphlab::atomic_add_vector<size_t> vec;
    graphlab::thread_pool threads;

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
      const bool success = vec.test_and_get(i, value);
      TS_ASSERT(success);
      TS_ASSERT_EQUALS(value, true_value);
    }
    for(size_t i = 0; i < vec.size(); ++i) {
      size_t value(-1);
      const bool success = vec.test_and_get(i, value);
      TS_ASSERT(!success);
      TS_ASSERT(vec.empty(i));
    }
  }

  void test_many_adds_and_gets() {
    
    graphlab::atomic_add_vector<size_t> vec;
    graphlab::atomic<size_t> shared_counter; 
    graphlab::thread_pool threads;
    const size_t num_threads = 32;
    const size_t count = 10000;
    const size_t vec_size = 10;
    threads.resize(num_threads);
    vec.resize(vec_size);
    for(size_t i = 0; i < threads.size(); ++i) {
      threads.launch(boost::bind(add_and_get_numbers, &vec, 
                                 &shared_counter, count));
    }
    threads.join();
    const size_t true_value = threads.size() * count * vec.size();
    TS_ASSERT_EQUALS(shared_counter.value, true_value);
  }
}; // end of test suite
