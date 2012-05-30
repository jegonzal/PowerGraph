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


#include <fstream>
#include <vector>
#include <graphlab/scheduler/scheduler_includes.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <cxxtest/TestSuite.h>


using namespace graphlab;

struct message_type {
  int value;
  double pr;

  explicit message_type(int value = 0, double pr = 0):value(value),pr(pr) { }
  
  double priority() const {
    return pr;
  }
  message_type& operator+=(const message_type& other) {
    value += other.value;
    pr = std::max(pr, other.pr);
    return *this;
  }
};


const size_t NCPUS = 4;
const size_t NUM_VERTICES = 101;
std::vector<atomic<int> > correctness_counter;


template <typename SchedulerType>
void test_scheduler_basic_functionality_single_threaded() {
  graphlab_options opts;
  opts.set_ncpus(NCPUS);
  SchedulerType sched(NUM_VERTICES, opts);
  const size_t target_value = 100;
  
  // inject a sequence of messages which will sum to 100 per vertex
  for (size_t c = 0;c < target_value; ++c) {
    for (size_t i = 0; i < NUM_VERTICES; ++i) {
      sched.schedule(i, message_type(1, 1.0));
    }
  }
  correctness_counter.clear();
  correctness_counter.resize(NUM_VERTICES, atomic<int>(0));
  sched.start();
  
  // pull stuff out
  bool allcpus_done = false; 
  while(!allcpus_done) {
    allcpus_done = true;
    for (size_t i = 0; i < NCPUS; ++i) {
      vertex_id_type v; message_type m;
      sched_status::status_enum ret = sched.get_next(i, v, m);
      if (ret == sched_status::NEW_TASK) {
        allcpus_done = false;
        correctness_counter[v].inc(m.value);
      }
    }
  }

  // check the counters
  for(size_t i = 0; i < NUM_VERTICES; ++i) {
    TS_ASSERT_EQUALS(correctness_counter[i].value, target_value);
  }
}






template <typename SchedulerType>
void test_scheduler_basic_functionality_parallel() {
  graphlab_options opts;
  opts.set_ncpus(NCPUS);
  SchedulerType sched(NUM_VERTICES, opts);
  const size_t target_value = 100;

  // inject a sequence of messages which will sum to 100 per vertex
  for (size_t c = 0;c < target_value; ++c) {
    for (size_t i = 0; i < NUM_VERTICES; ++i) {
      sched.schedule(i, message_type(1, 1.0));
    }
  }
  correctness_counter.clear();
  correctness_counter.resize(NUM_VERTICES, atomic<int>(0));
  sched.start();

  // pull stuff out
  bool allcpus_done = false;
  while(!allcpus_done) {
    allcpus_done = true;
    for (size_t i = 0; i < NCPUS; ++i) {
      vertex_id_type v; message_type m;
      sched_status::status_enum ret = sched.get_next(i, v, m);
      if (ret == sched_status::NEW_TASK) {
        allcpus_done = false;
        correctness_counter[v].inc(m.value);
      }
    }
  }

  // check the counters
  for(size_t i = 0; i < NUM_VERTICES; ++i) {
    TS_ASSERT_EQUALS(correctness_counter[i].value, target_value);
  }
}


class SerializeTestSuite : public CxxTest::TestSuite {
public:

  void test_scheduler_basic_single_threaded() {
    test_scheduler_basic_functionality_single_threaded<sweep_scheduler<message_type> >();
    test_scheduler_basic_functionality_single_threaded<fifo_scheduler<message_type> >();
    test_scheduler_basic_functionality_single_threaded<priority_scheduler<message_type> >();
    test_scheduler_basic_functionality_single_threaded<queued_fifo_scheduler<message_type> >();
  }

};

